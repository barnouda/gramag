/*
 *  primitive.c
 *  
 *  Created by olivier coutant on 05/07/10.
 *
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "declarations.h"

//=========================================================================================
//=========================================================================================
/*
 * given depth za, determine a depth zb such that in interval [za , zb]
 * volumic mass is varying as: rho=rho0+z*alpha
 *
 * since rho=rho(za) + (z-za)/(zb-za)*(rho(zb)-rho(za))
 *
 *       rho=rho(za)-za/(zb-za)*(rho(zb)-rho(za)) + z/(zb-za)*(rho(zb)-rho(za))
 *
 * then
 *
 *    rho0 = rho(za)-za/(zb-za)*(rho(zb)-rho(za))
 *
 *    alpha = (rho(zb)-rho(za))/(zb-za)
 *
 * When rho is constant in the interval, alpha=0
 */

double getZB(double za, double *rho0, double *alpha, double xs, double ys, int *outGrid)
{
	extern VMass *_vmass;
	VMass  *v;
	double zb,zup,zdown,rhoza,rhozb;
	
	v=_vmass;
	
	//
	// 1st case we are outside the grid
	//
	if (za<v->zmin || za>=v->zmax || xs<v->xmin || xs>v->xmax || ys<v->ymin || ys>v->ymax)
	{
		getLayerFromDepth(za,&zup,&zdown,v);
		zb=zdown;
		*rho0=getVmassFromDepth((za+zb)/2,v);
		*alpha=0.;
		*outGrid=1;

                // 2nd (sub)case we are however above the grid
                if (za<v->zmin && xs>=v->xmin && xs<=v->xmax && ys>=v->ymin && ys<=v->ymax)
                	zb=v->zmin;
	}
	//
	// 2nd case, we are inside the grid, interpolate
	//
 	else
	{int i;
		
		for (i=0;i<v->nz;i++)
		{
			if (za<zvmpos(i,v))
				break;
		}
		zb=zvmpos(i,v);
		rhoza=interpolateVMass0(v,xs,ys,za);
		rhozb=interpolateVMass0(v,xs,ys,zb);	
		*rho0=rhoza-za/(zb-za)*(rhozb-rhoza);
		*alpha=(rhozb-rhoza)/(zb-za);
		*outGrid=0;
	}
	return(zb);
}

/*
 * Integration analytique d'un prisme rectangulaire de la profondeur za a zb, 
 * rho peut varier lineairement avec la profondeur Z
 * mais est assumé constant selon X,Y
 *
 * input:
 *	za: top of the prism. Since Z is positive downward, this is the lowest
 *          z coordinate. Above, we have already process the topography.
 *	xr,yr,zr: receiver coordinates
 *	x1,x2,y1,y2: prism coordinates
 *			(x1,y2) *	* (x2,y2)
 *
 *			(x1,y1) *	* (x2,y1)
 *	res: contains the integration at topography
 *
 * output:
 *	res: updated with integration
 */
void 
prismIntegration(double *za, double *xr, double *yr, double *zr, \
                  double *x1, double *x2, double *y1, double *y2, double *res)
{
	int    outGrid,firstSegment;
	double rho0,alpha,ZM,Zm,YM,Ym,XM,Xm,Gc,Gl,za1,zb1,dS;
	double xs,ys,dh,GlUp,GlDown,GcUp,GcDown;
	extern VMass   *_vmass;
	VMass *v;
	
	v=_vmass;
	
	//
	// Integration over Z is performed on successive segments where rho (volumic mass)
	// varies like rho=a*z+b (a may be 0, but not b) 
	//
	dS=(*x2-*x1)*(*y2-*y1);
	xs=(*x1 + *x2)/2.,
	ys=(*y1 + *y2)/2.,

	// initialization
	za1=(*za);
	XM=(*xr)-(*x1);
	Xm=(*xr)-(*x2); 
	YM=(*yr)-(*y1);
	Ym=(*yr)-(*y2);	
	dh=((*xr)-(xs))*((*xr)-(xs))+((*yr)-(ys))*((*yr)-(ys));

	GlUp=0.;
	GcUp=0.;
	firstSegment=1;
	do                  	// as long as top of prism has not reach max depth
	{
		// given top depth, yield bottom depth for 
		// the prism volume where rho is constant or linear with z
		zb1=getZB(za1,&rho0,&alpha,xs,ys,&outGrid);
		ZM=(*zr)-zb1;	//upper limit of integration (but bottom of segment in terme of altitude)
		Zm=(*zr)-za1;	//lower limit of integration (but top of segment in term of altitude)

		// compute contribution for first integration segment
		if (firstSegment)
		{
			GcUp=(primJ1(XM,Zm,YM)-primJ1(Xm,Zm,YM)) - (primJ1(XM,Zm,Ym)-primJ1(Xm,Zm,Ym)); // contribution du bas de l'élément au dessus
			firstSegment=0;
		}

		if (outGrid)	// density is constant in an element
		{
//			Gc=(primJ1(XM,Zm,YM)-primJ1(Xm,Zm,YM)) - (primJ1(XM,Zm,Ym)-primJ1(Xm,Zm,Ym)) - \
//			   (primJ1(XM,ZM,YM)-primJ1(Xm,ZM,YM)) + (primJ1(XM,ZM,Ym)-primJ1(Xm,ZM,Ym));
//			Gc=-Gc*(rho0+alpha*(*zir));
			GcDown=GcUp;
  			GcUp=(primJ1(XM,ZM,YM)-primJ1(Xm,ZM,YM))-(primJ1(XM,ZM,Ym)-primJ1(Xm,ZM,Ym));

			Gc=(GcUp-GcDown)*(rho0+alpha*(*zr));
			Gl=0.;

		}
		else
		{
			// partie rho constante
//			Gc=(primJ1(XM,Zm,YM)-primJ1(Xm,Zm,YM)) - (primJ1(XM,Zm,Ym)-primJ1(Xm,Zm,Ym)) - \
//			   (primJ1(XM,ZM,YM)-primJ1(Xm,ZM,YM)) + (primJ1(XM,ZM,Ym)-primJ1(Xm,ZM,Ym));	
//			Gc=-Gc*(rho0+alpha*(*zir));
			GcDown=GcUp;
  			GcUp=(primJ1(XM,ZM,YM)-primJ1(Xm,ZM,YM))-(primJ1(XM,ZM,Ym)-primJ1(Xm,ZM,Ym));
			Gc=(GcUp-GcDown)*(rho0+alpha*(*zr));
			// partie rho lineaire
			if (alpha!=0.)
			{ 
//				Gl = ( I1(Xm,XM,Ym,YM,ZM) - I1(Xm,XM,Ym,YM,Zm) - \
//				       I2(Xm,XM,Ym,YM,ZM) + I2(Xm,XM,Ym,YM,Zm) );
				if (GlUp!=0.)
					GlDown=GlUp;
				else
					GlDown=I1(Xm,XM,Ym,YM,Zm)-I2(Xm,XM,Ym,YM,Zm);
				GlUp= I1(Xm,XM,Ym,YM,ZM)-I2(Xm,XM,Ym,YM,ZM);
				Gl = alpha *(GlUp-GlDown);
			}
			else 
			{
				GlUp=0.;
				GlDown=0.;
				Gl=0.;
			}

		}
//fprintf(stderr,"** %f %f %f %f\n",zb1,za1,(Gc+Gl),*res);
		*res+=(Gc+Gl);
		za1=zb1;
	} 
	while (za1<v->zzmax);	// as long as top of prism has not reach max depth



}
