/*
 *  primitive.c
 *  
 *  Created by olivier coutant on 05/07/10.
 *
 */

#include <math.h>
#include <stdio.h>
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
 * input:
 * 	(x,y)= center of the prism
 * 	vm= structure describing velocity (density) model
 * output:
 *      (z1,z2)= vertical limit of integration
 *      *id, alpha1,alpha2, rho1,rho2= interpolation and density coefficient related to nodes if any
 *      outgrid: logical 1 or 0 if we are (not) in the grid volume
 */
void interpolateVMass2Depth(VMass *vm,double x,double y,double z1, double z2,int *id, \
	                     double *alpha1, double *alpha2,double *rho1, double *rho2);
double getZB(double za, double *rho0, double *alpha, double xs, double ys, int *outGrid,int *id, double *RHO0,double*ALPHA)
{
	extern VMass *_vmass;
	VMass  *v;  
        double  Aa[8],Ab[8];
	double zb,zup,zdown,rhoza,rhozb,x1,x2;
	
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
//		rhoza=interpolateVMass(v,xs,ys,za,id,Aa);
//		rhozb=interpolateVMass(v,xs,ys,zb,id,Ab);
		// get Vmass at locations (xs,ys,za) and (xs,ys,zb)
		// and get linear coefficients (Aa, Ab) that link vmass 8 nodes value
		// to integrated gravity field
		interpolateVMass2Depth(v,xs,ys,za,zb,id,Aa,Ab,&rhoza,&rhozb);

		x2=1./(zb-za);
		x1=za*x2;
		*rho0=rhoza-x1*(rhozb-rhoza);
		*alpha=(rhozb-rhoza)*x2;

		for (i=0;i<8;i++)
		{
			RHO0[i]=Aa[i]-x1*(Ab[i]-Aa[i]);
			ALPHA[i]=(Ab[i]-Aa[i])*x2;
		}
		*outGrid=0;
	}
	return(zb);
}

/*
 * Integration analytique d'un prisme rectangulaire de la profondeur za a zb, 
 * rho peut varier lineairement avec la profondeur Z
 * mais est assumé constant selon X,Y
 *
 * contribution du prisme rectangulaire de hauteur za a zb (za<zb), centré en (xs,ys), de largeur dx*dy
 *
 * zir= profondeur du point d'evaluation
 * res= resultat de l'integration
 * xs,ys: coordonnees du centre de la colonne
 *
 * La condition za>=topographie est assurée avant appel
 *
 * resultat renvoyé dans res, ier garde par compatibilite
 */
void prismIntegration(double *za, double *xr, double *yr, double *zr, \
					   double *x1, double *x2, double *y1, double *y2, \
					   double *G1,double *G2)
{
	int    outGrid,id[8],i,im,firstSegment;
	double rho0,alpha,ZM,Zm,YM,Ym,XM,Xm,Gc,Gl,za1,zb1,dS;
	double RHO0[8],ALPHA[8],dh,GlUp,GlDown,GcUp,GcDown;
	double xs,ys;
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
	do {
		// define bottom depth for first segment
		zb1=getZB(za1,&rho0,&alpha,xs,ys,&outGrid,id,RHO0,ALPHA);
		Zm=(*zr)-za1;
		ZM=(*zr)-zb1;

		// if discontinuous, GcDown=(GcUp#0) * deltaRho	
		if (firstSegment)
		{
			GcUp=(primJ1(XM,Zm,YM)-primJ1(Xm,Zm,YM)) - (primJ1(XM,Zm,Ym)-primJ1(Xm,Zm,Ym));
			firstSegment=0;
		}
		
		if (outGrid)
		{
//			Gc=-((primJ1b(XM,Zm,YM)-primJ1b(Xm,Zm,YM)) - (primJ1b(XM,Zm,Ym)-primJ1b(Xm,Zm,Ym)) - \
//			     (primJ1b(XM,ZM,YM)-primJ1b(Xm,ZM,YM)) + (primJ1b(XM,ZM,Ym)-primJ1b(Xm,ZM,Ym)));
//
//			Gc=Gc*(rho0);		// constant part
//			Gl=0.;			// linear part
			GcDown=GcUp;
  			GcUp=(primJ1(XM,ZM,YM)-primJ1(Xm,ZM,YM))-(primJ1(XM,ZM,Ym)-primJ1(Xm,ZM,Ym));

			Gc=(GcUp-GcDown)*(rho0+alpha*(*zr));
			Gl=0.;

			(*G1)+=(Gc);

		}
		else // inside velocity grid
		{

			// partie constante

//			  Gc=-((primJ1b(XM,Zm,YM)-primJ1b(Xm,Zm,YM)) - (primJ1b(XM,Zm,Ym)-primJ1b(Xm,Zm,Ym)) - \
//			     (primJ1b(XM,ZM,YM)-primJ1b(Xm,ZM,YM)) + (primJ1b(XM,ZM,Ym)-primJ1b(Xm,ZM,Ym)));
			GcDown=GcUp;
  			GcUp=(primJ1(XM,ZM,YM)-primJ1(Xm,ZM,YM))-(primJ1(XM,ZM,Ym)-primJ1(Xm,ZM,Ym));

			//rhoalpha1=0;
			for (i=0;i<8;i++)
			{
				im=id[i];
				#pragma omp atomic
				G2[im] += (GcUp- GcDown)*(RHO0[i]+ALPHA[i]*(*zr));
				//rhoalpha1 += (RHO0[i]+ALPHA[i]*(*zr))*v->vm[im];
			}
			//Gc=(GcUp-GcDown)*rhoalpha1;//ModifOC

//			Gl = ( I1(Xm,XM,Ym,YM,ZM) - I1(Xm,XM,Ym,YM,Zm) - \
//				       I2(Xm,XM,Ym,YM,ZM) + I2(Xm,XM,Ym,YM,Zm) );
			if (GlUp!=0.)
				GlDown=GlUp;
			else
				GlDown=I1(Xm,XM,Ym,YM,Zm)-I2(Xm,XM,Ym,YM,Zm);
			GlUp= I1(Xm,XM,Ym,YM,ZM)-I2(Xm,XM,Ym,YM,ZM);
			Gl=(GlUp-GlDown);				

			for (i=0;i<8;i++)
			{
				im=id[i];
				#pragma omp atomic
				G2[im] += Gl*ALPHA[i];
				//alpha1+=ALPHA[i]*v->vm[im];
			}


			//Gl=Gl*alpha;//ModifOC
		}		
		za1=zb1;
	}
	while (za1<v->zzmax);	// as long as top of prism has not reach max depth
}
int
myisnan_(double *x)
{
if (isnan((*x)))
	return(1);
else
	return(0);
}
