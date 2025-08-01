/*
 *  primitive.c
 *  
 *  Created by olivier coutant on 05/07/10.
 *
 */

#include <math.h>
#include <stdio.h>
#include "declarations.h"
#include "mag.h"	//mag
#include "cubature.h"	//mag


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
		//ABif (za<v->zmin && xs>=v->xmin && xs<=v->xmax && ys>=v->ymin && ys<=v->ymax)
		//AB	zb=v->zmin;
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
		//if((zb-za)<10)	// BUG !
		//{
		//	printf("za %f zb %f\n",za,zb);
		//}
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

//
// Prism integration
// 	top of the prism = mean of the topography at the four corners
// 	magnetization = J = J0+beta*z <=> RHO0+alpha*z
// 	returns the three components (x,y,z) of the anomaly produced by the prism on (xr,yr,zr)
//

void prismIntegration(double *za, double *xr, double *yr, double *zr, \
					   double *x1, double *x2, double *y1, double *y2, \
					   double *G1, double *G2)
{
	int    outGrid,id[8],i,im,firstSegment;
	double rho0,alpha,ZM,Zm,YM,Ym,XM,Xm,Gc,Gl,za1,zb1,dS;
	double RHO0[8],ALPHA[8],dh,GlUp,GlDown,GcUp,GcDown;
	double xs,ys;
	extern VMass   *_vmass;
	VMass *v;
	v=_vmass;
	extern DMagnet *_dmagnet;	//mag
	DMagnet *u;			//mag
	u = _dmagnet;			//mag
	
	// variables for hcubature
	double xmin[2], xmax[2], valx, valy, valz;
	Argu argu;
	
	//
	// Integration over Z is performed on successive segments where rho (volumic mass)
	// varies like rho=a*z+b (a may be 0, but not b) 
	//
	xs=(*x1 + *x2)/2.,
	ys=(*y1 + *y2)/2.,

	// initialization
	za1=*za;

	argu.ux = u->ux;
	argu.uy = u->uy;
	argu.uz = u->uz;
	argu.zr = *zr;
	xmin[0] = *xr-*x1;
	xmin[1] = *yr-*y1;
	xmax[0] = *xr-*x2;
	xmax[1] = *yr-*y2;

	while(za1<v->zzmax) // as long as top prism has not reach max depth
	{
		// define bottom depth for first segment
		zb1=getZB(za1,&rho0,&alpha,xs,ys,&outGrid,id,RHO0,ALPHA);

		argu.za = *zr-za1;
		argu.zb = *zr-zb1;

		if(outGrid)
		{
			argu.J0 = rho0;
			argu.beta = alpha;

			integrate3Comp(xmin[0],xmax[0],xmin[1],xmax[1],argu,&valx,&valy,&valz);
			*G1 += project3Comp(valx,valy,valz);
		}
		else
		{
			for(i=0;i<8;i++)
			{
				im=id[i];
				
				argu.J0 = RHO0[i];
				argu.beta = ALPHA[i];

				integrate3Comp(xmin[0],xmax[0],xmin[1],xmax[1],argu,&valx,&valy,&valz);

				#pragma omp atomic
				G2[im] += project3Comp(valx,valy,valz);
			}
		}

		za1=zb1;
	}
}
