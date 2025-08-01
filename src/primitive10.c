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
#include "mag.h"		//mag
#include "cubature.h"		//mag

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

//
// Prism integration
// 	top of the prism = mean of the topography at the four corners
// 	magnetization = J = J0+beta*z <=> RHO0+alpha*z
// 	returns the three components (x,y,z) of the anomaly produced by the prism on (xr,yr,zr)
//
void prismIntegration(double *za, double *xr, double *yr, double *zr, \
			double *x1, double *x2, double *y1, double *y2, \
			double *resx, double *resy, double *resz)	//mag
{
	int outGrid;
	double rho0, alpha, za1, zb1;
	double xs, ys;
	extern VMass *_vmass;
	VMass  *v;
	v = _vmass;
	extern DMagnet *_dmagnet;
	DMagnet *u;
	u = _dmagnet;
	
	// variables for hcubature
	double xmin[2], xmax[2];
	double valx,valy,valz;
	Argu argu;

	// initialisation
	xs = (*x1+*x2)/2.;
	ys = (*y1+*y2)/2.;
	za1 = *za;

	argu.ux = u->ux;
	argu.uy = u->uy;
	argu.uz = u->uz;
	argu.zr = *zr;
	xmin[0] = *xr-*x1;
	xmin[1] = *yr-*y1;
	xmax[0] = *xr-*x2;
	xmax[1] = *yr-*y2;

	while(za1<v->zzmax)	// as long as top of prism has not reach max depth
	{
		// given top depth, yield bottom depth for
		// the prism volume where rho is constant or linear with z
		zb1 = getZB(za1,&rho0,&alpha,xs,ys,&outGrid);

		argu.J0 = rho0;
		argu.beta = alpha;
		argu.za = *zr-za1;
		argu.zb = *zr-zb1;

		integrate3Comp(xmin[0],xmax[0],xmin[1],xmax[1],argu,&valx,&valy,&valz);
		*resx+=valx;
		*resy+=valy;
		*resz+=valz;

		za1 = zb1;
	}
}

