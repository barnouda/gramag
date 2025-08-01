#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include "declarations.h"
#include "mag.h"	//mag
#include "cubature.h"	//mag
/*
 * integration sur un prisme de base parallepipedique (x1,x2,y1,y2)
 * et de hauteur z11, z12,z21,z22
 * de la fonction 
 * 1/Sqrt[(x - x0)^2 + (y - y0)^2 + (f(x,y) - z0)^2]
 * ou f(x,y) est l'interpolation bilineaire de z
 * entre les 4 points hauts du prisme.
 *
 *             (x1,y2,z12) *	* (x2,y2,z22)
 *
 *
 *	       (x1,y1,z11) *    * (x2,y1,z12)
 *
 * (x0,y0,z0)= point d'evaluation du champ de gravite.
 * Si (x0,y0) est dans le prisme, l'integral presente une 
 * fausse singularite levee par chgt variable si z0 est a la surface.
 *
 *
 *  topprismlowres_
 *	->topPrism
 *		->quadSing
 *		->quad
 *			->func
 *  topprismhighres_
 *	->topPrism
 *		->quadSing
 *		->quad
 *			->func
 */

/*
 *
 * version Low Res: on integre sur le prisme en interpolant
 * la topographie au 4 coins du prisme
 *
 * input:
 * xs,ys: centre du prisme
 * dx2,dy2: 1/2 largeur du prisme
 * x0,y0,z0: position recepteur
 * output:
 * G1,G2: resultat integration
 * za: base de l'integration
 *
 * geometrie:
 *             (x1,y2,z12) *    * (x2,y2,z22)\n *\n *\n");
 *             (x1,y1,z11) *    * (x2,y1,z12)\n\n\n");
 */
int
topPrismLowres(double *x1,double *x2,double*y1,double*y2,
				double *x0,double *y0,double *z0,
				double *G1,double *G2, double *za)
{
	int ier=0,id[8],i;
	double z11,z12,z21,z22,X,Y,Z,rho,GUp,GDown,ztop;
	double Xm,XM,Ym,YM,Zm,alpha[8];
	extern  Topo *_topo;
	extern  VMass *_vmass;

	extern DMagnet *_dmagnet;	//mag
	DMagnet *u;			//mag
	u = _dmagnet;			//mag
	
	// variables for hcubature				//mag
	double xmin[2], xmax[2], valx, valy, valz;		//mag
	Argu argu;						//mag

	// I) Contribution sommet segment (topographie)
	// profondeurs aux coins du prisme
	z11=elev(*x1,*y1,_topo);
	z12=elev(*x1,*y2,_topo);
	z21=elev(*x2,*y1,_topo);
	z22=elev(*x2,*y2,_topo);
	
	// top of prism
	ztop=z11;
	ztop=(z12<ztop)? z12: ztop;
	ztop=(z21<ztop)? z21: ztop;
	ztop=(z22<ztop)? z22: ztop;

	// base du prisme pour cette integration
	*za=z11;
	*za=(z12>*za)? z12: *za;
	*za=(z21>*za)? z21: *za;
	*za=(z22>*za)? z22: *za;

	// topo is flat, do not integrate
	if (ztop==*za)
		return(0);

	// densite moyenne
	X=((*x1)+(*x2))/2;	//BUG
	Y=((*y1)+(*y1))/2;
	Z=(z11+z12+z21+z22)/8+(*za)/2;
	rho=interpolateVMass(_vmass,X,Y,Z,id,alpha);

	argu.ux = u->ux;				//mag
	argu.uy = u->uy;				//mag
	argu.uz = u->uz;				//mag
	argu.zr = *z0;					//mag
	xmin[0] = *x0-*x1;				//mag
	xmin[1] = *y0-*y1;				//mag
	xmax[0] = *x0-*x2;				//mag
	xmax[1] = *y0-*y2;				//mag

	argu.za = *z0-ztop;				//mag
	argu.zb = *z0-*za;				//mag

	// III) total contribution
	if(id[0]<0)	// out of grid, fill G1
	{
		argu.J0 = rho;									//mag
		argu.beta = 0;									//mag

		integrate3Comp(xmin[0],xmax[0],xmin[1],xmax[1],argu,&valx,&valy,&valz);
		*G1 += project3Comp(valx,valy,valz);
	}
	else			// inside grid, fill G2 for each node
	{
		for (i=0;i<8;i++)
		{		
			argu.J0 = alpha[i];								//mag
			argu.beta = 0;									//mag
		
			integrate3Comp(xmin[0],xmax[0],xmin[1],xmax[1],argu,&valx,&valy,&valz);
			#pragma omp atomic
			G2[id[i]] += project3Comp(valx,valy,valz);
		}
	}
	
	return(ier);
}
