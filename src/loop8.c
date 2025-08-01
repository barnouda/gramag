/*
 *
 *  loop.c
 *
 *  Created by Olivier Coutant on 04/03/11.
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "declarations.h"

int
topPrismLowres(double *x1,double *x2,double*y1,double*y2,
	       double *x0,double *y0,double *z0,
	       double *res,double *za);
int
topPrismHighres(double *px1,double *px2,double*py1,double*py2,
				 double *x0,double *y0,double *z0,
				 double *res,double *pza);
void 
prismIntegration(double *za, double *xr, double *yr, double *zr, \
                  double *x1, double *x2, double *y1, double *y2, double *res);

/*
 *
 * loops to integrate over prisms and receiver
 *
 * X,Y,DXY: prism origin and width
 * nprims:  number of prisms
 * xr,yr,zr: receiver positions
 * nr: number of receiver
 * G: gravity field at receiver locations
 */
void loop(Prism *prism,double *xr,double *yr,double *zr,int nr,double *G)
{

	int      i,j,ir,ier;
	double   res,za,xs,ys,xir,yir,zir;
	double   x1,x2,y1,y2,dist2,dx;

	fprintf(stderr,"loop\n");

/*
!
! Integration from free surface to depth where density is assumed constant
!  Integration is performed analytically along z, then summed up
!  along x and y direction, assuming constant grid spacing
!
*/
#pragma omp parallel 	num_threads(NTHREAD) \
			default(none) \
            		shared(G,xr,yr,zr,nr,prism,stderr) \
            		private (i,j,ir,xs,x1,x2,ys,y1,y2,xir,yir,zir,dist2,ier,res,za,dx)
	{//parallel zone
#pragma omp for schedule(dynamic,10)
	for (ir=0;ir<nr;ir++)
	{
		printf("receiver %d/%d %d\r",ir,nr,prism->NPrism);
		xir=xr[ir];
		yir=yr[ir];
		zir=zr[ir];
		res=0;
// boucle prism
		for (i=0;i<prism->NPrism;i++)
		{double ress;
			if (prism->idxy[i]==DX_REFINED) 
				continue;
			else if (prism->idxy[i]==DX1)
				dx=prism->dx1;
			else
				dx=prism->dx2;

			x1=(double)prism->X[i];
			x2=x1+dx;
			xs=(x1+x2)/2.;
			y1=(double)prism->Y[i];
			y2=y1+dx;
			ys=(y1+y2)/2.;

			dist2=(ys-yir)*(ys-yir)+(xs-xir)*(xs-xir);
			// if receiver is on the side of the prism, shift it
#ifndef MIN
#define MIN(x,y) ((x)<(y)?(x):(y))
#endif
			if (xir==x1)
				xir+=MIN(0.01,((x2-x1)/100.));
			else if (xir==x2)
				xir-=MIN(0.01,((x2-x1)/100.));
			if (yir==y1)
				yir+=MIN(0.01,((y2-y1)/100.));
			else if (yir==y2)
				yir-=MIN(0.01,((y2-y1)/100.));


			/*
			* integration 
			* highres: integrate prism using free surface discretization that
			*          is (and should be) smaller than prism discretization
			* lowres: integrate using prism discretization and interpolate free surface
			*/
			/*
			 * integration on the volume is expressed as the summation of surface
			 * integration at each surface of discontinuity.
			 *
			 * I) The first surface is the topography and is integrated in topPrismLowres
			 *    or topPrismHighres
			 */
			ress=0;
			if (dist2>500.*500 || noTopo())
				ier=topPrismLowres(&x1,&x2,&y1,&y2,&xir,&yir,&zir,&ress,&za);
			else
				ier=topPrismHighres(&x1,&x2,&y1,&y2,&xir,&yir,&zir,&ress,&za);
			if (ier==3)
			{
				fprintf(stderr,"topprism integration error %d %d\n",ier,i);
				exit(-1);
			}
			/*
			 * II) The second (and following) surfaces is (are) layers flat boundaries
			 *     and is (are) integrated in prismIntegration, from za to bottom
			 */
			prismIntegration(&za, &xir, &yir, &zir, &x1,&x2,&y1,&y2,&ress);
			res+=ress;
#ifdef DEBUG
			fprintf(stderr,"%d %f %f\n",i,ress,res);
#endif
		}//end prism loop
		G[ir]+=res;
	}//end receiver loop
	}//end parallel zone
}
