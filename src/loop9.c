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
				double *G1,double *G2, double *za);
int
topPrismHighres(double *px1,double *px2,double*py1,double*py2,
				 double *x0,double *y0,double *z0,
				 double *G1,double *G2,double *pza);
void prismIntegration(double *za, double *xr, double *yr, double *zr, \
					   double *x1, double *x2, double *y1, double *y2, \
					   double *G1,double *G2);

/*
 *
 * loops to integrate over prisms and receiver
 *
 * X,Y,DXY: prism origin and width
 * nprims:  number of prisms
 * xr,yr,zr: receiver positions
 * nr: number of receiver
 * mxyz: number of nodes in density/velicity grid
 * G1,G2: vector and matrice to compute gravity field as a function of volumic mass in density/velicity grid
 */
void loop(Prism *prism,double *xr,double *yr,double *zr,int nr,int mxyz,double *G1,double *G2)
{

	int     i,ir,ier;
	double  res,za,xs,ys,zs,xir,yir,zir;
	double  x1,x2,y1,y2,dist2;
	double  dx;
// G1 is the part that does not depend on the heterogeneous velocity grid
// G2 is the part that depends on the heterogeneous velocity grid


#pragma omp parallel 	num_threads(NTHREAD) \
			default(none) \
			shared(G1,G2,xr,yr,zr,nr,mxyz,prism) \
			private (i,ir,xs,x1,x2,ys,y1,y2,xir,yir,zir,dist2,ier,res,za,dx)
	{
	double *G2t;
	G2t=malloc(sizeof(double)*mxyz);

#pragma omp for schedule(dynamic,10)
	for (ir=0;ir<nr;ir++)
	{
		printf("receiver %d/%d\n",ir,nr);
		xir=xr[ir];
		yir=yr[ir];
		zir=zr[ir];
		memset(G2t,0,sizeof(double)*mxyz);
		res=0;
        	for (i=0;i<prism->NPrism;i++)
        	{
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
			 ! integration 
			 ! highres: integrate prism using free surface discretization that
			 !          is (and should be) smaller than prism discretization
			 ! lowres: integrate using prism discretization and interpolate free surface
			 */
#ifdef DEBUG
res=0;
memset(G2t,0,sizeof(double)*mxyz);
#endif
			/*
			 * integration on the volume is expressed as the summation of surface
			 * integration at each surface of discontinuity.
			 *
			 * I) The first surface is the topography and is integrated in topPrismLowres
			 *    or topPrismHighres
			 */
			if (dist2>1000.*1000 || noTopo())
				ier=topPrismLowres(&x1,&x2,&y1,&y2,&xir,&yir,&zir,&res,G2t,&za);
			else
				ier=topPrismHighres(&x1,&x2,&y1,&y2,&xir,&yir,&zir,&res,G2t,&za);
			/*
			 * II) The second (and following) surfaces is (are) layers fat boundaries
			 *     ans is (are) integrated in prismIntegration
			 */
			prismIntegration(&za, &xir, &yir, &zir, &x1,&x2,&y1,&y2,&res,G2t);

#ifdef DEBUG
{double gg2=0.;
extern VMass *_vmass;
int jj;
static double resg=0;
			gg2=0;
			for (jj=0;jj<mxyz;jj++)
				gg2+=G2t[jj]*((double)_vmass->vm[jj]);
			resg+=(res+gg2);
			fprintf(stderr,"%d %f %f\n",i,res+gg2,resg);
}
#endif
				
		}// end prism loop
		#pragma omp atomic
		G1[ir]+=res;
		for (i=0;i<mxyz;i++)
			G2[ir*mxyz+i]+=G2t[i];


		
		
	}// end receiver loop
	
	free(G2t);
	}// end parallel region
}

