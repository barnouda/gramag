#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include "declarations.h"
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
typedef struct
{
	long double  a1,b1,a2,b2,denx;
	long double  x0,y0,z0;
	long double  x1,x2;
} Arg;

double quad(long double f(),double a,double b,double tol,int *neval,int*warni,Arg*arg);
double quadsing(double y0, double eps, double x0, double x1, double x2, Arg *arg);
double primJ1(double x, double a, double b);
/*
 * Integration for the part of the prism that is close to the topography
 *
 * input:
 *      x1,x2,y1,y2= horizontal locations of the prism corner
 *      z11,z12,z21,z22= vertical locations of the prism corner at topography
 *      x0,y0,z0= receiver location
 *
 * output: 
 *      ier: return error code
 */

double topPrism(double x1, double x2, double y1, double y2, 
                double z11, double z12, double z21, double z22,
	        double x0,double y0, double z0, int*ier)
{
	int    neval,inPrism;
	double deny,res,dz;
	long double func(long double, Arg*);
	long double epsabs,eps;
	Arg	arg;
extern  Topo  *_topo;
extern  int    projection;

// evaluation point is in the prism
	inPrism=(x0>=x1 && x0<=x2 && y0>=y1 && y0<=y2);

// elevation / topography
	if (!projection && inPrism)
		dz=z0-elev(x0,y0,_topo);
	else
		dz=0.;

// width of singular zone
        eps=1.e-4;
        eps=(y2-y1)/1000;

// variables globales pour la primitives
	arg.x0=x0;	
	arg.y0=y0;
	arg.z0=z0;

	arg.x1=x1;
	arg.x2=x2;

	arg.denx=1/(x2-x1);
	deny=1/(y2-y1);
	arg.a1=(z12-z11)*deny;
	arg.b1=(z11*y2-z12*y1)*deny;
	arg.a2=(z22-z21)*deny;
	arg.b2=(z21*y2-z22*y1)*deny;
	epsabs=1.e-6;

// integrate singularity
	if (inPrism && dz==0.)
	{
		eps=1.e-4;
		res=quadsing(y0,eps,x0,x1,x2,&arg);
        	res+=quad(func,y1,y0-eps,epsabs,&neval,ier,&arg);
        	res+=quad(func,y0+eps,y2,epsabs,&neval,ier,&arg);
		if (*ier>0 ||!isfinite(res))
			fprintf(stderr," recept in prism ier=%d %f\n",*ier,res);
	}
	else
	{
// integration numerique sur y 
		epsabs=1.e-6;
        	res=quad(func,y1,y2,epsabs,&neval,ier,&arg);
		if (*ier>0 ||!isfinite(res))
			fprintf(stderr,"recept hors prism ier=%d %f\n",*ier,res);
	}

	return(res); 
}
/*
 * compute the contribution of a small prism that contains the observation point
 * The singularity of func(à is removed by changing to cylindrical coordinates
 * input:
 *      x0,y0: observation point location (z0=a*x0+b)
 *      eps: width of the prism along y direction
 *      x1,x2: dimesnion along x direction
 *      Arg: temporary values
 *
 * !!!!! CAUTION !!!! the observation point must not be on the side 
 */
double quadsing(double y0, double eps, double x0, double x1, double x2, Arg *arg)
{
double  a,z1,z2,c,X1,X2,I,theta1M,theta2M;

        // compute slope a so that f(x,y)=z=a*x+b when y=y0
        z1=arg->a1*y0+arg->b1;
        z2=arg->a2*y0+arg->b2;
        a=(z2-z1)*arg->denx;

        //
        c=sqrtl(a*a+1);
        X1=c*(x0-x1);
        X2=c*(x2-x0);

// if X1 or X2 is zero, the evaluation point is on the side. Contribution is 0 from null surface
        theta1M=(X1==0.)? 0.: atan2(eps,X1);
        theta2M=(X2==0.)? 0.: atan2(eps,X2);

        I=X2*(log(fabsl(tan(theta2M/2+M_PI/4))))+X1*(log(fabsl(tan(theta1M/2.+M_PI/4.))));
        I=I*4/c;

        return(I);
}

/*
 * primitive de Integrate[1/Sqrt[(x - x0)^2 + (y - y0)^2 + (a*x + b - z0)^2], x] 
 * integre par rapport a x, a y fixe
 */
long double Primx(long double x,long double y,Arg *arg)
{
	long double a,a2,b,b2,z1,z2,X,Y,c,v,xx,yy;

	z1=arg->a1*y+arg->b1;
	z2=arg->a2*y+arg->b2;
	a=(z2-z1)*arg->denx;
	a2=a*a;
	b=(arg->x2*z1-arg->x1*z2)*arg->denx;
	b2=b*b;
	c=sqrtl((1+a2));
	Y = y - arg->y0;
//	Y=1.e-10;
//fprintf(stderr,"%.20f %.20f %.20f %.20f\n",Y,y,arg->y0,1.e-10);
	X=x-arg->x0;
	xx=b2+a2*x*x+Y*Y+X*X+2*b*(a*x-arg->z0)+arg->z0*arg->z0-2*x*a*arg->z0;
	yy=a*b+x+a2*x-arg->x0-a*arg->z0+c*sqrtl(fabsl(xx));
        v=logl(fabsl(yy))/c;

//fprintf(stderr,"%Lf %Le %Lf (%Le %Lg)\n",xx,yy,v,Y,X);
//fprintf(stderr,"%Lg %Lg %Lg %Lg %Lg %Lg %Lf //\n",a*b,x,a2*x,arg->x0,a*arg->z0,c*sqrt(fabs(xx)),a*b+x+a2*x-arg->x0-a*arg->z0+c*sqrt(fabs(xx)));
//      v=log(a*b+x+a2*x-_x0-a*_z0+c*sqrt(fabs(b2+(1+a2)*x*x+_x0*_x0+y*y-2*y*_y0+_y0*_y0+2*b*(a*x-_z0)+_z0*_z0-2*x*(_x0+a*_z0))))/c;
	return(v);
}
/*
 * fonction de y a integrer par quadrature
 */
long double func(long double y, Arg *arg)
{
	long double v;
	v=Primx(arg->x2,y,arg)-Primx(arg->x1,y,arg);
	return (v);
}

/*
 * interface fortran
 *
 * version High Res: on integre sur le prisme en decoupant
 * sur des sous-prismes qui suivent la discretisation de la topo
 *
 * input:
 *
 * px1,px2,py1,py2: coordonn�es x,y de la base du prisme
 * x0,y0,z0: position du recepteur
 *
 * output:
 * G1,G2: resultat integration
 * za: base de l'integration
 *
 * geometrie:
 *             (x1,y2,z12) *    * (x2,y2,z22)\n *\n *\n");
 *             (x1,y1,z11) *    * (x2,y1,z12)\n\n\n");
 */

int
topPrismHighres(double *px1,double *px2,double*py1,double*py2,
				 double *x0,double *y0,double *z0,
				 double *G1,double *G2,double *pza)
{
	int    i,j,i1,i2,j1,j2,nx,ny,ier,id[8];
	double z11,z12,z21,z22,x1,x2,y1,y2,rho,za,alpha[8];
	double X,Y,Z,GUp,GDown;
	double Xm,XM,Ym,YM,Zm;
	double *xx1,*xx2,*yy1,*yy2,*zz11,*zz12,*zz21,*zz22;
extern  Topo  *_topo;
	Topo  *t;
extern  VMass *_vmass;
	t=_topo;

	// get index for topo cell containing (x1,y1)
	ielev(*px1,*py1,t,&i1,&j1);
	// get index for topo cell containing (x2,y2)
	ielev(*px2,*py2,t,&i2,&j2);

	// check wether cell bounds match topo nodes
	nx=i2-i1+2;
	ny=j2-j1+2;
	if (*px1==t->xmin+i1*t->dx)
	{
		nx--;
		i1++;
	}
	if (*py1==t->ymin+j1*t->dy)
	{
		ny--;
		j1++;
	}
	if (*px2==t->xmin+i2*t->dx)
	{
		nx--;
		i2--;
	}
	if (*py2==t->ymin+j2*t->dy)
	{
		ny--;
		j2--;
	}
	// create arrays for integration cells	
	// each integration cell matches the topo cell
	xx1=(double*)malloc(sizeof(double)*nx);
	xx2=(double*)malloc(sizeof(double)*nx);
	yy1=(double*)malloc(sizeof(double)*ny);
	yy2=(double*)malloc(sizeof(double)*ny);
	zz11=(double*)malloc(sizeof(double)*nx*ny);
	zz12=(double*)malloc(sizeof(double)*nx*ny);
	zz21=(double*)malloc(sizeof(double)*nx*ny);
	zz22=(double*)malloc(sizeof(double)*nx*ny);

	xx1[0]=*px1;
	yy1[0]=*py1;
	xx2[nx-1]=*px2;
	yy2[ny-1]=*py2;
	for(i=1;i<nx;i++)
		xx2[i-1]=xx1[i]=t->xmin+(i-1+i1)*t->dx;
	for (j=1;j<ny;j++)
		yy2[j-1]=yy1[j]=t->ymin+(j-1+j1)*t->dy;

	// compute elevation at nodes	
	za=elev(xx1[0],yy1[0],t);
	for (i=0;i<nx;i++)
	{
		x1=xx1[i];
		x2=xx2[i];
		for (j=0;j<ny;j++)
		{
			y1=yy1[j];
			y2=yy2[j];
			z11=zz11[j+i*ny]=elev(x1,y1,t);
			z12=zz12[j+i*ny]=elev(x1,y2,t);
			z21=zz21[j+i*ny]=elev(x2,y1,t);
			z22=zz22[j+i*ny]=elev(x2,y2,t);
			// base du prisme pour cette integration
			za=(z11>za)? z11: za;
			za=(z12>za)? z12: za;
			za=(z21>za)? z21: za;
			za=(z22>za)? z22: za;
		}
	}
	//
	// estimate average density at center of mass
	//
	X=((*px1)+(*px2))/2.;
	Y=((*py1)+(*py2))/2.;
	Z=(elev(X,Y,t)+za)/2;
	rho=interpolateVMass(_vmass,X,Y,Z,id,alpha);


	// integration
	GDown=0;
	for (i=0;i<nx;i++)
	{
		x1=xx1[i];
		x2=xx2[i];
		for (j=0;j<ny;j++)
		{
			y1=yy1[j];
			y2=yy2[j];
			z11=zz11[j+i*ny];
			z12=zz12[j+i*ny];
			z21=zz21[j+i*ny];
			z22=zz22[j+i*ny];
			GDown += topPrism(x1,x2,y1,y2,z11,z12,z21,z22,*x0,*y0,*z0,&ier);
			if (ier>0 || !isfinite((GDown)))
			{FILE *fp;
				fprintf(stderr,"topprismhighres ier=%d %f\n",ier,GDown);
#define DEBUG 1
#ifdef DEBUG
#define wlog(x) (fwrite(&(x),sizeof(double),1,fp))
				fprintf(stderr,"%lf %lf %lf %lf\n",x1,x2,y1,y2);
				fprintf(stderr,"%lf %lf %lf %lf\n",z11,z12,z21,z22);
				fprintf(stderr,"%lf %lf %lf\n",*x0,*y0,*z0);
				fp=fopen("quad.arg","w");
				wlog(x1);wlog(x2);wlog(y1);wlog(y2);
				wlog(z11);wlog(z12);wlog(z21);wlog(z22);
				wlog(*x0);wlog(*y0);wlog(*z0);
				fclose(fp);
#endif
				if (!isfinite((GDown)))
				{	
					return(3);
				}
			}
		}
	}
	free(xx1);
	free(xx2);
	free(yy1);
	free(yy2);
	free(zz11);
	free(zz12);
	free(zz21);
	free(zz22);
	
	// II) contribution from bottom (in term of altitude)
	XM=(*x0)-(*px1);
	Xm=(*x0)-(*px2); 
	YM=(*y0)-(*py1);
	Ym=(*y0)-(*py2);	
	Zm=(*z0)-za;	
	GUp=((primJ1(XM,Zm,YM)-primJ1(Xm,Zm,YM)) - (primJ1(XM,Zm,Ym)-primJ1(Xm,Zm,Ym)));
	
	// III) total contribution
	if (id[0]<0)	// out of grid, fill G1
		*G1+=(GUp-GDown)*rho;
	else			// inside grid, fill G2 for each node
	{
		for (i=0;i<8;i++)
			#pragma omp atomic
			G2[id[i]]+=(GUp-GDown)*alpha[i];
	}

	*pza=za;
	return(ier);
}

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


	GDown=topPrism(*x1,*x2,*y1,*y2,z11,z12,z21,z22,*x0,*y0,*z0,&ier);

	if (ier>0)
	{FILE*fp;
		fprintf(stderr,"topprismlowres ier=%d\n",ier);
#define DEBUG 1
#ifdef DEBUG
#define wlog(x) (fwrite(&(x),sizeof(double),1,fp))
				fprintf(stderr,"%le %le %le %le\n",*x1,*x2,*y1,*y2);
				fprintf(stderr,"%le %le %le %le\n",z11,z12,z21,z22);
				fprintf(stderr,"%le %le %le\n",*x0,*y0,*z0);
				fp=fopen("quad.arg","w");
				wlog(x1);wlog(x2);wlog(y1);wlog(y2);
				wlog(z11);wlog(z12);wlog(z21);wlog(z22);
				wlog(*x0);wlog(*y0);wlog(*z0);
				fclose(fp);
#endif 
		if (ier==3)
			return(3);
	}
	// II) contribution bas segment
	XM=(*x0)-(*x1);
	Xm=(*x0)-(*x2); 
	YM=(*y0)-(*y1);
	Ym=(*y0)-(*y2);	
	Zm=(*z0)-(*za);	

	GUp=((primJ1(XM,Zm,YM)-primJ1(Xm,Zm,YM)) - (primJ1(XM,Zm,Ym)-primJ1(Xm,Zm,Ym)));

	// III) total contribution
	if (id[0]<0)	// out of grid, fill G1
		*G1+=(GUp-GDown)*rho;
	else			// inside grid, fill G2 for each node
	{
		for (i=0;i<8;i++)
			#pragma omp atomic
			G2[id[i]]+=(GUp-GDown)*alpha[i];
	}
	
	return(ier);
}
/*
 * integration par methode simpson adaptative
 * tiree de Matlab quad.m
 */
/*
function [Q,fcnt] = quad(funfcn,a,b,tol,trace,varargin)
%QUAD   Numerically evaluate integral, adaptive Simpson quadrature.
%   Q = QUAD(FUN,A,B) tries to approximate the integral of scalar-valued
%   function FUN from A to B to within an error of 1.e-6 using recursive
%   adaptive Simpson quadrature. FUN is a function handle. The function
%   Y=FUN(X) should accept a vector argument X and return a vector result
%   Y, the integrand evaluated at each element of X.
%
%   Q = QUAD(FUN,A,B,TOL) uses an absolute error tolerance of TOL
%   instead of the default, which is 1.e-6.  Larger values of TOL
%   result in fewer function evaluations and faster computation,
%   but less accurate results.  The QUAD function in MATLAB 5.3 used
%   a less reliable algorithm and a default tolerance of 1.e-3.
%
%   Q = QUAD(FUN,A,B,TOL,TRACE) with non-zero TRACE shows the values
%   of [fcnt a b-a Q] during the recursion. Use [] as a placeholder to
%   obtain the default value of TOL.
%
%   [Q,FCNT] = QUAD(...) returns the number of function evaluations.
%
%   Use array operators .*, ./ and .^ in the definition of FUN
%   so that it can be evaluated with a vector argument.
%
%   Notes:
%   QUAD may be most efficient for low accuracies with nonsmooth
%   integrands.
%   QUADL may be more efficient than QUAD at higher accuracies
%   with smooth integrands.
%   QUADGK may be most efficient for oscillatory integrands and any smooth
%   integrand at high accuracies. It supports infinite intervals and can
%   handle moderate singularities at the endpoints. It also supports
%   contour integration along piecewise linear paths.
%   QUADV vectorizes QUAD for array-valued FUN.
%
%   Example:
%      Q = quad(@myfun,0,2);
%   where the file myfun.m defines the function:
%      %-------------------%
%      function y = myfun(x)
%      y = 1./(x.^3-2*x-5);
%      %-------------------%
%
%   or, use a parameter for the constant:
%      Q = quad(@(x)myfun2(x,5),0,2);
%   where the file myfun2.m defines the function:
%      %----------------------%
%      function y = myfun2(x,c)
%      y = 1./(x.^3-2*x-c);
%      %----------------------%
%
%   Class support for inputs A, B, and the output of FUN:
%      float: double, single
%
%   See also QUADV, QUADL, QUADGK, QUAD2D, DBLQUAD, TRIPLEQUAD, TRAPZ, FUNCTION_HANDLE.

%   Based on "adaptsim" by Walter Gander.
%   http://www.inf.ethz.ch/personal/gander
%
%   Reference:
%   [1] W. Gander and W. Gautschi, Adaptive Quadrature - Revisited,
%       BIT Vol. 40, No. 1, March 2000, pp. 84-101.
%
%   Copyright 1984-2009 The MathWorks, Inc.
%   $Revision: 5.26.4.11 $  $Date: 2010/02/01 03:13:15 $
*/

long double quadstep(long double f(long double,Arg*),long double a,long double b,long double fa,
                long double fc,long double fb,long double tol,
                int *neval,long double hmin,int *warn, Arg*arg);
double quad(long double f(long double,Arg*),double a,double b,double tol,int *neval,int *warn, Arg*arg)
{
	int w,i;
	long double h,x[7],y[7],hmin,Q1,Q2,Q3,Q;
	double Qd;

*warn=0;
//Initialize with three unequal subintervals.
h = 0.13579*(b-a);
x[0] = a;
x[1] = a+h;
x[2] = a+2*h;
x[3] = (a+b)/2;
x[4] = b-2*h;
x[5] = b-h;
x[6] = b;
for (i=0;i<7;i++)
	y[i]=f(x[i],arg);
(*neval) = 7;

//Call the recursive core integrator.
hmin = __DBL_EPSILON__;
hmin = LDBL_EPSILON;
Q1= quadstep(f,x[0],x[2],y[0],y[1],y[2],tol,neval,hmin,&w,arg);
*warn=(*warn<w)? w: *warn;
Q2= quadstep(f,x[2],x[4],y[2],y[3],y[4],tol,neval,hmin,&w,arg);
*warn=(*warn<w)? w: *warn;
Q3= quadstep(f,x[4],x[6],y[4],y[5],y[6],tol,neval,hmin,&w,arg);
*warn=(*warn<w)? w: *warn;
Q = Q1+Q2+Q3;
Qd=(double)Q;

return(Qd);
}


long double quadstep(long double f(long double,Arg*),long double a,long double b,long double fa,long double fc,long double fb,long double tol,
                int *neval,long double hmin,int *warn,Arg*arg)
{
	int warn1,warn2;
	long double h,c,x1,x2,y1,y2,Q,Q1,Q2,Qac,Qcb;
	long double fd,fe;
 	int maxfcnt = 10000;	

//Evaluate integrand twice in interior of subinterval [a,b].
h = b - a;
c = (a + b)/2;
x1 = (a + c)/2;
x2 = (c + b)/2;
y1 = f(x1,arg);
y2 = f(x2,arg);
*neval = *neval + 2;
fd = y1;
fe = y2;

//Three point Simpson's rule.
Q1 = (h/6)*(fa + 4*fc + fb);

//Five point double Simpson's rule.
Q2 = (h/12)*(fa + 4*fd + 2*fc + 4*fe + fb);

//One step of Romberg extrapolation.
Q = Q2 + (Q2 - Q1)/15;

if (*neval > maxfcnt)
{
    *warn = 2;
    return(Q);
}
if (!isfinite(Q))
{
fprintf(stderr,"=== %Lf %Lf \n",fd,fe);
    *warn=3;
    return(Q);
}
if (fabsl(Q2 - Q) <= tol)
{
//    % Accuracy over this subinterval is acceptable.
    *warn = 0;
    return(Q);
}
if (fabsl(h) < hmin || c == a || c == b)
{
    // Minimum step size reached; singularity possible.
    *warn = 1;
}

//Subdivide into two subintervals.
Qac = quadstep(f,a,c,fa,fd,fc,tol,neval,hmin,&warn1,arg);
Qcb = quadstep(f,c,b,fc,fe,fb,tol,neval,hmin,&warn2,arg);
Q = Qac + Qcb;
*warn = (warn1>warn2)? warn1: warn2;
return(Q);
}
