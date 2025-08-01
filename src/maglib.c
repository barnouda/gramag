//
// CONTENT //
//
// void read_dir(DMagnet*dmagnet);
//
// double fxyz(double x, double y, double z, Argu *argu);
// int fxy(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval);
// void integrate3Comp(double xinf, double xsup, double yinf, double ysup, Argu argu, \
		double *valx, double *valy, double *valz);
//
// double project3Comp(double compx, double compy, double compz);
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mag.h"
#include "cubature.h"

//
// Read input file with directions for magnetism
//
// Input file dir.par
// 	first line: 	x, y and z component of magnetization of bodies
// 	second line:	x, y and z component of regional magnetic field
//
// Store components of unit vectors in variable dmagnet
//
// Input coordinates in the usual system for geomagnetism
// 	x northward
// 	y eastward
// 	z downward
// Output coordinates is the system used by the program
// 	x eastward
// 	y northward
// 	z downward
//
void read_dir(DMagnet*dmagnet)
{

	FILE	*fp;

	double u,U;
	double x,y,z;
	double X,Y,Z;
	int r;

	// read vector components from input file
	fp = fopen("dir.par","r");
	if (fp==NULL)
	{
		fprintf(stderr,"could not open file <dir.par>\n");
		exit(-1);
	}
	r = fscanf(fp,"%lf %lf %lf",&x,&y,&z);
	r = fscanf(fp,"%lf %lf %lf",&X,&Y,&Z);

	// compute vector lengths
	u = sqrt(x*x+y*y+z*z);
	U = sqrt(X*X+Y*Y+Z*Z);
	
	// compute components of unit vector of magnetization
	dmagnet->ux = y/u;
	dmagnet->uy = x/u;
	dmagnet->uz = z/u;

	// compute components of unit vector of region field
	dmagnet->Ux = Y/U;
	dmagnet->Uy = X/U;
	dmagnet->Uz = Z/U;
	
#ifdef DEBUGMAG
	printf("Magnetization direction %f %f %f\n",dmagnet->ux,dmagnet->uy,dmagnet->uz);
	printf("Geomagn field direction %f %f %f\n",dmagnet->Ux,dmagnet->Uy,dmagnet->Uz);
#endif
}


// Integrale analytique en z, fonctoin de x et y

// Primitive in z for the three components x, y and z
// We integrate the function f(X,Y,Z) = (J0+beta*z)*(AZ^2+BZ+C)/R^5
// 	X = xr-x
// 	Y = yr-y
// 	Z = zr-z
// 	R = sqrt(X^2+Y^2+Z^2)
// 	A, B and C
// 		are different for the three components x, y and z
//		depends on X and Y (numerical integration)
// The integral in Z takes the form
// 	intf = c0*int(1/R^5) + c1*int(Z/R^5) + c2*int(Z^2/R^5) + c3*int(Z^3/R^5)
// 	     = c0*    p0     + c1*    p1     + c2*      p2     + c3*    p3
double fxyz(double x, double y, double z, Argu *argu)
{
	double ux = argu->ux;
	double uy = argu->uy;
	double uz = argu->uz;
	double J0 = argu->J0;
	double beta = argu->beta;
	double zr = argu->zr;
	double za = argu->za;
	double zb = argu->zb;
	int comp = argu->comp;
	double A, B, C;
	double z2, dh2, r, r3;
	double p0, p1, p2, p3;
	double c0, c1, c2, c3;
	double p;

	if(comp==1)		// x-component
	{
		A = ux;
		B = -3*uz*x;
		C = ux*(-2*x*x+y*y)-3*uy*x*y;
	}
	else if(comp==2)	// y-component
	{
		A = uy;
		B = -3*uz*y;
		C = uy*(x*x-2*y*y)-3*ux*x*y;
	}
	else if(comp==3)	// z-component
	{
		A = -2*uz;
		B = -3*ux*x-3*uy*y;
		C = uz*(x*x+y*y);
	}

	z2 = z*z;
	dh2 = x*x+y*y;
	r = sqrt(dh2+z2);
	r3 = r*r*r;

	if(dh2!=0)
	{
		p0 = (3*dh2*z+2*z2*z)/(3*dh2*dh2*r3);
		p2 = (z2*z)/(3*dh2*r3);
	}
	else
	{
		p0 = -z/(4*r3*r*r);
		p2 = -z/(2*r3*r*r);
	}
	p1 = -1/(3*r3);
	p3 = -(2*dh2+3*z2)/(3*r3);

	c0 = (J0+beta*zr)*C;
	c1 = (J0+beta*zr)*B-beta*C;
	c2 = (J0+beta*zr)*A-beta*B;
	c3 = -beta*A;

#ifdef DEBUGMAG
	//c0 = 0;
	//c1 = 1000000;
	//c2 = 0;
	//c3 = 0;
#endif
	p = c0*p0 + c1*p1 + c2*p2 + c3*p3;
	return p;

}

//
// Function f called by hcubature:
// 	int f(unsigned ndim, const double *x, void *fdata,
//            unsigned fdim, double *fval);
//    	Input
//         	x       array of coordinates where the function has to be evaluated
//         	ndim    number of dimensions (ie length of x)
//         	fdata   for any additional information to compute the function
//    	Output
//         	fval    array of function values
//         	fdim    length of fval
//
int fxy(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval)
{
	Argu *argu;
	argu = ((Argu *) fdata);
	double fxyza = fxyz(x[0], x[1], argu->za, argu);
	double fxyzb = fxyz(x[0], x[1], argu->zb, argu);
	fval[0] = fxyza-fxyzb;
	return 0;
}

//
// Prism integration
// 	returns the three components (x,y,z) of the anomaly produced by the prism
// 	Input
// 		inf and sup limits of integration in x and y
// 		argu
// 	Output
// 		valx,valy,valz	result of integration for the 3 components
//
void integrate3Comp(double xinf, double xsup, double yinf, double ysup, Argu argu, \
		double *valx, double *valy, double *valz)
{
	double xmin[2], xmax[2], rrerr, err;
	xmin[0] = xinf;
	xmax[0] = xsup;
	xmin[1] = yinf;
	xmax[1] = ysup;

	rrerr = 1; // 1e-3 et 1e-6 : mm resultat

	//
	// Use of function hcubature:
	//	int hcubature(unsigned fdim, integrand f, void *fdata,
	//	              unsigned dim, const double *xmin, const double *xmax, 
	//	              size_t maxEval, double reqAbsError, double reqRelError, 
	//	              error_norm norm,
	//	              double *val, double *err);
	//	Input
	//		fdim	dimension of function f to integrate
	//		f	function to integrate
	//		fdata	for any additional information to compute the function
	//		dim	dimension of integration
	//		xmin	array of minimum boundaries of integration, of length dim
	//		xmax	array of maximum boundaries of integration, of length dim
	//		maxEval	maximum number of function evaluations (0 for no limit)
	//		reqAbsError	absolute error requested if maxEval!=0 (0 to ignore)
	//		reqRelError	relative error requested if maxEval!=0 (0 to ignore)
	//		norm	norm used to measure the error and determine convergence properties
	//			(irrelevant and ignored for fdim<=1)
	//	Output
	//		val	array of computed itegral values, of length fdim
	//		err	array of estimated errors, of length fdim
	//
	//

#ifdef DEBUGMAG
	time_t time0, time1, time2, time3;
	time0 = time(NULL);
#endif
	
	// x-component
	argu.comp = 1;
	hcubature(1,fxy,&argu,2,xmin,xmax,1,0,rrerr,ERROR_INDIVIDUAL,valx,&err);

#ifdef DEBUGMAG
	time1 = time(NULL);
	if(difftime(time1,time0)>1)
		fprintf(stderr,"prism %.3f %.3f %.3f %.3f %.3f %.3f\n",xinf,xsup,yinf,ysup,argu.za,argu.zb);
#endif

	// y-component
	argu.comp = 2;
	hcubature(1,fxy,&argu,2,xmin,xmax,1,0,rrerr,ERROR_INDIVIDUAL,valy,&err);

#ifdef DEBUGMAG
	time2 = time(NULL);
	if(difftime(time2,time1)>1)
		fprintf(stderr,"prism %.3f %.3f %.3f %.3f %.3f %.3f\n",xinf,xsup,yinf,ysup,argu.za,argu.zb);
#endif
	// z-component
	argu.comp = 3;
	hcubature(1,fxy,&argu,2,xmin,xmax,1,0,rrerr,ERROR_INDIVIDUAL,valz,&err);
#ifdef DEBUGMAG
	time3 = time(NULL);
	if(difftime(time3,time2)>1)
		fprintf(stderr,"prism %.3f %.3f %.3f %.3f %.3f %.3f\n",xinf,xsup,yinf,ysup,argu.za,argu.zb);
#endif
}

//
// Project anomalous field on the direction of the regional field
//
double project3Comp(double compx, double compy, double compz)
{
	double proj;

	extern DMagnet *_dmagnet;
	DMagnet *u;
	u = _dmagnet;
	
	proj = compx*u->Ux+compy*u->Uy+compz*u->Uz;
	return(proj);
}

