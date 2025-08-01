// Example from http://ab-initio.mit.edu/wiki/index.php/Cubature
// Compilation : icc -o example hcubature.c example.c -lm
// Vérification avec matlab OK ! (prog verif_example2d.m)

#include <stdio.h>
#include <math.h>
#include "cubature.h"

typedef struct Argu
{
	double sigma;
}Argu;

// The integrand function F should be of the form:
// 	int f(unsigned ndim, const double *x, void *fdata,
//            unsigned fdim, double *fval);
// Input
// 	x	array of evaluation points
// 	ndim	length of x
// 	fdata	= fdata argument passed to hcubature (for additional information)
// Output
// 	fval	array of function values at the point x
// 	fdim	length of fval

// Rq : ici, ndim=3 et fdim=1
int f(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval)
{
	//double sigma = *((double *) fdata); // we can pass σ via fdata argument
	Argu *argu;
	argu = ((Argu *) fdata);
	double sigma = argu->sigma;
	double sum = 0;
	unsigned i;
	//printf("%d %d\n",ndim,fdim);
	for(i = 0; i < ndim; ++i)
	{
		sum += x[i] * x[i];
		//printf("%f\n",x[i]);
	}
	// compute the output value: note that fdim should == 1 from below
	fval[0] = exp(-sigma * sum);
	return 0; // success
}

// Call for h-adaptative cubature:
// 	int hcubature(unsigned fdim, integrand f, void *fdata,
// 	              unsigned dim, const double *xmin, const double *xmax, 
// 	              size_t maxEval, double reqAbsError, double reqRelError, 
// 	              error_norm norm,
// 	              double *val, double *err);

int main(int argc, char **argv)
{
	double xmin[2] = {-2,-2}, xmax[2] = {2,2}, sigma = 0.5, val, err;
	//double xmin[3] = {-2,-2,-2}, xmax[3] = {2,2,2}, sigma = 0.5, val, err;
	Argu argu;
	argu.sigma = 0.5;
	//hcubature(1, f, &sigma, 3, xmin, xmax, 0, 0, 1e-4, ERROR_INDIVIDUAL, &val, &err);
	hcubature(1, f, &argu, 2, xmin, xmax, 0, 0, 1e-4, ERROR_INDIVIDUAL, &val, &err);
	printf("Computed integral = %0.10g +/- %g\n", val, err);
	return 0;
}

