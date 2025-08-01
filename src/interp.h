#ifndef interp_h
#define interp_h
/*
* linear interpolation in 1(?), 2 and 3 dimension
*/

double interp1(double x,double *X,double *V);
double interp2(double x,double y,double *X,double *Y,double *V);
double interp3(double x,double y,double z,double *X,double *Y,double *Z,double *V);
double interp3b(double x,double y,double z,double *X,double *Y,double *Z,double *V,double *alpha);

#endif
