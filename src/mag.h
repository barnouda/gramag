typedef struct
{
	double ux,uy,uz;	// direction of magnetization
	double Ux,Uy,Uz;	// direction of magnetic field (IGRF)
}DMagnet;

typedef struct Argu
{
	double ux,uy,uz;
	double J0,beta;
	double zr,za,zb;
	int comp; // component x=1, y=2, z=3
}Argu;

void read_dir(DMagnet*dmagnet);

double fxyz(double x, double y, double z, Argu *argu);
int fxy(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval);
void integrate3Comp(double xinf, double xsup, double yinf, double ysup, Argu argu, \
		double *valx, double *valy, double *valz);

double project3Comp(double compx, double compy, double compz);

