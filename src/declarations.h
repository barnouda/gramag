typedef unsigned char uchar;

// code used to describe which integration step must be used
#define DX_REFINED 0 	// initial prism  inside velocity/density grid, grid step to be refined
#define DX1 1		// outside velocity/density grid
#define DX2 2		// inside velocity/density grid, refined grid step

typedef struct _Prism {
	double dx1,dx2;		// prism width(s) out and in the mesh volume
	float  *X,*Y;
	int    NPrism;		// number of prisms
	uchar *idxy;		// for each prism, OUT: don't use, DX1 use dx1; DX2 use DX2
} Prism;

typedef struct _Par {
	double	ymin,ymax,xmin,xmax,zmin,zmax;		// limites du modele
//	int 	nx,ny,nz;				// discretisation pour l'integration (nz inutile)
//	double	dxy1,dxy2;				// par d'integration zone 1 (modele vitesse) et zone2 (hors modele) 
//	double 	dx,dy,dz;				// pas d'integration 
//	double	lat,lon;				// ?
	char	topoFile[132];				// fichier bianire decrivant la topo
	char	gridparamFile[132];			// fichier binaire decrivant les parametres de la grille du modele de vitesse
	char	velocityFile[132];			// fichier binaire decrivant le modele de vitese
	char	receiverFile[132];			// fichier ascii contenant les points d'evaluation
	int      nlayer;				// outside grid, density varies with depth, number of layers
	double	*zlayer;
	double   *vmz;					// array of vmass in depth layers
	char    prefix[120];				//prefix for output files
} Par;
//
// topography information
// Convention is elev >0 above water
//
typedef struct {
	double	xmin,xmax,ymin,ymax,tmin,tmax;
	int	nx,ny;
	double	dx,dy;
	double	*elev;
} Topo;
//
// receiver or measurement points
//
typedef struct {
	int	nr;
	double	*x,*y,*z;
} Rec;
//
// topography on grid nodes
//
typedef struct {
	int	nx,ny;
	double *topo;
} TopoGrid;

//
// velocity model
//
typedef struct {
double	xmin,ymin,zmin;	// origin of velocity model
int	nx,ny,nz;	// number of grid point
double	dx,dy,dz;	// grid step
double	*vel;		// array of velocities
} Vel;

//
// volumic mass 
//
typedef struct {
double  xxmin,xxmax,yymin,yymax,zzmin,zzmax;	// model limit for integration
double	xmin,ymin,zmin,xmax,ymax,zmax;		// grid limit for interpolation
double  rhoatzlim;				// volumic mass above zlim is constant
int	    nx,ny,nz;				// number of nodes
double	dx,dy,dz;				// grid step
double	*vm;					// array of values
int      nlayer;				// outside grid, density varies with depth, number of layers
double	*zlayer;
double   *vmz;					// array of vmass in depth layers
} VMass;
//
// node position given its indices
//
static inline double xvmpos(int i, VMass *vm)
{
	return(vm->xmin+i*vm->dx);
}
static inline double yvmpos(int i, VMass *vm)
{
	return(vm->ymin+i*vm->dy);
}
static inline double zvmpos(int i, VMass *vm)
{
	return(vm->zmin+i*vm->dz);
}
// ===================================================================================
//
// Prototypes (from gravilib)
//
void 	usage_and_exit(char *argv);
void    getOptions(int argv, char **argc, int *assignLayer,int *projection, int*checkReceivers, char** inputFile);
double 	elev(double x, double y, Topo*t);
Topo	*topography(Par*par);
Rec 	*readReceiver(Par*par);
Vel	*readVelocity(char *parfile, char *file);
void 	adjustReceiverPosition(Rec *r,Topo *t, Par *p);
VMass 	*createVMassGrid(Par *par,Vel *vel, Topo *t);
double 	interpolateVMass(VMass *vm,double x,double y,double z,int *id,double *alpha);
double 	interpolateVMass0(VMass *vm,double x,double y,double z);
void 	read_param(Par*par, FILE*);
double 	getVmassFromDepth(double z, VMass *v);
void  	getLayerFromDepth(double z, double *zup, double *zdown,VMass *vG);
double 	NafeAndDrake(double);
void    checkGrids(VMass*v, Topo*t, Prism *prism);
int     noTopo(void);



double 	primJ1(double x, double a, double b);
double 	primJ4(double x, double a);
double 	primJ3( double x, double a, double b);
double 	J1(double Xm,double XM,double Y,double Z);
double 	J2(double Xm,double XM,double Y,double Z);
double 	J3(double Xm,double XM,double Y,double Z);
double 	J4(double Xm,double XM,double Y);
double 	I1(double Xm, double XM, double Ym , double YM, double Z);
double 	I2(double Xm, double XM, double Ym , double YM, double Z);

void dumpVMassSlice(Par *p, VMass *v);
void ielev(double x, double y, Topo*t,int*i, int*j);

//
// CONSTANTS
//
#define NOLAYER -999999.
#define BELOW  99999.
#define HALFSPACE 9999999
#define ABOVE -99999.

//
// default values for global
//
#define ASSIGNLAYER 0 // default value for assignLayer. read density(velocity) model from file if 0, or use layer values if 1
#define PROJECTION 1 // default value for projection. project receiver altitude on topography if 1, not if 0
#define CHECKRECEIVERS 1 // default 1 = check if receivers are inside the topography boundaries //AB

