
// Content
/*
 * void    usage_and_exit()
 * void    getOptions(argc, argv, &assignLayer,&projection,&checkReceivers);
 * double  elev(double x, double y, Topo*t)
 * void    ielev(double x, double y, Topo*t,int*i, int*j)
 * Topo   *topography(char*file)
 * Vel    *readVelocity(char *parfile, char *file)
 * Rec    *readReceiver(Par*par)
 * void    adjustReceiverPosition(Rec *r,Topo *t, Par *p)
 * static  inline int vmindex(int i, int j, int k, VMass *vm)
 * double  getVmassFromDepth(double z, VMass *v)
 * void    getLayerFromDepth(double z, double *zup, double *zdown, VMass *v)
 * VMass  *createVMassGrid(Par *par,Vel *vel, Topo *t)
 * double  interpolateVMass(VMass *vm,double x,double y,double z,int *id,double *alpha)
 * void    interpolateVMass2Depth(....)
 * double  interpolateVMass0(VMass *vm,double x,double y,double z)
 * void    read_param(Par*par, FILE *fpi)
 * void    dumpVMassSlice(Par *p, VMass *v)
 *
 * double primJ1(double x, double a, double b)
 * double primJ4(double x, double a)
 * double primJ3( double x, double a, double b)
 * double J1(double Xm,double XM,double Y,double Z)
 * double J2(double Xm,double XM,double Y,double Z)
 * double J3(double Xm,double XM,double Y,double Z)
 * double J4(double Xm,double XM,double Y)
 * double I1(double Xm, double XM, double Ym , double YM, double Z)
 * double I2(double Xm, double XM, double Ym , double YM, double Z)
 */

/*
 Density grid:

nx=3
         <-dx ->
        |   *   |   *   |   *   |
       xmin                    xmax

x:
i=0		i=1		i=3
xmin+dx/2, 	xmin+dx+dx/2  	xmin+2*dx+dx/2



*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
//#include <omp.h>
#include "interp.h"
#include "declarations.h"
#include "version.h"

#define NODES 8	// for linear interpolation in 3D
#define G0 2700 // kg per m3
#define Gg 6.67428e-11	// m^3 kg^-1 s^-2

#define ASMALL 1.e-10
#define ERRDEPTH -9999999.

// ===================================================================================

// 			USAGE_AND_EXIT and GETOPTIONS

// ===================================================================================
void usage_and_exit(char *pgm)
{
	fprintf(stderr,"  | %s: version %s\n",pgm,VERSION);
	fprintf(stderr,"  | (Compute vertical component of gravity field due to heterogeneous distribution of density)\n");
	fprintf(stderr,"  | usage:\n");
	fprintf(stderr,"  | %s: [-l][-np][-nc] input_file\n\n",pgm);
	fprintf(stderr,"  | [-l]:\t density values at the (optional) grid nodes are set according to density layers\n");
	fprintf(stderr,"  | [-np]:\t do not project measurement elevation on topography (default: project)\n");
	fprintf(stderr,"  | [-nc]:\t allow to compute gravity field at locations outside of integration domain\n");
	exit(0);
}

void    getOptions(int argc, char **argv, int *assignLayer,int *projection, int* checkReceivers, char **inputFile)
{
        argc--;
        argv++;
        while (argc>0)
        {
                if (argv[0][0]=='-')
		{
			// density values at grid nodes at forced to values given by the density depth profile
                        if (argv[0][1]=='l')
                                *assignLayer=1;
			// do not project measurement point elevation to topography (default is project)
			else if(argv[0][1]=='n' && argv[0][2]=='p' )
				*projection=0;
			// allow to compute gravity field for locations outside of integration domain
			else if(argv[0][1]=='n' && argv[0][2]=='c')
				*checkReceivers=0;

                	argc--;
                	argv++;
		}
		else
		{
			*inputFile=argv[0];
			break;
		}
        }
}

// ===================================================================================

// 			TOPOGRAPHY

// ===================================================================================

// ===================================================================================
/*
 * return elevation corresponding to longitude x and latitude y
 * (considering the z downward convention, we should rather say depth instead of elevation)
 *
 * input:
 *	x= longitude UTM meter
 *	y= latitude UTM meter
 *	t= Topo structure describing topography
 * return: elevation in meter, (z pointing downward, negative if above sea level
 *         positive otherwise)
 *	   if position is out of topo lateral definition, output ERRDEPTH and issue a message
 */
double elev(double x, double y, Topo*t)
{
int	i,j;
double	z,X[4],Y[4],Z[4];
extern VMass *_vmass;
extern Par     *_par;

	// no topo, the elevation is constant and set to the upper integration limit
	// given in parameters: zmin
	if (t->elev==NULL)
		return _par->zmin;

	if (x<t->xmin || x>t->xmax || y<t->ymin || y>t->ymax)
	{
		fprintf(stderr," (x,y) position is out of topography domain bounds\n");
		fprintf(stderr," x %f [%f %f ]\n",x,t->xmin,t->xmax);
		fprintf(stderr," y %f [%f %f ]\n",y,t->ymin,t->ymax);
		return(ERRDEPTH);
	}
	i=(int)floor((x-t->xmin)/t->dx);
	j=(int)floor((y-t->ymin)/t->dy);
	i=(i==(t->nx-1))? i-1: i;
	j=(j==(t->ny-1))? j-1: j;
	/*
	 * 	0 -- 3
	 * 	|    |     ^
	 *	1 -- 2     |_ 
	 */
	X[0]=X[1]=t->xmin+i*t->dx;
	X[2]=X[3]=X[0]+t->dx;
	Y[1]=Y[2]=t->ymin+j*t->dy;
	Y[0]=Y[3]=Y[1]+t->dy;

	Z[1]=t->elev[i+j*t->nx];
	Z[2]=t->elev[i+1+j*t->nx];
	Z[0]=t->elev[i+(j+1)*t->nx];
	Z[3]=t->elev[i+1+(j+1)*t->nx];
/*
	Z[1]=t->elev[j+i*t->ny];
	Z[2]=t->elev[j+(i+1)*t->ny];
	Z[0]=t->elev[(j+1)+i*t->ny];
	Z[3]=t->elev[(j+1)+(i+1)*t->ny];
*/
	z=interp2(x,y,X,Y,Z);

	return(-z);
}
void ielev(double x, double y, Topo*t,int*i, int*j)
{
	*i=(int)floor((x-t->xmin)/t->dx);
	*j=(int)floor((y-t->ymin)/t->dy);
	*i=(*i==(t->nx-1))? *i-1: *i;
	*j=(*j==(t->ny-1))? *j-1: *j;
}

int noTopo(void)
{
extern Topo *_topo;
	if (_topo->elev==NULL)
		return 1;
	else
		return 0;
}
// -----------------------------------------------------------------------------------

double pgcd(double a, double b)
{
int i,j;
	i=lrint(a);
	j=lrint(b);
	while (i!=j)
		if (i>j)
			i-=j;
		else
			j-=i;
	return ((double)i);
}
/*
 * read the topography file created by createMNT.m
 */
Topo*topography(Par *par)
{
char *file;
int     nd,i;
FILE	*fp;
Topo    *t;
double	dxy;
float   f,*tmp;
int     r;

	file=par->topoFile;
	t=(Topo*)malloc(sizeof(Topo));

	if (strcmp(file,"none")==0 || strcmp(file,"None")==0)
	{
		t->elev=NULL;
		t->dx=1;//pgcd(par->xmax-par->xmin,par->ymax-par->ymin);
fprintf(stderr,"pas de topo mais pas de topo fixé à %f\n",t->dx);
		return t;
	}

	fp=fopen(file,"r");
	if (fp==NULL)
	{
		fprintf(stderr,"could not open topography file <%s>\n",file);
		exit(-1);
	}
	r = fread(&f,sizeof(float),1,fp);
	t->ymin=(double)f;
	r = fread(&f,sizeof(float),1,fp);
	t->ymax=(double)f;
	r = fread(&f,sizeof(float),1,fp);
	t->xmin=(double)f;
	r = fread(&f,sizeof(float),1,fp);
	t->xmax=(double)f;
	fprintf(stderr,"topography range: [%.2f %.2f][%.2f %.2f]\n",t->xmin,t->xmax,t->ymin,t->ymax);
	r = fread(&f,sizeof(float),1,fp);	
	dxy=(double)f;
	fprintf(stderr,"topography grid step %f  ",dxy);
	t->dx=t->dy=dxy;
	r = fread(&t->nx,sizeof(int),1,fp);		// number of node along x direction
	r = fread(&t->ny,sizeof(int),1,fp);		// number of node along y direction
	fprintf(stderr,"read topography grid [%d,%d]  ",t->nx,t->ny);
	tmp=(float*)malloc(sizeof(float)*t->nx*t->ny);
	t->elev=(double*)malloc(sizeof(double)*t->nx*t->ny);
	if (t->elev==NULL || tmp==NULL)
	{
		fprintf(stderr,"could not allocate memory for topography\n");
		exit(-1);
	}
	nd=fread(tmp,sizeof(float),t->nx*t->ny,fp);// matrix of elevations
	if (nd!=t->nx*t->ny)
	{
		fprintf(stderr,"could not read topography data\n");
	}
	fclose(fp);
	fprintf(stderr,"done\n");
	t->tmin=t->tmax=(double)tmp[0];
	for (i=0;i<t->nx*t->ny;i++)
	{
		t->elev[i]=(double)tmp[i];
		t->tmin=(t->elev[i]<t->tmin)? t->elev[i]: t->tmin;	// minimum elevation
		t->tmax=(t->elev[i]>t->tmax)? t->elev[i]: t->tmax;	// maximum elevation
	}
	fprintf(stderr,"topo: min elevation: %.2f; max elevation: %.2f (upward, origin at sea level)\n",t->tmin,t->tmax);
	t->tmin=-t->tmin;
	t->tmax=-t->tmax;
	free(tmp);
	return(t);	
}

// ===================================================================================

// 			VELOCITY MODEL

// ===================================================================================
/*
 * read velocity file that contains the velocity model obtain from 
 * TLR3
 * file .par contains model parameters:
 * 1
 *  642600.000000 1773800.000000   -1500.000000
 *  15  15  11
 *     100.000000     100.000000     100.000000
 *  modelP
 *          1           1           1
 *         13          13          11
 *
 * !!!!!! convention pour les Z: positif vers le bas, <0 si audessus de l'eau, >0 en profondeur
 */
Vel*readVelocity(char *parfile, char *file)
{
int     nd,i;
FILE	*fp;
Vel	*v;
float   *tmp;
extern int assignLayer;
int     return_value;

	if (strcmp(file,"none")==0 || strcmp(file,"None")==0)
	{// No grid, use constant or density varying with depth only
		v=NULL;
		fprintf(stderr,"no velocity grid\n");
		return(v);
	}
	v=(Vel*)malloc(sizeof(Vel));
	// read parameters
	fp=fopen(parfile,"r");
	if (fp==NULL)
	{
		fprintf(stderr,"could not open velocity file <%s>\n",parfile);
		exit(-1);
	}
	return_value = fscanf(fp,"%*d");
	return_value = fscanf(fp,"%lf %lf %lf",&v->xmin,&v->ymin,&v->zmin);
	return_value = fscanf(fp,"%d %d %d",&v->nx,&v->ny,&v->nz);
	return_value = fscanf(fp,"%lf %lf %lf",&v->dx,&v->dy,&v->dz);
	fclose(fp);
	fprintf(stderr,"velocity model:\n");
	fprintf(stderr,"%lf %lf %lf\n",v->xmin,v->ymin,v->zmin);
	fprintf(stderr,"%d %d %d\n",v->nx,v->ny,v->nz);
	fprintf(stderr,"%lf %lf %lf\n",v->dx,v->dy,v->dz);	

	// read model
	// If the density are set equal to the layer model at nodes (-l option) we actually
	// do not need to read the model //AB : seg fault dans outputDerivative sinon
	//if (assignLayer)
	//	return(v);
	// otherwise, read model
	fp=fopen(file,"r");
	if (fp==NULL)
	{
		fprintf(stderr,"could not open velocity file <%s>\n",file);
		exit(-1);
	}
	tmp=(float*)malloc(sizeof(float)*v->nx*v->ny*v->nz);
	v->vel=(double*)malloc(sizeof(double)*v->nx*v->ny*v->nz);
	nd=fread(tmp,sizeof(float),v->nx*v->ny*v->nz,fp);// matrix of velocities
	if (nd!=v->nx*v->ny*v->nz)
	{
		fprintf(stderr,"could not read velocity data\n");
	}
	fclose(fp);
	for (i=0;i<v->nx*v->ny*v->nz;i++)
		 v->vel[i]=tmp[i];
	free(tmp);
	return(v);		
}
static inline int velindex(int i, int j, int k, Vel *v)
{
// fast axis is k, then j, then i
//	return(k + (j*v->nz) + i*(v->nz*v->ny));
// fast axis is i, then j, then k
	return(i+j*v->nx+k*v->nx*v->ny);
}
// ===================================================================================

// 			RECEIVER (or EVALUATION) points

// ===================================================================================
/*
 * read ascii file that contains measurement points
 * format:
 * 		x,y,z
 */
Rec *readReceiver(Par*par)
{
int	i,nr;
double   x,y,z;
Rec	*rec;
FILE 	*fp;
extern int checkReceivers;
int     r;
	fp=fopen(par->receiverFile,"r");
	if (fp==NULL)
	{
		fprintf(stderr,"could not open receiver file <%s>\n",par->receiverFile);
		exit(-1);
	}
	// Count the number of receiver points that are located inside the model limit 
	nr=0;
	while(fscanf(fp,"%lf %lf %lf",&x,&y,&z)!=EOF)
	{
		if(checkReceivers)
		{
			if ( x<par->xmin || x>par->xmax || \
		     	     y<par->ymin || y>par->ymax )
				continue;
		}
		nr++;
	}
	if (nr==0)
	{
		fprintf(stderr,"no receiver read, check model bound or receiver file\n");
		exit(-1);
	}
	rewind(fp);
	rec=(Rec*)malloc(sizeof(Rec));
	rec->nr=nr;
	rec->x=(double*)malloc(sizeof(double)*nr);
	rec->y=(double*)malloc(sizeof(double)*nr);
	rec->z=(double*)malloc(sizeof(double)*nr);
	fprintf(stderr,"reading positions for %d receivers\n",nr);
	for (i=0;i<nr;)
	{
		r = fscanf(fp,"%lf %lf %lf",&rec->x[i],&rec->y[i],&rec->z[i]);
		if(checkReceivers)
		{
			if ( rec->x[i]<par->xmin || rec->x[i]>par->xmax || \
		             rec->y[i]<par->ymin || rec->y[i]>par->ymax )
				continue;
		}
//		fprintf(stderr,"%f %f %f\n",rec->x[i],rec->y[i],rec->z[i]);
		i++;
	}
	fclose(fp);
	return(rec);
}

// ===================================================================================
/*
 * adjust the measurement points to the density and topography grid
 * Assuming thet the grid step is small enough, positions are moved to the center
 * of the element. We introduce an error that is at most dx/2, dy/2
 * This is done to make gravity integration easier
 * 
 * !!!!! convention pour Z en entree: 
 *						negatif toward top,    <0 if below sea level
 *						positif toward bottom, >0 if below sea level
 */
void
adjustReceiverPosition(Rec *r,Topo *t, Par *p)
{
int	i;
double	z;
FILE    *fp;
extern int projection;
extern int checkReceivers;
	fp=fopen("receiverMoved.asc","w");
	for (i=0;i<r->nr;i++)
	{
	
		if(checkReceivers)
		{	
			// x position, check for bounds
			if (r->x[i]>p->xmax || r->x[i]<p->xmin)
			{
				fprintf(stderr,"receiver %d out of X (lon) bound (%.2f)\n",i,r->x[i]);
				exit(-1);
			}

//			ix=floor((r->x[i]-p->xmin)/p->dx);
//			r->x[i]=ix*p->dx+p->dx/2+p->xmin;

			// y position, check for bounds
			if (r->y[i]>p->ymax || r->y[i]<p->ymin)
			{
				fprintf(stderr,"receiver %d out of Y (lat) bound (%.2f)\n",i,r->y[i]);
				exit(-1);
			}
//			iy=floor((r->y[i]-p->ymin)/p->dy);
//			r->y[i]=iy*p->dy+p->dy/2+p->ymin;
		}

		if (!projection)
		{
			z=r->z[i];
			r->z[i]*=-1;
		}
		else
		{
			// z position, project on topography
			// and reverse z axis convention, from positive upward to positive downward
			z=r->z[i];
			r->z[i]=elev(r->x[i],r->y[i],t);
		}

//		fprintf(stderr,"receiver %d moved to %.2f %.2f %.2f (%.2f)\n",i,r->x[i],r->y[i],-r->z[i],z);
		fprintf(fp,"%f %f %f\n",r->x[i],r->y[i],-r->z[i]);
if (1==2)
{
FILE *fp2=fopen("topo.xyz","w");
int i,j;
double x,y;
for (i=505382;i<818922;i+=500)
   for (j=1611303;j<1929803;j+=500)
       {x=i;y=j;
           fprintf(fp2,"%.2f %.2f %.2f\n",x,y,elev(x,y,t));
       }
fclose(fp2);       
exit(-1);
}
	}
	fclose(fp);
}
// ===================================================================================

// 			VOLUMIC MASS GRID

// ===================================================================================
/*
 * create the volumic mass grid based on velocity grid
 * This grid will interpolate density on a volume corresponding to the velocity grid
 * Outside this velocity grid, density is constant or obeys simple relationship (eg with respect to depth)
 * defined in "par" structure
 *
 * The initial grid is as follow, with constant velocity inside each element
 * nx=3
 *              xmin                    xmax
 *              ||  V(1) |  V(2) |  V(3) ||
 *
 *
 * From that, we create the volumic mass grid where volumic mass is prescribed at node position.
 * Two more nodes are added in order to interpolate between outside medium and initial grid
 *
 * nx=5
 *      |   *   ||   *   |   *   |   *   ||   *   |
 *
 * node position is xmin+dx/2, xmin+dx+dx/2,...
 *
 * and interpolated linearly between the nodes
 */
/*
 * Nafe and Drake, 1963 law between density and VP, taken from Onizawa, Mikada, Watanabe et Sakashita, Earth Planet Space, VOL54, 2002
 */

//
// 3D addressing in vmass->vmass array
//
static inline int vmindex(int i, int j, int k, VMass *vm)
{
// fast axis is i, then j, then k
	return(i+j*vm->nx+k*vm->nx*vm->ny);
}

/*
 * volumic mass outside grid may be constant, or varying according to depth
 * These two routines return vmass and layers corresponding to depth z
 */
double
getVmassFromDepth(double z, VMass *v)
{
	int i;
	if (v->nlayer==0)	// no layers, return constant volumic mass
		return(v->vmz[0]);
	else			// layers, return volumic mass in the corresponding layer
	{
		i=0;
		while (z>v->zlayer[i])
		{
			i++;
			if (i==v->nlayer)
				break;
		}
		return (v->vmz[i]);
	}
}
/*
 * return layer information such that:   *zdown <= z < *zup
 */
void
getLayerFromDepth(double z, double *zup, double *zdown, VMass *v)
{
	int i;
	if (v->nlayer==0)
	{
		*zup=*zdown=NOLAYER;//ModifOC1
		*zup=NOLAYER;
		*zdown=v->zzmax;
	}
	else
	{
		for (i=0;i<v->nlayer;i++)
			if (z<v->zlayer[i])
				break;
		if (i==0)			// in first layer
		{
			*zup=v->zzmin;
			*zdown=v->zlayer[i];
		} 
		else if (i==(v->nlayer-1))	// in last layer
		{
			*zup=v->zlayer[i-1];
			*zdown=v->zzmax;
		}
		else
		{
			*zup=v->zlayer[i-1];
			*zdown=v->zlayer[i];
		}		
	}

}

/*
 * create the vmass grid
 */
VMass *createVMassGrid(Par *par,Vel *vel, Topo *t)
{
int	i,j,k;
double   dx2,dy2,dz2,z0;
double  x,y,z;
VMass   *vm;
extern int assignLayer;

	vm=(VMass*)malloc(sizeof(VMass));
	// copy depth varying vmass information
	vm->nlayer=par->nlayer;
	vm->zlayer=par->zlayer;
	vm->vmz=par->vmz;

	// if no velocity model exist, fill in required fields
	if (vel==NULL)
	{
		// these are the integration limits
		// zzmin may also represent the topography constant elevation 
		// in case no topography files are given
		vm->xxmax=par->xmax;
		vm->xxmin=par->xmin;
		vm->yymax=par->ymax;
		vm->yymin=par->ymin;
		vm->zzmax=par->zmax;
		vm->zzmin=par->zmin;
		// define virtual grid limit so that a test for point
		// inside grid always return false
//		vm->xmax=vm->xxmin-par->dx;
//		vm->xmin=vm->xxmax+par->dx;
//		vm->ymax=vm->yymin-par->dy;
//		vm->ymin=vm->yymax+par->dy;
//		vm->zmax=vm->zzmin-par->dz;
//		vm->zmin=vm->zzmax+par->dz;
		vm->vm=NULL;
		return(vm);
	}
	// keep the same grid spacing as velocity model
	vm->dx=vel->dx;
	vm->dy=vel->dy;
	vm->dz=vel->dz;
	dx2=vm->dx/2.;
	dy2=vm->dy/2.;
	dz2=vm->dz/2.;

	// extend velocity grid one node outside to interpolate toward constant model
	// xmin, xmax = position of nodes at each ends
	vm->xmin=vel->xmin                     - vel->dx;
	vm->xmax=vel->xmin+(vel->nx-1)*vel->dx + vel->dx;

	vm->ymin=vel->ymin                     - vel->dy;
	vm->ymax=vel->ymin+(vel->ny-1)*vel->dy + vel->dy;
	
	vm->zmin=vel->zmin                     - vel->dz;
	vm->zmax=vel->zmin+(vel->nz-1)*vel->dz + vel->dz;

	// number of node in each direction
	vm->nx=(vm->xmax-vm->xmin)/vm->dx+1;
	vm->ny=(vm->ymax-vm->ymin)/vm->dy+1;
	vm->nz=(vm->zmax-vm->zmin)/vm->dz+1;
	fprintf(stderr,"density grid extended to (x):%.2f %.2f\n",vm->xmin,vm->xmax);
	fprintf(stderr,"density grid extended to (y):%.2f %.2f\n",vm->ymin,vm->ymax);
	fprintf(stderr,"density grid extended to (z):%.2f %.2f\n",vm->zmin,vm->zmax);
	fprintf(stderr,"density grid size: %d x %d x %d\n",vm->nx,vm->ny,vm->nz);

	// size of complete density model as defined by parameters
	vm->xxmax=par->xmax;
	vm->xxmin=par->xmin;
	vm->yymax=par->ymax;
	vm->yymin=par->ymin;
	vm->zzmax=par->zmax;
	vm->zzmin=par->zmin;

	// create matrix of mass volumic values
	vm->vm=(double*)malloc(sizeof(double)*vm->nx*vm->ny*vm->nz);

	// initialize volumic mass at each node to the value deduced from 
	// the model out of grid (constant density, or density constant per layer)
	// Node located above topography are initialized to the value at topography
	for (k=0;k<vm->nz;k++)		//depth
	{
		z=vm->zmin+k*vm->dz;
		for (j=0;j<vm->ny;j++)	//latitude
		{
			y=vm->ymin+j*vm->dy;
			for (i=0;i<vm->nx;i++)	//longitude
			{
				x=vm->xmin+i*vm->dx;
				z0=elev(x,y,t);
				if (z0<ERRDEPTH+1)
				{
					fprintf(stderr,"Abort: vmass grid has been extended out of topography limits\n");
					fprintf(stderr,"       Increase topography limits or reduce grid limits\n");
					exit(-1);

				}
				if (z<z0)	// above topography (!! z>O downward)
					vm->vm[vmindex(i,j,k,vm)]=getVmassFromDepth(z0,vm);//=0	// A modifier OC 21/05/2013
				else
					vm->vm[vmindex(i,j,k,vm)]=getVmassFromDepth(z,vm);
			}
		}
	}
#ifdef DEBUG
{
FILE *fp=fopen("grid.dat","w");
for (j=0;j<vm->ny;j++)	//latitude
{
	y=vm->ymin+j*vm->dy;
	for (i=0;i<vm->nx;i++)	//longitude
	{
		x=vm->xmin+i*vm->dx;
		fprintf(fp,"%.2f %.2f %d %d\n",x,y,i,j);
	}
}
fclose(fp);
}
#endif

	vm->rhoatzlim=G0;//Usefull?

	if (assignLayer==0)
	// copy values deduced from density model
	// volumic mass grid has 1 more node on each side
	{
 	double d,v;

	for (k=1;k<vm->nz-1;k++)			//depth
	{	
		for (j=1;j<vm->ny-1;j++)		//latitude
		{	
			for (i=1;i<vm->nx-1;i++)	//longitude
			{
				d=vel->vel[velindex(i-1,j-1,k-1,vel)];
				vm->vm[vmindex(i,j,k,vm)]=d;

			}
		}
	}	
	}
	return(vm);
}


// ===================================================================================
/*
 * Interpolate volumic mass from grid at position x,y,z
 * We should never arrive in this routine if no grid is defined (vm->vm==NULL)
 * input:
 *	vm: structure describing volumic mass grid
 *	x,y,z: position at which we interpolate
 *	id: array of node indices used for interpolation
 *	alpha: array of interpolation coefficient for indices id
 * output:
 * 	nnodes: number of nodes used for interpolation (eg 8 for linear)
 * return:
 *	volumic mass interpolated at (x,y,z)
 */	
double interpolateVMass(VMass *vm,double x,double y,double z,int *id,double *alpha)
{
// use linear interpolation
#define NNODES 8
int	i,j,k;
double 	X[NNODES],Y[NNODES],Z[NNODES],VM[NNODES],rho;	

	// check wether (x,y,z) is in/outside of grid
	// outside grid, vmass varies according to depth
//fprintf(stderr,"%.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f\n",x,y,z,vm->xmin,vm->xmax,vm->ymin,vm->ymax,vm->zmin,vm->zmax);
	if (x<vm->xmin || x>vm->xmax || y<vm->ymin || y>vm->ymax || z<vm->zmin || z>vm->zmax )
	{
		id[0]=-1; //means that point was out of grid
		return(getVmassFromDepth(z,vm));
	}
	// if point is located on an exterior face, move it slighly inside the grid
	if (z==vm->zmax)
		z-=vm->dz*ASMALL;
	if (y==vm->ymax)
		y-=vm->dy*ASMALL;
	if (x==vm->xmax)
		x-=vm->dx*ASMALL;

	// inside grid, vmass is interpolated between nodes
	i=(int)floor((x-vm->xmin)/vm->dx);
	j=(int)floor((y-vm->ymin)/vm->dy);
	k=(int)floor((z-vm->zmin)/vm->dz);

	id[0]=vmindex(i,j,k,vm);
	id[1]=vmindex(i,j,k+1,vm);
	id[2]=vmindex(i,j+1,k,vm);
	id[3]=vmindex(i,j+1,k+1,vm);
	id[4]=vmindex(i+1,j,k,vm);
	id[5]=vmindex(i+1,j,k+1,vm);
	id[6]=vmindex(i+1,j+1,k,vm);
	id[7]=vmindex(i+1,j+1,k+1,vm);

	X[0]=X[1]=X[2]=X[3]=xvmpos(i,vm);
	X[4]=X[5]=X[6]=X[7]=xvmpos(i+1,vm);
	Y[0]=Y[1]=Y[4]=Y[5]=yvmpos(j,vm);
	Y[2]=Y[3]=Y[6]=Y[7]=yvmpos(j+1,vm);
	Z[0]=Z[2]=Z[4]=Z[6]=zvmpos(k,vm);
	Z[1]=Z[3]=Z[5]=Z[7]=zvmpos(k+1,vm);

        for (i=0;i<8;i++)
                VM[i]=vm->vm[id[i]];

	rho=interp3b(x,y,z,X,Y,Z,VM,alpha);

	return(rho);
}
//same as above but for 2 points, sharing same (x,y) coordinates, different z but in the same cell
void interpolateVMass2Depth(VMass *vm,double x,double y,double z1, double z2,int *id, \
	                     double *alpha1, double *alpha2,double *rho1, double *rho2)
{
// use linear interpolation
#define NNODES 8
int	i,j,k;
double 	X[NNODES],Y[NNODES],Z[NNODES],VM[NNODES];	


	// if point is located on an exterior face, move it slighly inside the grid
	if (z1==vm->zmax)
		z1-=vm->dz*ASMALL;
	if (y==vm->ymax)
		y-=vm->dy*ASMALL;
	if (x==vm->xmax)
		x-=vm->dx*ASMALL;

	// inside grid, vmass is interpolated between nodes
	i=(int)floor((x-vm->xmin)/vm->dx);
	j=(int)floor((y-vm->ymin)/vm->dy);
	k=(int)floor((z1-vm->zmin)/vm->dz);

	id[0]=vmindex(i,j,k,vm);
	id[1]=vmindex(i,j,k+1,vm);
	id[2]=vmindex(i,j+1,k,vm);
	id[3]=vmindex(i,j+1,k+1,vm);
	id[4]=vmindex(i+1,j,k,vm);
	id[5]=vmindex(i+1,j,k+1,vm);
	id[6]=vmindex(i+1,j+1,k,vm);
	id[7]=vmindex(i+1,j+1,k+1,vm);

	X[0]=X[1]=X[2]=X[3]=xvmpos(i,vm);
	X[4]=X[5]=X[6]=X[7]=xvmpos(i+1,vm);
	Y[0]=Y[1]=Y[4]=Y[5]=yvmpos(j,vm);
	Y[2]=Y[3]=Y[6]=Y[7]=yvmpos(j+1,vm);
	Z[0]=Z[2]=Z[4]=Z[6]=zvmpos(k,vm);
	Z[1]=Z[3]=Z[5]=Z[7]=zvmpos(k+1,vm);

        for (i=0;i<8;i++)
                VM[i]=vm->vm[id[i]];


	*rho1=(double)interp3b(x,y,z1,X,Y,Z,VM,alpha1);
	*rho2=(double)interp3b(x,y,z2,X,Y,Z,VM,alpha2);


}
double interpolateVMass0(VMass *vm,double x,double y,double z)
{
// use linear interpolation
#define NNODES 8
int	i,j,k;
double 	X[NNODES],Y[NNODES],Z[NNODES],VM[NNODES],rho,alpha[8];
int     id[8];


	// check wether (x,y,z) is in/outside of grid
	if (x<vm->xmin || x>vm->xmax || y<vm->ymin || y>vm->ymax || z<vm->zmin || z>vm->zmax )
	{
		return(getVmassFromDepth(z,vm));
	}
	// if point is located on an exterior face, move it slighly inside the grid
	if (z==vm->zmax)
	{
		z-=vm->dz*ASMALL;
	}
	if (y==vm->ymax)
		y-=vm->dy*ASMALL;
	if (x==vm->xmax)
		x-=vm->dx*ASMALL;
	
	i=(int)floor((x-vm->xmin)/vm->dx);
	j=(int)floor((y-vm->ymin)/vm->dy);
	k=(int)floor((z-vm->zmin)/vm->dz);

	id[0]=vmindex(i,j,k,vm);
	id[1]=vmindex(i,j,k+1,vm);
	id[2]=vmindex(i,j+1,k,vm);
	id[3]=vmindex(i,j+1,k+1,vm);
	id[4]=vmindex(i+1,j,k,vm);
	id[5]=vmindex(i+1,j,k+1,vm);
	id[6]=vmindex(i+1,j+1,k,vm);
	id[7]=vmindex(i+1,j+1,k+1,vm);

	X[0]=X[1]=X[2]=X[3]=xvmpos(i,vm);
	X[4]=X[5]=X[6]=X[7]=xvmpos(i+1,vm);
	Y[0]=Y[1]=Y[4]=Y[5]=yvmpos(j,vm);
	Y[2]=Y[3]=Y[6]=Y[7]=yvmpos(j+1,vm);
	Z[0]=Z[2]=Z[4]=Z[6]=zvmpos(k,vm);
	Z[1]=Z[3]=Z[5]=Z[7]=zvmpos(k+1,vm);

        for (i=0;i<8;i++)
                VM[i]=vm->vm[id[i]];

	rho=interp3b(x,y,z,X,Y,Z,VM,alpha);

	return(rho);
}

// ===================================================================================

//			INPUT PARAMETER

// ===================================================================================
/*
 * read computation parameters
 * grid parameter used for inegration
 *
 * 
 */
void read_param(Par*par, FILE *fpi)
{
	int r;
	printf("prefix:  ");
	r = fscanf(fpi,"%s",par->prefix);
	printf("%s\n",par->prefix);
	printf("topoFile:  ");
	r = fscanf(fpi,"%s",par->topoFile);
	printf("%s\n",par->topoFile);
	printf("gridparamFile:  ");
	r = fscanf(fpi,"%s",par->gridparamFile);
	printf("%s\n",par->gridparamFile);
	printf("velocityFile:  ");
	r = fscanf(fpi,"%s",par->velocityFile);
	printf("%s\n",par->velocityFile);
	printf("receiverFile:  ");
	r = fscanf(fpi,"%s",par->receiverFile);
	printf("%s\n",par->receiverFile);
	printf("xmin, xmax  ");
	r = fscanf(fpi,"%lf %lf",&par->xmin,&par->xmax);
	printf("%f %f\n",par->xmin,par->xmax);
	printf("ymin, ymax  ");
	r = fscanf(fpi,"%lf %lf",&par->ymin,&par->ymax);
	printf("%f %f\n",par->ymin,par->ymax);
	printf("zmin, zmax  ");
	r = fscanf(fpi,"%lf %lf",&par->zmin,&par->zmax);
	printf("%f %f\n",par->zmin,par->zmax);
//	printf("dx, dy  ");
//	fscanf("%lf %lf",&par->dx,&par->dy);
//	printf("%f %f\n",par->dx,par->dy);
//	printf("lat, lon  ");
//	fscanf("%lf %lf",&par->lat,&par->lon);
//	printf("%f %f\n",par->lat,par->lon);
	printf("number of layer (0= none)\n");
	r = fscanf(fpi,"%d",&par->nlayer);
	printf("%d\n",par->nlayer);
	par->vmz=(double*)malloc(sizeof(double)*(par->nlayer));
	par->zlayer=(double*)malloc(sizeof(double)*(par->nlayer));
	if (par->nlayer==0)
	{
		printf("default volumic mass (kg/m3): \n");
		r=fscanf(fpi,"%lf",&par->vmz[0]);
		if (r!=1)
		{
			fprintf(stderr,"read_param: abort: expect volumic mass value\n");
			exit(-1);
		}
		par->zlayer=NULL;
	} else {
		int i;
		
		printf("enter volumic mass (kg/m3) and bottom depth (0 for last layer): \n");
		for (i=0;i<par->nlayer;i++)
		{ 
			printf("layer %d\n",i);
			r=fscanf(fpi,"%lf %lf",&par->vmz[i],&par->zlayer[i]);
			if (r!=2)
			{
				fprintf(stderr,"read_param: abort: expect volumic mass value and bottom depth (0 for last layer)\n");
				exit(-1);
			}
			if (par->zlayer[i]>par->zmax || par->zlayer[i]<par->zmin)
			{
				fprintf(stderr,"read_param: Abort: layer #%d bottom depth (%f) is out of integration domain bounds [%f %f]!!\n",i,par->zlayer[i],par->zmin,par->zmax);
				exit(-1);
			}
		}
	}
	if (par->nlayer>0)
		par->zlayer[par->nlayer-1]=HALFSPACE;

//	par->nx=(par->xmax-par->xmin)/par->dx;
//	par->ny=(par->ymax-par->ymin)/par->dy;
}
// ===================================================================================

//			DUMP VOLUMIC MASS GRID SLICES

// ===================================================================================
/*
 * write plane of density value at constant depth.
 * read input ascii file "density.map" to get parameters
 * input format is depth, dx, dy
 * all in meters
 */
void
dumpVMassSlice(Par *p, VMass *v)
{
FILE *fp,*fpo;
char  file[132];
int   i,j,k,nx,ny;
double dx,dy,alt,x,y,rho;

	fp=fopen("density.map","r");
	if (fp==NULL)
		return;
	k=1;
	while (fscanf(fp,"%lf %lf %lf",&alt,&dx,&dy)!=EOF)
	{
		fprintf(stderr,"plot density section at depth %f\n",alt);
		sprintf(file,"density.%d",k);
		fpo=fopen(file,"w");
		nx=(p->xmax-p->xmin)/dx;
		ny=(p->ymax-p->ymin)/dy;
		fwrite(&nx,sizeof(int),1,fpo);
		fwrite(&ny,sizeof(int),1,fpo);
		fwrite(&alt,sizeof(double),1,fpo);
		for (j=0;j<ny;j++)
		{
			y=p->ymin+j*dy;
			for (i=0;i<nx;i++)
			{
				x=p->xmin+i*dx;
				rho=interpolateVMass0(v,x,y,alt);

				fwrite(&x,sizeof(double),1,fpo);
				fwrite(&y,sizeof(double),1,fpo);
				fwrite(&rho,sizeof(double),1,fpo);
			}
		}
		fclose(fpo);
		k++;
	}
	fclose(fp);
	fp=fopen("density.vol","w");
	fwrite(v->vm,sizeof(double),v->nx*v->ny*v->nz,fp);
	fclose(fp);
}

/*
 * compute grid step for integration.
 * Two steps are set, outside density/velocity grid, and inside
 * dx1: integration grid step out of density grid
 * dx2: integration grid step inside density grid
*/
// return 1 if out of grid, 0 otherwise
double outGrid(double x, double y, double dx, VMass *v)
{
	if (v->vm==NULL)
		return(1);
	if ((x+dx<v->xmin)||(x>v->xmax)||(y+dx<v->ymin)||(y>v->ymax))
		return(1);
	return(0);
}
void
checkGrids(VMass*v, Topo*t, Prism *prism)
{
int	nx,ny,i,j,k,l,nprism,mprism,m,mxy;
double  dx1,dx2;
float  *x,*y,*ptx,*pty;
uchar  *grid,*ptg;
#ifdef DEBUG
int    *index;
#endif

	// integration grid step outside density grid
	if (noTopo())
		dx1=t->dx;///10;
	else
		dx1=20.*t->dx;				//20 times the topography discretization step
	nx=floor((v->xxmax-v->xxmin)/dx1);
	ny=floor((v->yymax-v->yymin)/dx1);
	fprintf(stderr,"integration step outside density/velocity grid: %.2f\n",dx1);

	fprintf(stderr,"checkGrid xmin,xmax %.2f %.2f [%.2f %.2f]\n",v->xxmin,(nx-1)*dx1+v->xxmin+dx1,v->xxmin,v->xxmax);
	fprintf(stderr,"checkGrid ymin,ymax %.2f %.2f [%.2f %.2f]\n",v->yymin,(ny-1)*dx1+v->yymin+dx1,v->yymin,v->yymax);

	// loop for prisms outside density/velocity grid
	mprism=nprism=nx*ny;
	x=(float*)malloc(sizeof(float)*nprism);
	y=(float*)malloc(sizeof(float)*nprism);
	grid=(uchar*)malloc(sizeof(uchar)*nprism);
#ifdef DEBUG
	index=(int*)malloc(sizeof(int)*nprism);
#endif
	for (j=0,l=0;j<ny;j++)
		for (i=0;i<nx;i++,l++)
		{
			//lower left coordinate of prism
			x[l]=i*dx1+v->xxmin;	
			y[l]=j*dx1+v->yymin;
			if (outGrid(x[l],y[l],dx1,v))
				grid[l]=DX1;
			else
 				grid[l]=DX_REFINED;
#ifdef DEBUG
			index[l]=l;
#endif
		}

	// loop for prisms inside density/velocity grid
	if (v->vm!=NULL)
	{
		// integration grid size inside density grid
		// is a fraction of integration grid size outside
		// At least, 1/third of grid step
		m=1;
		dx2=dx1/m;
		while (dx2>0.25*v->dx)
		{
			m++;
			dx2=dx1/m;
		}
		//m--;	//BUG??
		fprintf(stderr,"integration step inside density/velocity grid: %.2f [1/%d]\n",dx2,m);
		mxy=m*m;
		mprism=nprism;	
		for (l=0;l<nprism;l++)
		{
			if (grid[l] != DX_REFINED)
				continue;
			mprism+=mxy;	//add mxy more prisms to integrate
			x=(float*)realloc(x,mprism*sizeof(float));			
			y=(float*)realloc(y,mprism*sizeof(float));			
			grid=(uchar*)realloc(grid,mprism*sizeof(uchar));			
#ifdef DEBUG
			index=(int*)realloc(index,sizeof(int)*mprism);
#endif
			ptx=x+mprism-mxy;
			pty=y+mprism-mxy;
			ptg=grid+mprism-mxy;
			for (j=0,k=0;j<m;j++)
				for (i=0;i<m;i++,k++)
				{
					ptx[k]=x[l]+i*dx2;	// pointer to x array
					pty[k]=y[l]+j*dx2;	// pointer to y array
					ptg[k]=DX2;		// pointer to grid array, DX2 means use dx2
#ifdef DEBUG
					index[mprism-mxy+k]=l;
#endif
				}
		}
	}
/*
	for (l=0;l<mprism;l++)
	{double xx,yy;
		if (grid[l]==0.)
			continue;
		xx=x[l];
		yy=y[l];
		if (xx<v->xxmin || xx+grid[l]>v->xxmax || yy<v->yymin || yy+grid[l]>v->yymax)
		{
			fprintf(stderr,"%d out of bounds %.2f %.2f %.2f\n",l,xx,yy,grid[l]);
			exit(-1);
		}
	}
*/
	prism->X=x;	
	prism->Y=y;
	prism->idxy=grid;
	prism->dx1=dx1;
	prism->dx2=dx2;
	prism->NPrism=mprism;
//
#ifdef DEBUG
{
FILE *fp=fopen("integrationGrid.dat","w");
	for (l=0;l<mprism;l++)
		fprintf(fp,"%.2f %.2f %.2f \n",x[l],y[l],grid[l]);
	fclose(fp);
	fp=fopen("indexPrism.dat","w");
	        for (l=0;l<mprism;l++)
                fprintf(fp,"%d %d %f\n",l,index[l],grid[l]);
        fclose(fp);
	free(index);
}
#endif
//

	return;
}

static int _isnan_=0;
#ifndef isnan
#define isnan(x) ((x) != (x))
#endif

/*
 * primitive de log(sqrt(x^2+y^2+z^2)+y)
 * avec log(sqrt(x^2+y^2+dz^2)+y) = primitive de 1/sqrt(x^2+y^2+dz^2)
 */
double primJ1(double x, double a, double b)
{
	double r,R;
	
		r=sqrt(a*a + b*b + x*x);
		if (a!=0)
			R= -x + a*(atan(x/a) - atan((b*x)/(a*r))) + x*log(b + r) + b*log(2*(x + r));
		else
			R= -x + a*(atan2(a,x) - atan2(a*r,b*x)) + x*log(b + r) + b*log(2*(x + r));

	return(R);
}
	
/*
 * primtive de x*atan(a/x)
 */
double primJ4(double x, double a)
{
	double R;
	
	R=(a*x + (a*a + x*x)*atan2(a,x))/2;
	
	return(R);
}
/*
 * primitive de x*atan(a*b/x/sqrt(x^2+a^2+b^2))
 *
 */
double
primJ3( double x, double a, double b)
{
        double  sgn,r,a2,b2;
        double  two,eigth,four;
        double  RR;

        if (b==0.)
                b=0.00001;
        two=2.;
        eigth=8.;
        four=4.;
        a2=a*a;
        b2=b*b;
        r=sqrt(a2 + b2 + x*x);
        sgn=a*b/fabs(a*b);
        a=fabs(a);
        b=fabs(b);

	if (x!=0)
        	RR=sgn*(x*x*atan((a*b)/(x*r)) + 2*a*b* log(2*(x+r)) +           \
                          b2*atan(b/x) + a2*atan(a/x) -                 \
                          b2*atan((a2+b2+a*r)/(b*x)) - a2*atan((a2+b2+b*r)/(a*x)) )/2;
	else
        	RR=sgn*(x*x*atan2((x*r),(a*b)) + 2*a*b* log(2*(x+r)) +           \
                          b2*atan2(x,b) + a2*atan2(x,a) -                 \
                          b2*atan2((b*x),(a2+b2+a*r)) - a2*atan2((a*x),(a2+b2+b*r)) )/2;
		
        return(RR);

}
#ifdef VERSIONFORT
double primJ3(double x, double a, double b)
{
double r1,r2;
	double primj3_(double*x, double*a, double*b);
	r1=primj3_(&x, &a, &b);
	return(r1);
}
#endif

double J1(double Xm,double XM,double Y,double Z)
{
	double R;
	
	R=Y*(primJ1(XM,Y,Z)-primJ1(Xm,Y,Z));
//fprintf(stderr,"J1: %f %f (%f %f %f %f)\n",primJ1(XM,Y,Z),primJ1(Xm,Y,Z),Xm,XM,Y,Z);	
	return(R);
}
double J2(double Xm,double XM,double Y,double Z)
{
	double R;
	
	R=Z*(primJ1(XM,Z,Y)-primJ1(Xm,Z,Y));
	
	return(R);
}
double J3(double Xm,double XM,double Y,double Z)
{
	double R;
	double primJ3(double,double,double);
	R=primJ3(XM,Y,Z)-primJ3(Xm,Y,Z);
	
if (_isnan_)
fprintf(stderr,"J3 %f %g %g %g\n",primJ3(XM,Y,Z),XM,Y,Z);
	return(R);
}
double J4(double Xm,double XM,double Y)
{
	double R;
	
    R=primJ4(XM,Y)-primJ4(Xm,Y);
	
	return(R);
}
double I1(double Xm, double XM, double Ym , double YM, double Z)
{
	double R,J31,J32;
    J31=J3(Xm,XM,YM,Z);
    J32=J3(Xm,XM,Ym,Z);
    R=(J1(Xm,XM,YM,Z) - J1(Xm,XM,Ym,Z)) + \
      (J2(Xm,XM,YM,Z) - J2(Xm,XM,Ym,Z)) - \
      (J31 - J32) + \
      (J4(Xm,XM,YM)   - J4(Xm,XM,Ym)) - \
      (YM-Ym)*(XM-Xm);
if (_isnan_)
fprintf(stderr,"I1 %f %g %g\n",R,J31,J32);

	return(R);
}
double I2(double Xm, double XM, double Ym , double YM, double Z)
{
	double R;
	
	R=Z*((primJ1(XM,Z,YM)-primJ1(Xm,Z,YM)) - (primJ1(XM,Z,Ym)-primJ1(Xm,Z,Ym)));
//fprintf(stderr,"I2 %f\n",R);	
	
	return(R);
}
