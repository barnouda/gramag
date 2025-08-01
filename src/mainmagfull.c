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
#include <omp.h>
#include "interp.h"

#define NODES 8	// for linear interpolation in 3D
#define G0 2700 // kg per m3
#define Gg 6.67428e-11	// m^3 kg^-1 s^-2
#define Cm 1e-7 // henry/m => mag induction in Tesla //mag
//
// general parameters
//

#include "declarations.h"
#include "mag.h" //mag
//local prototypes
void loop(Prism *prism,double *xr,double *yr,double *zr,int nr,int mxyz,double *G1,double *G2);

Par     *_par;
Topo    *_topo;
VMass   *_vmass;
DMagnet *_dmagnet; //mag

// global
int 	assignLayer=ASSIGNLAYER;
int     projection=PROJECTION;
int     checkReceivers=CHECKRECEIVERS;

void outputDerivatives(double *G1,double *G2,int nr,int mxyz,Vel *vel,VMass *vm,char *prefix);
//			MAIN
// ===================================================================================

int main(int argc, char**argv)
{

char	file[132],*inputFile;
int	i,j;
int	mxyz;		//number of nodes in vmass grid
double	*G1,*G2,G;
Prism	prism;
Par	par;
Rec 	*rec;
Topo	*topo;
Vel	*vel;
VMass	*vmass;
DMagnet	dmagnet;	//mag
FILE    *fp;

        if (argc==1)
                usage_and_exit(argv[0]);

	// read options
	getOptions(argc, argv, &assignLayer,&projection,&checkReceivers,&inputFile);
	if (projection)
	              fprintf(stderr,"warning: receiver elevation are projected on topography\n");
	else
	              fprintf(stderr,"warning: receiver elevation are NOT projected on topography\n");
        if (inputFile==NULL)
                usage_and_exit(argv[0]);

	// read density grid information
        fp=fopen(inputFile,"r");
        if (fp==NULL)
        {
                fprintf(stderr,"wrong input file <%s>\n",inputFile);
                exit(-1);
        }
        read_param(&par,fp);
        fclose(fp);
        _par=&par;

	fprintf(stderr,"read magnetization and field directions\n");	//mag
	read_dir(&dmagnet);						//mag
	_dmagnet=&dmagnet;						//mag

	fprintf(stderr,"read receiver positions\n");
	// read receiver information
	rec=readReceiver(&par);

	fprintf(stderr,"read topography model\n");
	// read topography grid
	topo=topography(&par);
	_topo=topo;

	fprintf(stderr,"read velocity model\n");
	// read velocity file
	vel=readVelocity(par.gridparamFile,par.velocityFile);

	fprintf(stderr,"create mass grid\n");
	// create volumic mass grid with same grid step as velocity model
	vmass=createVMassGrid(&par,vel,topo);
	_vmass=vmass;

	fprintf(stderr,"adjust receiver position\n");
	// adjust receiver positionq to the density and topography grid
	adjustReceiverPosition(rec,topo,&par);

        fprintf(stderr,"compute integration grid\n");
        // check integration grid that must be < velocity grid
        checkGrids(vmass,topo, &prism);

	fprintf(stderr,"dump volumic mass slice\n");
	// create volumic mass grid with same grid step as velocity model
	dumpVMassSlice(&par,vmass);

// loop for integration		
	fprintf(stderr,"start loop\n");
	fflush(stderr);
	mxyz=vmass->nx*vmass->ny*vmass->nz;
	G1= (double*)malloc(sizeof(double)*rec->nr);
	if (G1==NULL)
	{
		fprintf(stderr,"could not allocate memory for G1 (%.0f Kbytes)\n",(double)sizeof(double)*rec->nr/1000);
		exit(-1);
	}
	memset(G1,0,sizeof(double)*rec->nr);
	G2= (double*)malloc(sizeof(double)*rec->nr*mxyz);
	if (G2==NULL)
	{
		fprintf(stderr,"could not allocate memory for G2 (%.0f Kbytes)\n",(double)sizeof(double)*rec->nr*mxyz/1000);
		exit(-1);
	}
	memset(G2,0,sizeof(double)*rec->nr*mxyz);

//use fortran to achieve better openmp parallelisation 
	loop(&prism,rec->x,rec->y,rec->z,rec->nr,mxyz,G1,G2);

// output results
// G: kg m^-2
//mag nT
printf(" DELTA T T | nT | G\n");			//mag
	sprintf(file,"magfull.%s.res",par.prefix);	//mag
	fp=fopen(file,"w");

	for (i=0;i<rec->nr;i++)
	{
		G=G1[i];

	// add contribution from homogeneous and heterogeneous density
		for (j=0;j<mxyz;j++)
		{
			G+= G2[j+i*mxyz]*((double)vmass->vm[j]);
		}

	// output result
		printf("DELTA T %f T %f nT %f G\n",G*Cm,G*Cm*1e9,G);		//mag
		fprintf(fp,"%f %f %f %f\n",rec->x[i],rec->y[i],rec->z[i],G*Cm*1e9);	//mag
	}
	fprintf(stderr,"end\n");
	fclose(fp);

// ouput derivatives
	outputDerivatives(G1,G2,rec->nr,mxyz,vel,vmass,par.prefix);

return 0;
}

