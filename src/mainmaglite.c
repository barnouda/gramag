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
//local protopotypes
void loop(Prism *prism,double *xr,double *yr,double *zr,int nr,double *G);


//
// global
//
Par     *_par;
Topo    *_topo;
VMass   *_vmass;
DMagnet *_dmagnet; //mag

int assignLayer=ASSIGNLAYER;
int projection=PROJECTION;
int checkReceivers=CHECKRECEIVERS;


// ===================================================================================

//			MAIN

// ===================================================================================
int main(int argc, char **argv)
{

char	file[132],*inputFile;
int	i;
double	*G;
Prism   prism;
Par	par;
Rec 	*rec;
Topo	*topo;
Vel	*vel;
VMass	*vmass;
DMagnet dmagnet; //mag
FILE    *fp;

	if (argc==1)
		usage_and_exit(argv[0]);

	// read options
	inputFile=NULL;
	getOptions(argc, argv, &assignLayer,&projection,&checkReceivers, &inputFile);
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
	read_dir(&dmagnet); 						//mag
	_dmagnet=&dmagnet; 						//mag

	fprintf(stderr,"read receiver positions\n");
	// read receiver information
	rec=readReceiver(&par);
	G= (double*)malloc(sizeof(double)*rec->nr);
	memset(G,0,sizeof(double)*rec->nr);

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
		
	fprintf(stderr,"start loop\n");

        loop(&prism,rec->x,rec->y,rec->z,rec->nr,G);



// G: kg m^-2
//mag nT
printf(" DELTA T  T | nT | G\n");			//mag
	sprintf(file,"maglite.%s.res",par.prefix);	//mag
	fp=fopen(file,"w");
	for (i=0;i<rec->nr;i++)
	{
		printf("DELTA T  %f T %f nT %f G\n",G[i]*Cm,G[i]*Cm*1e9,G[i]); //mag
		fprintf(fp,"%f %f %f %f\n",rec->x[i],rec->y[i],rec->z[i],G[i]*Cm*1e9); //mag
	}
	fprintf(stderr,"end\n");
	fclose(fp);

return 0;
}
