#include <stdio.h>
#include "declarations.h"

//
// 3D addressing in vmass->vmass array
//
static inline int vmindex(int i, int j, int k, VMass *vm)
{
// fast axis is i, then j, then k
        return(i+j*vm->nx+k*vm->nx*vm->ny);
}
static inline int velindex(int i, int j, int k, Vel *v)
{
// fast axis is k, then j, then i
//      return(k + (j*v->nz) + i*(v->nz*v->ny));
// fast axis is i, then j, then k
        return(i+j*v->nx+k*v->nx*v->ny);
}

void outputDerivatives(double *G1,double *G2,int nr,int mxyz,Vel *vel,VMass *vm,char *prefix)
{
char   file[132];
int    i,j,k,ijk,ir;
double rho;
FILE  *fp;
	
	// density grid contains more nodes than velocity grid
	// This is because density grid introduced new nodes on the boundaries for the transition between
	// surrounding medium and the velocity grid
	// remove these nodes by copying their contribution from G2 to G1
   	for (ir=0;ir<nr;ir++)
	{
		i=0;
		for (j=0;j<vm->ny;j++)
			for (k=0;k<vm->nz;k++)
			{
				ijk=vmindex(i,j,k,vm);
				G1[ir]+= G2[ijk+ir*mxyz]*vm->vm[ijk];
			}
		i=vm->nx-1;
		for (j=0;j<vm->ny;j++)
			for (k=0;k<vm->nz;k++)
			{
				ijk=vmindex(i,j,k,vm);
				G1[ir]+= G2[ijk+ir*mxyz]*vm->vm[ijk];
			}
		j=0;
		for (i=1;i<vm->nx-1;i++)
			for (k=0;k<vm->nz;k++)
			{
				ijk=vmindex(i,j,k,vm);
				G1[ir]+= G2[ijk+ir*mxyz]*vm->vm[ijk];
			}
		j=vm->ny-1;
		for (i=1;i<vm->nx-1;i++)
			for (k=0;k<vm->nz;k++)
			{
				ijk=vmindex(i,j,k,vm);
				G1[ir]+= G2[ijk+ir*mxyz]*vm->vm[ijk];
			}
		k=0;
		for (i=1;i<vm->nx-1;i++)
			for (j=1;j<vm->ny-1;j++)
			{	
				ijk=vmindex(i,j,k,vm);
				G1[ir]+= G2[ijk+ir*mxyz]*vm->vm[ijk];
			}
		k=vm->nz-1;
		for (i=1;i<vm->nx-1;i++)
			for (j=1;j<vm->ny-1;j++)
			{
				ijk=vmindex(i,j,k,vm);
				G1[ir]+= G2[ijk+ir*mxyz]*vm->vm[ijk];
			}
	}
	// output results

	// contribution from medium surrounding the grid
	sprintf(file,"G1.%s.bin",prefix);
	fp=fopen(file,"w");
	fwrite(&nr,sizeof(int),1,fp);
	fwrite(G1,sizeof(double),nr,fp);
	fclose(fp);

	if (vel==NULL)
		return;
	// contribution from the grid
	sprintf(file,"G2.%s.bin",prefix);
	fp=fopen(file,"w");
	fwrite(&vel->nx,sizeof(int),1,fp);
	fwrite(&vel->ny,sizeof(int),1,fp);
	fwrite(&vel->nz,sizeof(int),1,fp);
	fwrite(&nr,sizeof(int),1,fp);
	for (ir=0;ir<nr;ir++)
		for (k=1;k<vm->nz-1;k++)
			for (j=1;j<vm->ny-1;j++)
				for (i=1;i<vm->nx-1;i++)
					fwrite(G2+vmindex(i,j,k,vm)+ir*mxyz,sizeof(double),1,fp);
	fclose(fp);	

	// density deduced from velocities

	sprintf(file,"rho.%s.bin",prefix);
	fp=fopen(file,"w");
	fwrite(&vel->nx,sizeof(int),1,fp);
	fwrite(&vel->ny,sizeof(int),1,fp);
	fwrite(&vel->nz,sizeof(int),1,fp);
	for (k=1;k<vm->nz-1;k++)
		for (j=1;j<vm->ny-1;j++)
			for (i=1;i<vm->nx-1;i++)
			{
				rho=(double)(vm->vm[vmindex(i,j,k,vm)]);
				fwrite(&rho,sizeof(double),1,fp);
			}
	fclose(fp);
	//velocities
	sprintf(file,"vel.%s.bin",prefix);
	fp=fopen(file,"w");
	fwrite(&vel->nx,sizeof(int),1,fp);
	fwrite(&vel->ny,sizeof(int),1,fp);
	fwrite(&vel->nz,sizeof(int),1,fp);
	for (k=0;k<vel->nz;k++)
		for (j=0;j<vel->ny;j++)
			for (i=0;i<vel->nx;i++)
			{
				rho=(double)(vel->vel[velindex(i,j,k,vel)]);
				fwrite(&rho,sizeof(double),1,fp);
			}
	fclose(fp);		
}
void outputDerivatives2(double *G1,double *G2,int nr,int mxyz,Vel *vel,VMass *vm,char *prefix)
{
char   file[132];
FILE  *fp;
	
	// output results

	// contribution from medium surrounding the grid
	sprintf(file,"G1.%s.bin",prefix);
	fp=fopen(file,"w");
	fwrite(&nr,sizeof(int),1,fp);
	fwrite(G1,sizeof(double),nr,fp);
	fclose(fp);

	if (vel==NULL)
		return;
	// contribution from the grid
	sprintf(file,"G2.%s.bin",prefix);
	fp=fopen(file,"w");
	fwrite(&vm->nx,sizeof(int),1,fp);
	fwrite(&vm->ny,sizeof(int),1,fp);
	fwrite(&vm->nz,sizeof(int),1,fp);
	fwrite(&nr,sizeof(int),1,fp);
	fwrite(G2,sizeof(double),nr*vm->nx*vm->ny*vm->nz,fp);
	fclose(fp);	

	// density deduced from velocities

	sprintf(file,"rho.%s.bin",prefix);
	fp=fopen(file,"w");
	fwrite(&vm->nx,sizeof(int),1,fp);
	fwrite(&vm->ny,sizeof(int),1,fp);
	fwrite(&vm->nz,sizeof(int),1,fp);
	fwrite(vm->vm,sizeof(float),vm->nx*vm->ny*vm->nz,fp);
	fclose(fp);
	//velocities
	sprintf(file,"vel.%s.bin",prefix);
	fp=fopen(file,"w");
	fwrite(&vel->nx,sizeof(int),1,fp);
	fwrite(&vel->ny,sizeof(int),1,fp);
	fwrite(&vel->nz,sizeof(int),1,fp);
	fwrite(vel->vel,sizeof(float),vel->nx*vel->ny*vel->nz,fp);
	fclose(fp);		
}

