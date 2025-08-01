/*
* linear interpolation in 1(?), 2 and 3 dimension
*/

double interp1(double x,double *X,double *V)
{
#define NNODE1 2
int      i,j;
int      nd=NNODE1;    //node number
double    v,vi,num,den;
	v=0;
	for (i=0;i<nd;i++)
	{
	   vi=V[i];
	   for (j=0;j<nd;j++)
	   {
	       if (j==i)
	           continue;
	       num=(x-X[j]);
 	      den=(X[i]-X[j]);
 	      vi*=num/den;
	   }
	   v+=vi;
	}
   	return(v);
}
/*
 *      0 -- 3
 *      |    |     ^
 *      1 -- 2     |_ 
 */

double interp2(double x,double y,double *X,double *Y,double *V)
{
#define NNODE2 4
double   v,num,den,numden;

	v=0;
	num=(x-X[2])*(y-Y[1]);
	den=(X[0]-X[2])*(Y[0]-Y[1]);
	numden=num/den;
	v+=V[0]*numden;
	num=(x-X[2])*(y-Y[0]);
	den=(X[1]-X[2])*(Y[1]-Y[0]);
	numden=num/den;
	v+=V[1]*numden;
	num=(x-X[0])*(y-Y[0]);
	den=(X[2]-X[0])*(Y[2]-Y[0]);
	numden=num/den;
	v+=V[2]*numden;
	num=(x-X[0])*(y-Y[1]);
	den=(X[3]-X[0])*(Y[3]-Y[1]);
	numden=num/den;
	v+=V[3]*numden;

	return(v);
}
/*
 * 3D linear interpolation in a rectangular prism
 */
double interp3(double x,double y,double z,double *X,double *Y,double *Z,double *V)
{
#define NNODE3 8
int      i,j;
int      nd=NNODE3;    //node number
double    v,num,den;
double	xmin,xmax,ymin,ymax,zmin,zmax,xc[NNODE3],yc[NNODE3],zc[NNODE3];

	xmax=xmin=X[0];
	ymin=ymax=Y[0];
	zmin=zmax=Z[0];
	for (i=1;i<nd;i++)
	{
		xmin=(xmin>X[i])? X[i]: xmin;
		ymin=(ymin>Y[i])? Y[i]: ymin;
		zmin=(zmin>Z[i])? Z[i]: zmin;
		xmax=(xmax<X[i])? X[i]: xmax;
		ymax=(ymax<Y[i])? Y[i]: ymax;
		zmax=(zmax<Z[i])? Z[i]: zmax;
	}
	for (i=0;i<nd;i++)
	{
		xc[i]=(X[i]==xmin)? xmax: xmin;
		yc[i]=(Y[i]==ymin)? ymax: ymin;
		zc[i]=(Z[i]==zmin)? zmax: zmin;
	}
	v=0;
	for (i=0;i<nd;i++)
	{
	   num=(x-xc[i])*(y-yc[i])*(z-zc[i]);
	   den=(X[i]-xc[i])*(Y[i]-yc[i])*(Z[i]-zc[i]);
	   v+=V[i]*num/den;
	}
   	return(v);
}
/*
 * 3D linear interpolation in a rectangular prism
 * same as above but faster coding
 */
double interp3b(double x,double y,double z,double *X,double *Y,double *Z,double *V,double *alpha)
{
double    v,num,den;

	v=0;
//nodes
//0
	den=(X[0]-X[4])*(Y[0]-Y[2])*(Z[0]-Z[1]);
	num=(x-X[4])*(y-Y[2])*(z-Z[1]);
	alpha[0]=num/den;
	v+=V[0]*alpha[0];
//1
	den=(X[1]-X[4])*(Y[1]-Y[2])*(Z[1]-Z[0]);
	num=(x-X[4])*(y-Y[2])*(z-Z[0]);
	alpha[1]=num/den;
	v+=V[1]*alpha[1];
//2
	den=(X[2]-X[4])*(Y[2]-Y[0])*(Z[2]-Z[1]);
	num=(x-X[4])*(y-Y[0])*(z-Z[1]);
	alpha[2]=num/den;
	v+=V[2]*alpha[2];
//3
	den=(X[3]-X[4])*(Y[3]-Y[0])*(Z[3]-Z[0]);
	num=(x-X[4])*(y-Y[0])*(z-Z[0]);
	alpha[3]=num/den;
	v+=V[3]*alpha[3];
//4
	den=(X[4]-X[0])*(Y[4]-Y[2])*(Z[4]-Z[1]);
	num=(x-X[0])*(y-Y[2])*(z-Z[1]);
	alpha[4]=num/den;
	v+=V[4]*alpha[4];
//5
	den=(X[5]-X[0])*(Y[5]-Y[2])*(Z[5]-Z[0]);
	num=(x-X[0])*(y-Y[2])*(z-Z[0]);
	alpha[5]=num/den;
	v+=V[5]*alpha[5];
//6
	den=(X[6]-X[0])*(Y[6]-Y[0])*(Z[6]-Z[1]);
	num=(x-X[0])*(y-Y[0])*(z-Z[1]);
	alpha[6]=num/den;
	v+=V[6]*alpha[6];
//7
	den=(X[7]-X[0])*(Y[7]-Y[0])*(Z[7]-Z[0]);
	num=(x-X[0])*(y-Y[0])*(z-Z[0]);
	alpha[7]=num/den;
	v+=V[7]*alpha[7];

   	return(v);
}
