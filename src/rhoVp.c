/*
 * Nafe and Drake, 1963 law between density and VP, taken from Onizawa, Mikada, Watanabe et Sakashita, Earth Planet Space, VOL54, 2002
 */
double NafeAndDrake(double vp)
{
	double rho, delta, a, b, c, A, B, C;
#define NAFE 1
#define LIN1 2
#define LIN2 3
//AB#define LCST 4		//AB
#define LOI LIN2	//AB
//AB#define LOI LCST	//AB
#if LOI == NAFE
	if (vp==0)
		return(0);
	else if (vp<1500.)	//2000
		return(vp);	//1880
	vp/=1000;
	a=6.86;
	b=-7.55;
	c=2.64;
	//vp=a + b*rho + c*rho²
	A=c;
	B=b;
	C=a-vp;
	delta=B*B-4*A*C;
	rho=(-B+sqrt(delta))/2/A*1000;
#elif LOI == LIN1
	if (vp>3500)
		return(2700);
	else if (vp<2000)
		return(1500);
	else
		rho=1500+(2700-1500)*(vp-2000)/(3500-2000);
#elif LOI == LIN2
	if (vp<500)
		rho=1000;
	else if (vp <1100)
		rho=1000*(1100-vp)/(1100-500) + 2175*(vp-500)/(1100-500);
	else if (vp<2900)
		rho=2175*(2900-vp)/(2900-1100) + 2700*(vp-1100)/(2900-1100);
	else if (vp<4000)//max 4000
		rho=2700*(4000-vp)/(4000-2900) + 2800*(vp-2900)/(4000-2900);
	else
		rho=2800;
//AB#elif LOI == LCST	//AB
//AB	rho=1000;	//AB
#endif
//fprintf(stderr,"Birch %f %f\n",vp,rho);
	return(rho);
}
