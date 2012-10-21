/* 
 * File:   integrate.h
 * Author: fredrik
 *
 * Created on 19. oktober 2012, 11:20
 */
#include <cstdlib>

/*Integration program. Solves a 6 dimensional integral by Gauss-Legendre quadrature
  and by Monte Carlo integration*/
#include <armadillo>
#include <cmath>
#include <iostream>
#include <time.h>
using namespace std;
using namespace arma;

#define EPS 3.0e-14
#define MAXIT 10

#ifndef INTEGRATE_H
#define	INTEGRATE_H

void gauleg(double x1, double x2, vec &x, vec &w, int n);
double f(double x1,double x2,double y1,double y2, double z1, double z2);
double gammln( double xx);
double timediff(double time1, double time2);
double f_sub(double r1,double r2,double theta1,double theta2,double phi1,double phi2);
void gauss_laguerre(vec &x, vec &w, int n, double alpha);
double ran0(long *idum);
#endif	/* INTEGRATE_H */



void gauleg(double x1, double x2, vec &x, vec &w, int n)
{
   int         m,j,i,s;
   double      z1,z,xm,xl,pp,p3,p2,p1;
   double      const  pi = 3.14159265359; 
   s=n-1;
   m  = (n + 1)/2;
   xm = 0.5 * (x2 + x1);
   xl = 0.5 * (x2 - x1);

   for(i = 1; i <= m; i++) {                             // loops over desired roots
      z = cos(pi * (i - 0.25)/(n + 0.5));

           /*
	   ** Starting with the above approximation to the i'th root
           ** we enter the main loop of refinement by Newtons method.
           */

      do {
         p1 =1.0;
	 p2 =0.0;

   	   /*
	   ** loop up recurrence relation to get the
           ** Legendre polynomial evaluated at x
           */

	 for(j = 1; j <= n; j++) {
	    p3 = p2;
	    p2 = p1;
	    p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3)/j;
	 }

	   /*
	   ** p1 is now the desired Legrendre polynomial. Next compute
           ** ppp its derivative by standard relation involving also p2,
           ** polynomial of one lower order.
           */
 
	 pp = n * (z * p1 - p2)/(z * z - 1.0);
	 z1 = z;
	 z  = z1 - p1/pp;                   // Newton's method
      } while(fabs(z - z1) > 1e-10);

          /* 
	  ** Scale the root to the desired interval and put in its symmetric
          ** counterpart. Compute the weight and its symmetric counterpart
          */
      x(i-1) = xm - xl * z;
      x(s) = xm + xl * z;
      w(i-1) = 2.0 * xl/((1.0 - z * z) * pp * pp);
      w(s) = w(i-1);
      s--;
      
   }
} // End_ function gauleg()
double f(double x1,double x2,double y1,double y2, double z1, double z2){
    double denom = (fabs(sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2))));
    double fexp = -4*(sqrt(x1*x1+y1*y1+z1*z1)+sqrt(x2*x2+y2*y2+z2*z2));
    double func =  exp(fexp)/denom;
    return (denom>1e-12)?func:0;
}
double timediff(double time1, double time2){
    // This function returns the elapsed time in milliseconds
    return ((time2 - time1)*1000)/CLOCKS_PER_SEC;
}
double f_sub(double r1,double r2,double theta1,double theta2,double phi1,double phi2){
    double denom = (r1*r1 + r2*r2 -\
r1*r2*(cos(theta1+theta2)*(1-cos(phi1-phi2))+cos(theta1-theta2)*(1+cos(phi1-phi2))));
    //double pow = -1*(r1+r2);
    double kapow = sin(theta1)*sin(theta2)/sqrt(denom);
    return (denom>1e-12)?kapow:0;
}
void gauss_laguerre(vec &x, vec &w, int n, double alpha)
{
	int i,its,j;
	double ai;
	double p1,p2,p3,pp,z,z1;

	for (i=1;i<=n;i++) {
		if (i == 1) {
			z=(1.0+alpha)*(3.0+0.92*alpha)/(1.0+2.4*n+1.8*alpha);
		} else if (i == 2) {
			z += (15.0+6.25*alpha)/(1.0+0.9*alpha+2.5*n);
		} else {
			ai=i-2;
			z += ((1.0+2.55*ai)/(1.9*ai)+1.26*ai*alpha/
				(1.0+3.5*ai))*(z-x(i-2))/(1.0+0.3*alpha);
                        
		}
		for (its=1;its<=MAXIT;its++) {
			p1=1.0;
			p2=0.0;
			for (j=1;j<=n;j++) {
				p3=p2;
				p2=p1;
				p1=((2*j-1+alpha-z)*p2-(j-1+alpha)*p3)/j;
			}
			pp=(n*p1-(n+alpha)*p2)/z;
			z1=z;
			z=z1-p1/pp;
			if (fabs(z-z1) <= EPS) break;
		}
		if (its > MAXIT) cout << "too many iterations in gaulag" << endl;
		x(i)=z;
		w(i) = -exp(gammln(alpha+n)-gammln((double)n))/(pp*n*p2);
	}
}
// end function gaulag

double gammln( double xx)
{
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
	int j;

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define MASK 123459876

double ran0(long *idum)
{
   long     k;
   double   ans;

   *idum ^= MASK;
   k = (*idum)/IQ;
   *idum = IA*(*idum - k*IQ) - IR*k;
   if(*idum < 0) *idum += IM;
   ans=AM*(*idum);
   *idum ^= MASK;
   return ans;
}
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef MASK
