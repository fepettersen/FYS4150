/* 
 * File:   integrate.h
 * Author: fredrik
 *
 * Created on 2. november 2012, 10:30
 */

#include <cstdlib>
#include <omp.h>
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

mat init(int high,int n);
//double f(double x1,double x2,double y1,double y2, double z1, double z2);
double timediff(double time1, double time2);
//double f_sub(double r1,double r2,double theta1,double theta2,double phi1,double phi2);
double ran0(long *idum);
void metropolis(int n,mat spinmatrix,double &E,double &M, double temp);
void update_ghosts(mat &spinmatrix, int n);
#endif	/* INTEGRATE_H */



mat init(int high, int n)
{

	int a,b;
	long idum = time(0);
	mat spinmatrix;
	spinmatrix.ones(n+2,n+2);
	if(high){
		//cout<<"hei"<<endl;
		for (int i =0;i<n;i++){
			a = int (1 + n*ran0(&idum));
			b = int (1 + n*ran0(&idum));
			//cout<<"a = "<<a<<" b = "<<b<<endl;
			spinmatrix(a,b) *= -1;
		}
  }
  return spinmatrix;
}
void update_ghosts(mat &spinmatrix, int n){
	
	for(int i =1; i <= n; i++){
		spinmatrix(0,i) = spinmatrix(n,i);
		spinmatrix(n+1,i) = spinmatrix(1,i);
		spinmatrix(i,0) = spinmatrix(i,n);
		spinmatrix(i,n+1) = spinmatrix(i,1);
	}
	return ;
}
void metropolis(int n,mat spinmatrix,double &E,double &M, double temp){
	long idum = time(0);
	int a,b,energy,s,dE;
	double r;
	for (int x = 1;x<=n;x++){
		for(int y=1;y<=n;y++){
			//cout<<"hei, y  ="<<y<<endl;
			a = int (1 + n*ran0(&idum));
			b = int (1 + n*ran0(&idum));
			//cout<<"a = "<<a<<" b = "<<b<<endl;
			s = spinmatrix(a,b);
			energy = spinmatrix(a+1,b)+spinmatrix(a-1,b)+spinmatrix(a,b+1)+spinmatrix(a,b-1);
			spinmatrix(a,b) *= -1;
			dE = 2*energy;
			if(dE>0){
				r = ran0(&idum);
				if(r>exp(-dE/temp)){
					spinmatrix(a,b)*=-1;
				}
			}
			E += dE;
			M += spinmatrix(x,y);
		}
	}
}
/*
double f(double x1,int high,double y1,double y2, double z1, double z2){
	double denom = (fabs(sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2))));
	double fexp = -4*(sqrt(x1*x1+y1*y1+z1*z1)+sqrt(x2*x2+y2*y2+z2*z2));
	double func =  exp(fexp)/denom;
	return (denom>1e-12)?func:0;
}
*/
double timediff(double time1, double time2){
	// This function returns the elapsed time in milliseconds
	return ((time2 - time1)*1000)/CLOCKS_PER_SEC;
}
/*
double f_sub(double r1,double r2,double theta1,double theta2,double phi1,double phi2){
	double denom = (r1*r1 + r2*r2 -\
r1*r2*(cos(theta1+theta2)*(1-cos(phi1-phi2))+cos(theta1-theta2)*(1+cos(phi1-phi2))));
	//double pow = -1*(r1+r2);
	double kapow = sin(theta1)*sin(theta2)/sqrt(denom);
	return (denom>1e-12)?kapow:0;
}
*/
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
