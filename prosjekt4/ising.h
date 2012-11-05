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
void metropolis(int n,mat spinmatrix,double &E,double &M, vec w);
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
void metropolis(int n,mat spinmatrix,double &E,double &M, vec w){
	long idum = time(0);
	int a,b,energy,s,dE;
	double r;
	for (int x = 1; x <= n; x++){
		for(int y=1; y <= n; y++){
			a = int (1 + n*ran0(&idum));
			b = int (1 + n*ran0(&idum));
			s = spinmatrix(a,b);
			energy = spinmatrix(a,b)*(spinmatrix(a+1,b)+spinmatrix(a-1,b)+spinmatrix(a,b+1)+spinmatrix(a,b-1));
			//spinmatrix(a,b) *= -1;
			dE = 2*energy;
			//cout<<"  dE = "<<dE<<"  dksnf  "<<w(dE/4 + 2);
			/*
			if(dE>0){
				r = ran0(&idum);
				if(r > w(dE%4 + 2)){
					spinmatrix(a,b)*=-1;
					cout<<"hei"<<endl;
				}
			}
			*/
			if(ran0(&idum) <= w(dE/4 + 2)){
				cout<<"hei"<<endl;
				spinmatrix(a,b) *= -1;
				E += (double) dE;
				M += (double) 2*spinmatrix(a,b);
			}
			
			//spinmatrix.print("...");
			update_ghosts(spinmatrix,n);
		}
	}
}

double timediff(double time1, double time2){
	// This function returns the elapsed time in milliseconds
	return ((time2 - time1)*1000)/CLOCKS_PER_SEC;
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
