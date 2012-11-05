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

mat init(int high,int n,long idum);
double timediff(double time1, double time2);
double ran0(long *idum);
void metropolis(int n,mat spinmatrix,double &E,double &M, vec w, long idum);
void update_ghosts(mat &spinmatrix, int n);
double ran3(long *idum);
double ran2(long *idum);
#endif	/* INTEGRATE_H */



mat init(int high, int n,long idum)
{

	int a,b;
	mat spinmatrix;
	spinmatrix.ones(n+2,n+2);
	if(high){
		//cout<<"hei"<<endl;
		for (int i =0;i<n;i++){
			a = int (1 + n*ran2(&idum));
			b = int (1 + n*ran2(&idum));
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
void metropolis(int n,mat spinmatrix,double &E,double &M, vec w, long idum){
	int a,b,dE;
	for (int x = 1; x <= n; x++){
		for(int y=1; y <= n; y++){
			double r1 = ran2(&idum);
			double r2 = ran2(&idum);
			a = (int) (1 + (n-1)*r1);
			b = (int) (1 + (n-1)*r2);
			//cout<<"a = "<<a<<" b = "<<b<<endl;
			//cout<<"r1 = "<<r1<<" r2 = "<<r2<<endl;
			dE = 2*spinmatrix(a,b)*(spinmatrix(a+1,b)+spinmatrix(a-1,b)+spinmatrix(a,b+1)+spinmatrix(a,b-1));
			cout<<"dE = "<<dE<<" w(dE/4 + 2) "<<w(dE/4 + 2)<<endl;
			if(ran3(&idum) <= w(dE/4 + 2)){
				//
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
#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

double ran3(long *idum)
{
   static int        inext, inextp;
   static long       ma[56];      // value 56 is special, do not modify
   static int        iff = 0;
   long              mj, mk;
   int               i, ii, k;

   if(*idum < 0 || iff == 0) {                 // initialization
      iff    = 1;

      mj     = MSEED - (*idum < 0 ? -*idum : *idum);
      mj    %= MBIG;
      ma[55] = mj;                            // initialize ma[55] 

      for(i = 1, mk = 1; i <= 54; i++) {      // initialize rest of table 
         ii     = (21*i) % 55;
	 ma[ii] = mk;
	 mk     = mj - mk;
	 if(mk < MZ) mk += MBIG;
	 mj = ma[ii];
      }

      for(k = 1; k <= 4; k++) {   // randimize by "warming up" the generator
         for(i = 1; i <= 55; i++) {
	    ma[i] -= ma[1 + (i + 30) % 55];
	    if(ma[i] < MZ) ma[i] += MBIG;
	 }
      }

      inext  =  0;              // prepare indices for first generator number
      inextp = 31;              // 31 is special
      *idum  = 1;
   }

   if(++inext == 56)  inext  = 1;
   if(++inextp == 56) inextp = 1;
   mj = ma[inext] - ma[inextp];
   if(mj < MZ) mj += MBIG;
   ma[inext] = mj;
   return mj*FAC;
}
#undef MBIG
#undef MSEED
#undef MZ
#undef FAC
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double ran2(long *idum)
{
   int            j;
   long           k;
   static long    idum2 = 123456789;
   static long    iy=0;
   static long    iv[NTAB];
   double         temp;

   if(*idum <= 0) {
      if(-(*idum) < 1) *idum = 1;
      else             *idum = -(*idum);
      idum2 = (*idum);
      for(j = NTAB + 7; j >= 0; j--) {
         k     = (*idum)/IQ1;
	 *idum = IA1*(*idum - k*IQ1) - k*IR1;
	 if(*idum < 0) *idum +=  IM1;
	 if(j < NTAB)  iv[j]  = *idum;
      }
      iy=iv[0];
   }
   k     = (*idum)/IQ1;
   *idum = IA1*(*idum - k*IQ1) - k*IR1;
   if(*idum < 0) *idum += IM1;
   k     = idum2/IQ2;
   idum2 = IA2*(idum2 - k*IQ2) - k*IR2;
   if(idum2 < 0) idum2 += IM2;
   j     = iy/NDIV;
   iy    = iv[j] - idum2;
   iv[j] = *idum;
   if(iy < 1) iy += IMM1;
   if((temp = AM*iy) > RNMX) return RNMX;
   else return temp;
}
// End: function ran3()