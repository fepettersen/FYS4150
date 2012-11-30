/*
Headerfile for diffusion porject
*/

#include <cstdlib>
#include <omp.h>
#include <armadillo>
#include <cmath>
#include <iostream>
#include <time.h>
#include <fstream>

/*from lib.h*/
#include <new>
#include <cstdio>
#include <cstdlib>
#include <cstring>


using namespace std;
using namespace arma;


#ifndef DIFFUSION_H
#define	DIFFUSION_H

double timediff(double time1, double time2);
char *make_filename(int n, int scheme);
void tridiag(double a, double b, double c, vec &v, vec &f, int n);
void make_uprev(vec &uprev, double a, double c, double b1, int n);
void output(ofstream* outfile, vec &u, int n, int scheme, int N);
void output2D(ofstream* outfile, mat &u, int n, int scheme, int N);
#endif	/* DIFFUSION_H */


double timediff(double time1, double time2){
	// This function returns the elapsed time in milliseconds
	return ((time2 - time1)*1000)/CLOCKS_PER_SEC;
}
char *make_filename(int n, int scheme){
        //Returns a filename saying something about the particular run.
        char* buffer = new char[60];
		if(scheme == 0){
        	sprintf(buffer,"results_FE_n%d.txt",n);
    	}
    	else if(scheme ==1){
    		sprintf(buffer,"results_BE_n%d.txt",n);	
    	}
    	else if(scheme ==2){
    		sprintf(buffer,"results_CN_n%d.txt",n);
        }
        else if(scheme ==3){
            sprintf(buffer,"results_FE2D_n%d.txt",n);
        }
        else{
            sprintf(buffer,"results_LF2D_n%d.txt",n);
        }
        return buffer;
}


void tridiag(double a, double b, double c, vec &v, vec &f, int n){
    vec temp = zeros<vec>(n+1);
    double btemp = b;
    v(0) = f(0)/btemp;
    //v(0) = 1;

    for(int i=1;i<n;i++){
        //forward substitution withoutvectors
        temp[i] = c/btemp;
        btemp = b -a*temp[i];
        v[i] = (f[i] -a*v[i-1])/btemp;
    }
    
    for(int i=n-2;i>=0;i--){
        //Backward substitution 
        v[i] -= temp[i+1]*v[i+1];
    }
    //v.print("asdf");
    //return v;
}

void make_uprev(vec &uprev, double a, double c, double b1, int n){
    uprev(0) =b1*uprev(0) +c*uprev(1);
    for(int i =1;i< n; i++){
        uprev(i) = a*uprev(i-1) +c*uprev(i+1) - b1*uprev(i);
    }
    uprev(n) = a*uprev(n-2)+b1*uprev(n-1);
}

void output(ofstream* outfile, vec &u, int n, int scheme, int N){
    /*outfile is an ofstram-object letting us open a file
    **u is an armadillo-object containing the solution at time n
    **n is the timestep number
    **scheme is an integer telling what scheme is used to obtain the solution
    **N is the size of the array*/
	outfile->open(make_filename(n,scheme));
	for(int i=0;i<N;i++){
		*outfile <<u(i)<<endl;
	}
	outfile->close();
}

void output2D(ofstream* outfile, mat &u, int n, int scheme, int N){
    /*outfile is an ofstram-object letting us open a file
    **u is an armadillo-object containing the solution at time n
    **n is the timestep number
    **scheme is an integer telling what scheme is used to obtain the solution
    **N is the size of the array (in one direction)*/
    outfile->open(make_filename(n,scheme));
    for(int i=0;i<=N;i++){
        for(int j=0;j<=N;j++){
            *outfile <<u(i,j)<<"  ";
            }
        *outfile <<endl;
    }
    outfile->close();
}