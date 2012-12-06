/*
Headerfile for diffusion porject
*/

#include <cstdlib>
#include <omp.h>
#include <armadillo>
#include <cmath>
#include <iostream>
#include <iomanip>
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
void make_uprev(vec &uprev,vec &u, double a, double c, double b1, int n);
void output(ofstream* outfile, vec &u, int n, int scheme, int N);
void output2D(ofstream* outfile, mat &u, int n, int scheme, int N);
void initial_condition(mat &u, double dx, int n);
void update_boundarys(mat &u, double t, double dx, int n);
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

void make_uprev(vec &uprev, vec &u, double a, double c, double b, int n){
    uprev(0) =b*u(0) +c*u(1);
    for(int i =1;i< n; i++){
        uprev(i) = a*u(i-1) +c*u(i+1) + b*u(i);
    }
    uprev(n) = a*u(n-1)+b*u(n);
}

void output(ofstream* outfile, vec &u, int n, int scheme, int N){
    /*outfile is an ofstram-object letting us open a file
    **u is an armadillo-object containing the solution at time n
    **n is the timestep number
    **scheme is an integer telling what scheme is used to obtain the solution
    **N is the size of the array*/
	outfile->open(make_filename(n,scheme));
	for(int i=0;i<=N;i++){
		*outfile <<u(i)<<setprecision(12)<<endl;
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
            *outfile <<u(i,j)<<setprecision(12)<<"  ";
            }
        *outfile <<endl;
    }
    outfile->close();
}

void initial_condition(mat &u, double dx, int n){
    for(int i=0;i<=n;i++){
        for(int j=0;j<=n;j++){
            u(i,j) = (1- j*dx)*exp(i*dx);
        }
    }
}

void update_boundarys(mat &u, double t, double dx, int n){
    for(int i=0; i<=n;i++){
        for(int j=0; j<=n;j++){
            u(0,j) = (1- j*dx)*exp(t);
            u(n,j) = (1- j*dx)*exp(1+t);
            u(i,0) = exp(i*dx+t);
            u(i,n) = 0;
        }
    }
}