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
    vec bv = zeros<vec>(n+1);
    double temp = 0;
    bv(1)=b;
    v(1) = f(1)/b;
    //v(0) = 1;

    for(int i=2;i<=n-1;i++){
        //forward substitution withoutvectors
        temp = a/bv(i-1);
        bv(i) = b -c*temp;
        f(i) -= f(i-1)*temp;
    }
    v(n-1)= f(n-1)/bv(n-2);
    //v(n)=0;

    for(int i=n-2;i>=1;i--){
        //Backward substitution 
        v[i] = (f(i)-c*v(i+1))/bv(i);
    }
}

void make_uprev(vec &uprev, vec &u, double a, double c, double b, int n){
    uprev(1) =b*u(1) +c*u(2);
    for(int i =2;i< n-1; i++){
        uprev(i) = a*u(i-1) +c*u(i+1) + b*u(i);
    }
    uprev(n-1) = a*u(n-2)+b*u(n-1);
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
        if(i<N){*outfile <<endl;}
    }
    outfile->close();
}

void initial_condition(mat &u, double dx, int n){
    for(int i=0;i<=n;i++){
        for(int j=0;j<=n;j++){
            u(j,i) = (1- j*dx)*exp(i*dx);
        }
    }
}

void update_boundarys(mat &u, double t, double dx, int n){
    for(int i=0; i<=n;i++){
        for(int j=0; j<=n;j++){
            u(j,0) = (1- j*dx)*exp(t);
            u(j,n) = (1- j*dx)*exp(1+t);
            u(0,i) = exp(i*dx+t);
            u(n,i) = 0;
        }
    }
}