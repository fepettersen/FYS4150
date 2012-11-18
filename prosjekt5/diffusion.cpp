/* 
 * File:   main.cpp
 * Author: fredrik
 *
 * Created on 15. oktober 2012, 09:58
 */

//#include "diffusion.h"
#include <armadillo>
#include <iostream>

using namespace std;
using namespace arma;

#define PI 3.14159


void tridiag(double a, double b, double c, vec &v, vec &f, int n){
    vec temp = zeros<vec>(n+1);
    double btemp = 0;


    for(int i=1;i<n;i++){
        //forward substitution without vectors
        temp[i] = c/btemp;
        btemp = b -a*temp[i];
        v[i] = (f[i] -a*v[i-1])/btemp;
    }
    
    for(int i=n-2;i>=0;i--){
        //Backward substitution 
        v[i] -= temp[i+1]*v[i+1];
    }
    //return v;
}

void make_uprev(vec &uprev, double a, double c, double b1, int n){
    for(int i =1;i< n; i++){
        uprev(i) = a*uprev(i-1) +c*uprev(i+1) -b1*uprev(i);
    }
}
int main(int argc, char** argv){
    //--------Common variables--------------
    int N = atoi(argv[1]);
    
    //########################################
    //--------Forward Euler scheme----------##
    //########################################
    
    int N_x = 100;
    vec u_new = zeros<vec>(N_x+1);
    vec u_prev = u_new;
    double dx = 1.0/(N_x+1);
    double dt = dx*dx/4.0;      //Stability criterion dt <= dx*dx/2
    double dtdx2 = dt/(dx*dx);
    int N_t = 1/dt;
    
    for (int n = 0;n<N_t;n++){
        for(int i=1;i<N_x;i++){
            u_new(i) = dtdx2*(u_prev(i+1)-2*u_prev(i) + u_prev(i-1)) + u_prev(i);
        }
        u_prev = u_new;
        /*write to file for plotting*/
    }
    cout<<"Explicit scheme finished. "<<endl;
    
    //#########################################
    //-------Backward Euler scheme-----------##
    //#########################################
    
    double a = -dtdx2;
    double c = a;
    double b = 2+dtdx2;
    for(int n = 0;n<N_t;n++){
       tridiag(a,b,c, u_new, u_prev,N_x); 
       /*Write to file for plotting*/
    }
    cout<<"Implicit scheme finished"<<endl;
    
    //##########################################
    //------Crank Nicolson scheme-------------##
    //##########################################
    
    vec uprev = zeros<vec>(N_x+1);
    a = c = dtdx2;
    double b1 = 2-2*dtdx2;
    double b2 = 2+2*dtdx2;
    
    for(int n = 0; n < N_t; n++){
       make_uprev(uprev,a,c,b1,N_x); 
       tridiag(a,b2,c, u_new, u_prev,N_x); 
       /*Write to file for plotting*/
    }
    cout<<"Crank Nicolson scheme finished"<<endl;
    return 0;
}        
