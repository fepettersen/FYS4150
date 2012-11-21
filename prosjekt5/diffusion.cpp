/* 
 * File:   main.cpp
 * Author: fredrik
 *
 * Created on 15. oktober 2012, 09:58
 */

#include "diffusion.h"
/*
#include <armadillo>
#include <iostream>

using namespace std;
using namespace arma;
*/
#define PI 3.14159


int main(int argc, char** argv){
    //--------Common variables--------------
    //int N = atoi(argv[1]);
    ofstream outfile;
    int tofile = atoi(argv[1]);
    int spacing = atoi(argv[2]);
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
        u_new(0) = 1; u_new(N_x) = 0;
        u_prev = u_new;
        /*write to file for plotting*/
        if(tofile && (n%spacing)==0){output(&outfile,u_prev,n,0,N_x);}

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
