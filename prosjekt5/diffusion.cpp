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
int main(int argc, char** argv){
    //--------Common variables--------------
    N = atoi(argv[1]);
    
    //--------Forward Euler scheme----------
    vec u_new = zeros<vec>(N+1);
    vec u_prev = u_new;
    double dx = 1.0/(N+1);
    double dt = dx*dx/4.0;      //Stability criterion dt <= dx*dx/2
    double dtdx2 = dt/(dx*dx);
    for (int n = 0;n<N_t;n++){
        for(int i=1;i<N_x;i++){
            u_new(i) = dtdx2*(u_prev(i+1)-2*u_prev(i) + u_prev(i-1)) + u_prev(i);
        }
        u_prev = u_new;
        /*write to file for plotting*/
    }
    cout<<"Explicit scheme finished. ";
    
    return 0;
}        


