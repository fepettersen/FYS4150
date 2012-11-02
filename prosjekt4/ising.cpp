/* 
 * File:  ising.cpp
 * Author: fredrik
 *
 * Created on 2. november 2012, 10:30

 * Simulate the development of a system of magnetic dipoles by the Ising model.
 */

#include "ising.h"

#define PI 3.14159

int main(int argc, char** argv) {
  
    int n = atoi(argv[1]);
    int N = 10000;
    double crude_mc,variance;
    crude_mc=0;
    double g_sigma = 0;
    int num_cores = 0;
    spinmatrix = init(0,n);
    double start = clock();
    double max_temp = 4.0;
    double temp_step = 0.5;

#pragma omp parallel 
    {
        double l_sum = 0;
        double r1,r2,theta1,theta2,phi1,phi2;
        double fx = 0;
        double sum_sigma = 0;
        double E,M;
#pragma omp for
    for(int temp = 1; temp < max_temp; temp += temp_step){
        /*Loop over temperatures*/
        E = M = 0;
        spinmatrix = init();
        for(int j = 0; j<N;j++){
            update_ghosts(&spinmatrix,n);
            /*Loop over Monte Carlo cycles*/
             metropolis(n,spinmatrix,&E,&M);
       }
      
    }
#pragma omp critical
    {
        crude_mc += l_sum; 
        g_sigma += sum_sigma;
        num_cores = omp_get_num_threads();
    }
    }

    double stop = clock();
    double diff = timediff(start,stop)*0.25;

    /*Print some results to terminal*/
    crude_mc /= ((double) N);
    g_sigma /= ((double) N);
    variance = g_sigma - crude_mc*crude_mc;
    crude_mc *= lambda*lambda*4*pow(PI,4);
    cout<<"Monte Carlo simulation with N = "<<N<<" gives "<<crude_mc<<endl;
    cout<<"The variance is "<<variance<<" and standard deviation "<<sqrt(variance)<<endl;
    cout<<"This took "<<diff<<" ms on "<<num_cores<<" cores"<<endl;

    return 0;
}
