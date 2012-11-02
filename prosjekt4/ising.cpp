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
    mat  A = zeros<mat>(5,6);
    mat spinmatrix = zeros<mat>(n+2,n+2);
    double start = clock();
    double max_temp = 4.0;
    double temp_step = 0.5;
    vec averages = zeros<vec>(5);

#pragma omp parallel 
    {
        double l_sum = 0;
        double fx = 0;
        double sum_sigma = 0;
        double E,M;
#pragma omp for
    for(double temp = 1; temp < max_temp; temp += temp_step){
        /*Loop over temperatures*/
        cout<<"temp = "<<temp<<" out of "<<max_temp<<endl;
        E = M = 0;
        spinmatrix = init(1,n);
        for(int j = 0; j<N;j++){
            /*Loop over Monte Carlo cycles*/
            update_ghosts(spinmatrix,n);
            metropolis(n,spinmatrix, E, M, temp);
            averages(0)+=E; averages(1)+=E*E; averages(2)+=fabs(M);
            averages(3) +=M;    averages(4) +=M*M;
       }
       variance = averages(1)/((double)N) -(averages(0)/((double)N))*(averages(0)/((double)N));
       crude_mc = averages(4)/((double)N) -(averages(3)/((double)N))*(averages(3)/((double)N));
      cout<<"average energy "<<averages(0)/((double)N)<<" variance_E "<<variance/N;
      cout<<" average magnetization "<<averages(3)/((double)N)<<" variance_M "<<crude_mc/N<<endl;
    }
#pragma omp critical
    {
        crude_mc += l_sum; 
        g_sigma += sum_sigma;
        //num_cores = omp_get_num_threads();
    }
    }

    double stop = clock();
    double diff = timediff(start,stop)*0.25;

    /*Print some results to terminal*/
    /*
    crude_mc /= ((double) N);
    g_sigma /= ((double) N);
    variance = g_sigma - crude_mc*crude_mc;
    crude_mc *= lambda*lambda*4*pow(PI,4);
    cout<<"Monte Carlo simulation with N = "<<N<<" gives "<<crude_mc<<endl;
    cout<<"The variance is "<<variance<<" and standard deviation "<<sqrt(variance)<<endl;
    cout<<"This took "<<diff<<" ms on "<<num_cores<<" cores"<<endl;
*/
    return 0;
}
