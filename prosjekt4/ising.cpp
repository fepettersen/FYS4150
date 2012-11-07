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
    int N = 500000;
    double variance_M,variance_E;
    variance_E=variance_M = 0;
    //double g_sigma = 0;
    //int num_cores = 0;
    mat spinmatrix = zeros<mat>(n+2,n+2);
    //double start = clock();
    
    double temp_step = 0.1;
    double start_temp = 1.6;
    double max_temp = 2.0;
    vec averages = zeros<vec>(5);
    long idum = -1*time(0);
    ofstream outfile;
    outfile.open("results.txt");
#pragma omp parallel 
    {
        double average_E = 0;
        double average_M = 0;
        double average_E2 = 0;
        double average_M2 = 0;
        double E,M;
        vec w = zeros<vec>(5);
#pragma omp for
    for(double temp = start_temp; temp < max_temp; temp += temp_step){
        /*Loop over temperatures*/
        averages(0) = 0; averages(1) = 0; averages(2) = 0;
        averages(3) = 0; averages(4) = 0;
        average_E = 0;
        average_M = 0;
        average_E2 = 0;
        average_M2 = 0;
        cout<<"temp = "<<temp<<" out of "<<max_temp<<endl;
        for(int p=0;p<5;p++){
            w(p) = exp(-(4*p-8)/temp);
        }
        E = M = 0;
        spinmatrix = init(1,n,E,M);

        for(int j = 0; j<N;j++){
            /*Loop over Monte Carlo cycles*/
             
            metropolis(n, spinmatrix, E, M, w, idum);
            //cout<<"Initial E = "<<E<<" initial M = "<<M<<endl;
            if (j>N/10){
                averages(0) += E; averages(1) += E*E; averages(2) += fabs(M);
                averages(3) += M; averages(4) += M*M;
            }
        }
        //cout<<"Cumulative energy: "<<averages(0)<<" compared to 4*8*N = "<<4*8*N<<endl;
        average_E = averages(0)/((double) N);
        average_M = averages(2)/((double) N);    //Note the use of abs(M)
        average_E2 = averages(1)/((double) N);
        average_M2 = averages(4)/((double) N);
        outfile<<average_E<<"  "<< average_M<<"  "<< N<<endl;
        variance_E = (average_E2 - average_E*average_E)/n/n;
        variance_M = (average_M2 - average_M*average_M)/n/n;
        cout<<"average energy "<<average_E/n/n<<" variance_E "<<variance_E<<endl;
        cout<<"average magnetization "<<average_M/n/n<<" variance_M "<<variance_M<<endl;
    
    }
#pragma omp critical
    {
        //crude_mc += l_sum; 
        //g_sigma += sum_sigma;
        //num_cores = omp_get_num_threads();
    }
    }
    outfile.close();
    //double stop = clock();
    //double diff = timediff(start,stop)*0.25;

    return 0;
}
