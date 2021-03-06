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
    int N = 10000000;
    //double g_sigma = 0;
    //int num_cores = 0;
    //double start = clock();
    
    ofstream outfile;
   
    double start_temp = 1.0;
    double max_temp = 1.0001;
    double temp_step = 1.0;
    int ntemps = ((max_temp-start_temp)/temp_step)+1;
#pragma omp parallel 
    {
        double average_E = 0;
        double average_M = 0;
        double average_E2 = 0;
        double average_M2 = 0;
        double average_Mabs = 0;
        double E,M;
        double variance_M,variance_E;
        variance_E=variance_M = 0;
        long idum = -1*time(0);
        mat spinmatrix = zeros<mat>(n+2,n+2);
        vec averages = zeros<vec>(5);
        int accepted_flips = 0;
        vec w = zeros<vec>(5);
        vec temp = linspace<vec>(start_temp,max_temp,ntemps);
        double prob[n*n+1];
        for(int i = 0; i<=n*n;i++){
            prob[i] = 0;
        }
  

#pragma omp for
    for(int t = 0; t < ntemps; t++){
        /*Loop over temperatures*/
        averages(0) = 0; averages(1) = 0; averages(2) = 0;
        averages(3) = 0; averages(4) = 0;
        average_E = 0;
        average_M = 0;
        average_E2 = 0;
        average_M2 = 0;
        average_Mabs = 0;
        for(int p=0;p<5;p++){
            w(p) = exp(-(4*p-8)/temp(t));
        }
        E = M = 0;
        spinmatrix = init(0,n,E,M);
        for(int j = 0; j<N;j++){
            /*Loop over Monte Carlo cycles*/
            accepted_flips = metropolis(n, spinmatrix, E, M, w, &idum);
            if (j>0.1*N){
                averages(0) += E; averages(1) += E*E; averages(2) += fabs(M);
                averages(3) += M; averages(4) += M*M;
                int j = (E+2*n*n)/4;
                prob[j] +=1;
                /*
                outfile<<averages(0)/((double)j)<<"               "<< averages(2)/((double)j)<<\
                        "      "<<accepted_flips<<"   "<<probendl;
                  */      
            }
        }
        /*Handling the results*/
        //outfile.open("probability2.txt");
        
        average_E = averages(0)/(0.9*((double) N));
        average_M = averages(3)/(0.9*((double) N));    //Note the use of abs(M)
        average_Mabs = averages(2)/(0.9*((double) N));
        average_E2 = averages(1)/(0.9*((double) N));
        average_M2 = averages(4)/(0.9*((double) N));
        variance_E = (average_E2 - average_E*average_E)/n/n;
        variance_M = (average_M2 - average_M*average_M)/n/n;
        cout<<"-----------------------------"<<endl;
        cout<<"temp = "<<temp(t)<<" out of "<<max_temp<<endl;
        cout<<"average energy "<<average_E/n/n<<" heat capacity "<<variance_E/(temp(t)*temp(t))<<endl;
        cout<<"average magnetization "<<average_Mabs/n/n<<" magnetic suceptibility "<<\
                variance_M/(temp(t))<<endl;
        /*        
        outfile<<average_E/n/n<<"      "<<variance_E/(temp(t)*temp(t))<<\
                "        "<<average_Mabs/n/n <<"           "<<\
                variance_M/(temp(t))<<"      "<<temp(t)<<endl;
        
        for(int i=0;i<=n*n;i++){
            outfile<<prob[i]/((double)0.9*N)<<endl;
        }
        outfile.close();
        */
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
