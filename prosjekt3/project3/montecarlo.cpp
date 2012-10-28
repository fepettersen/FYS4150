/* 
 * File:   montecarlo.cpp
 * Author: fredrik
 *
 * Created on 21. oktober 2012, 18:30
 */

#include "integrate.h"

#define PI 3.14159

int main(int argc, char** argv) {
    double lambda = atof(argv[2]);
 
    int N = atoi(argv[1]);
    double crude_mc,variance;
    crude_mc=0;
    double g_sigma = 0;
    int num_cores = 0;
 #if 0
    double start = clock();
#pragma omp parallel 
    {
        double l_sum = 0;
        double r1,r2,theta1,theta2,phi1,phi2;
        double fx = 0;
        double sum_sigma = 0;
        long idum1,idum2,idum3,idum4,idum5,idum6;
        idum1 = time(0)-123; idum2 = time(0)-383;idum3 = time(0)-73;idum4 = time(0)-239; 
        idum5 = time(0)-830;idum6 = time(0)-955;
#pragma omp for
    for(int i=0;i<N;i++){
        r1 = lambda*ran0(&idum1);
        r2 = lambda*ran0(&idum2);
        theta1 = PI*ran0(&idum3);
        theta2 = PI*ran0(&idum4);
        phi1 = 2*PI*ran0(&idum5);
        phi2 = 2*PI*ran0(&idum6);
       
        fx = r1*r1*r2*r2*f_sub(r1,r2,theta1,theta2,phi1,phi2)*exp(-4*(r1+r2));
        l_sum += fx;
        sum_sigma += fx*fx;
    }
#pragma omp critical
    {crude_mc += l_sum; g_sigma += sum_sigma;num_cores = omp_get_num_threads();}
    }
    double stop = clock();
    double diff = timediff(start,stop)*0.25;
    /*
    double var = ((pow(PI,4)*lambda*lambda*4)/N)*(g_sigma-((pow(PI,4)*lambda*lambda*4)/N)*crude_mc*crude_mc);
    cout<<"test: "<<var<<" sdv: "<<sqrt(var)<<" comparisson: "<<var/sqrt(N)<<endl;*/
    crude_mc /= ((double) N);
    g_sigma /= ((double) N);
    variance = g_sigma - crude_mc*crude_mc;
    crude_mc *= lambda*lambda*4*pow(PI,4);
    cout<<"Monte Carlo simulation with N = "<<N<<" gives "<<crude_mc<<endl;
    cout<<"The variance is "<<variance<<" and standard deviation "<<sqrt(variance)<<endl;
    cout<<"This took "<<diff<<" ms on "<<num_cores<<" cores"<<endl;
#else
    double start = clock();
#pragma omp parallel 
    {
        /*Set up private variables for the paralell execution of this loop*/
        double l_sum = 0;
        double l_sum2 = 0;
        double r1,r2,theta1,theta2,phi1,phi2;
        double fx = 0;
        double sum_sigma = 0;
        double sum_sigma2 = 0;
        long idum1,idum2,idum3,idum4,idum5,idum6;
        idum1 = time(0)-123; idum2 = time(0)-383;idum3 = time(0)-73;idum4 = time(0)-239; 
        idum5 = time(0)-830;idum6 = time(0)-955;
#pragma omp for
    
    for(int i=0;i<N;i++){
        r1 = -log(1-ran0(&idum1));
        r2 = -log(1-ran0(&idum2));
        theta1 = PI*ran0(&idum3);
        theta2 = PI*ran0(&idum4);
        phi1 = 2*PI*ran0(&idum5);
        phi2 = 2*PI*ran0(&idum6);
        fx = (1./1024)*r1*r1*r2*r2*f_sub(r1,r2,theta1,theta2,phi1,phi2);
        if (l_sum > 1e4){
            l_sum2 +=fx;
            sum_sigma2 +=fx*fx;
        }
        else{
                l_sum += fx;
                sum_sigma += fx*fx;
        }
    }
#pragma omp critical
        {
        crude_mc +=l_sum+l_sum2; 
        g_sigma += sum_sigma+sum_sigma2; 
        num_cores = omp_get_num_threads();
        }
    }
    double stop = clock();
    double diff = timediff(start,stop)*0.25;
    crude_mc /= ((double) N);
    g_sigma /= ((double) N);
    variance = g_sigma - crude_mc*crude_mc;
    crude_mc *= 4*pow(PI,4);
    cout<<"Monte Carlo simulation with (importance sampling) N = "<<N<<" gives "<<crude_mc<<endl;
    cout<<"correct result: "<<5*PI*PI/(16*16)<<endl;
    cout<<"The variance is "<<variance<<" and standard deviation "<<pow(4,PI)*sqrt(variance/((double) N))<<endl;
    cout<<"This took "<<diff<<" ms on "<<num_cores<<" cores"<<endl;
#endif
    return 0;
}
