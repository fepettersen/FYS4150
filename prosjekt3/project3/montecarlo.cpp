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
 
    //idum1 = ;
    //cout<<"time: "<<idum1<<endl;
    int N = atoi(argv[1]);
    double crude_mc,variance;
    crude_mc=0;
    //sum_sigma = 0;
    //int k = 4;
    double g_sigma = 0;
    int num_cores = 0;
 #if 1
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
        
        //cout<<" f(x) "<<fx<<endl;
    }
        //cout<<"l_sum "<<l_sum<<endl;
#pragma omp critical
    {crude_mc += l_sum; g_sigma += sum_sigma;num_cores = omp_get_num_threads();}
    }
    double stop = clock();
    double diff = timediff(start,stop)*0.25;
    //cout<<"sum = "<<crude_mc<<endl;
 
    crude_mc /= ((double) N);
    crude_mc *= lambda*lambda*4*pow(PI,4);
    g_sigma /= ((double) N);
    g_sigma *= (lambda*lambda*4*pow(PI,4));     //*(lambda*lambda*4*pow(PI,4))
    variance = g_sigma - crude_mc*crude_mc;
    cout<<"Monte Carlo simulation with N = "<<N<<" gives "<<crude_mc<<endl;
    cout<<"The variance is "<<variance<<" and standard deviation "<<sqrt(variance)<<endl;
    cout<<"This took "<<diff<<" ms on "<<num_cores<<" cores"<<endl;
#else
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
        r1 = -log(1-ran0(&idum1));
        r2 = -log(1-ran0(&idum2));
        theta1 = PI*ran0(&idum3);
        theta2 = PI*ran0(&idum4);
        phi1 = 2*PI*ran0(&idum5);
        phi2 = 2*PI*ran0(&idum6);
        //cout<<"x1 = "<<x1<<"  y2 = "<<y2<<"  z1 = "<<z1<<endl;
        fx = r1*r1*r2*r2*f_sub(r1,r2,theta1,theta2,phi1,phi2)*exp(-4*(r1+r2));
        l_sum += fx;
        sum_sigma += fx*fx;
    }
#pragma omp critical
        {}
    }
    cout<<"sum = "<<crude_mc<<" sum_sigma = "<<g_sigma<<endl;
    crude_mc = crude_mc/((double) N);
    crude_mc *= 4*pow(PI,4)*lambda*lambda;
    g_sigma = g_sigma*pow((2*lambda),6)/((double) N);
    cout<<"sum_sigma = "<<g_sigma<<endl;
    variance = g_sigma - crude_mc*crude_mc;
    cout<<"Monte Carlo simulation with N = "<<N<<" gives "<<crude_mc<<endl;
    cout<<"The variance is "<<variance<<" and standard deviation "<<sqrt(variance)<<endl;
#endif
    return 0;
}
