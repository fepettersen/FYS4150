/* 
 * File:   montecarlo.cpp
 * Author: fredrik
 *
 * Created on 21. oktober 2012, 18:30
 */

#include "integrate.h"

#define PI 3.14159

int main(int argc, char** argv) {
    double lambda = 3.0;
    long idum1,idum2,idum3,idum4,idum5,idum6;
    idum1 = 12; idum2 = 83;idum3 = 723;idum4 = 2389; idum5 = 234;idum6 = 935;
    int N = atoi(argv[1]);
    double crude_mc,x1,x2,y1,y2,z1,z2,sum_sigma,fx,variance;
    crude_mc=0;
    sum_sigma = 0;
    
    for(int i=0;i<N;i++){
        x1 = -lambda +2*lambda*ran0(&idum1);
        x2 = -lambda +2*lambda*ran0(&idum2);
        y1 = -lambda +2*lambda*ran0(&idum3);
        y2 = -lambda +2*lambda*ran0(&idum4);
        z1 = -lambda +2*lambda*ran0(&idum5);
        z2 = -lambda +2*lambda*ran0(&idum6);
        //cout<<"x1 = "<<x1<<"  y2 = "<<y2<<"  z1 = "<<z1<<endl;
        fx = f(x1,x2,y1,y2,z1,z2);
        crude_mc += fx;
        sum_sigma += fx*fx;
    }
    cout<<"sum = "<<crude_mc<<endl;
    crude_mc = crude_mc*pow((2*lambda),6)/((double) N);
    sum_sigma = sum_sigma*pow((2*lambda),6)/((double) N);
    variance = sum_sigma - crude_mc*crude_mc;
    cout<<"Monte Carlo simulation with N = "<<N<<" gives "<<crude_mc<<endl;
    cout<<"The variance is "<<variance<<" and standard deviation "<<sqrt(variance)<<endl;
    return 0;
}
