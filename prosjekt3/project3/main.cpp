/* 
 * File:   main.cpp
 * Author: fredrik
 *
 * Created on 15. oktober 2012, 09:58
 */

#include "integrate.h"

#define PI 3.14159
int main(int argc, char** argv){
    int n =atoi(argv[1]);
    int k = atoi(argv[2]);      /*Read number of cores to use from commandline*/
#if 1
    vec x1(n),x2(n),y1(n),y2(n),z1(n),z2(n),w1(n),w2(n),w3(n),w4(n),w5(n),w6(n);
    x1.zeros();x2.zeros();y1.zeros();y2.zeros();z1.zeros();z2.zeros();
    w1.zeros();w2.zeros();w3.zeros();w4.zeros();w5.zeros();w6.zeros();
    
    double lambda = 3.0;
    gauleg(-lambda,lambda,x1,w1,n);
    gauleg(-lambda,lambda,x2,w2,n);
    gauleg(-lambda,lambda,y1,w3,n);
    gauleg(-lambda,lambda,y2,w4,n);
    gauleg(-lambda,lambda,z1,w5,n);
    gauleg(-lambda,lambda,z2,w6,n);
    int num_cores = 0;
    double sum = 0.0;
    u_long counter=0;
    double start = clock();
#pragma omp parallel shared(x1,x2,y1,y2,z2,z1,n,w1,w2,w3,w4,w5,w6,sum,counter) \
num_threads(k)
    {
        u_long l_counter = 0;
        double l_sum = 0;
#pragma omp for 

    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            for(int k=0;k<n;k++){
                for(int l=0;l<n;l++){
                    for(int m=0;m<n;m++){
                        for(int p=0;p<n;p++){
                            l_sum+=w1(i)*w2(j)*w3(k)*w4(l)*w5(m)*w6(p)*\
                                    f(x1(i),x2(j),y1(k),y2(l),z1(m),z2(p));
                            l_counter++;
                        }
                    }
                }
            }
        }
    }
    
    #pragma omp critical
        {
        sum += l_sum; 
        counter += l_counter;
        num_cores = omp_get_num_threads();
        }
    
    }
    double stop = clock();
    double diff = timediff(start,stop)/((double) k);
    u_long n6 = n*n*n*n*n*n;
    cout<<"sum = "<<sum<<" in "<<diff<<" ms, steps: "<<counter<<" out of "<<n6<<endl;
    cout<<"correct sum  = "<<((5*PI*PI)/(16*16))<<" number of cores used: "<<num_cores<<endl;
 #else 
    vec r(n+1),theta(n),phi(n),w1(n+1),w2(n),w3(n);
    r.zeros();theta.zeros();phi.zeros();
    w1.zeros();w2.zeros();w3.zeros();

    gauleg(0,PI,theta,w2,n);
    gauleg(0,2*PI,phi,w3,n);
    gauss_laguerre(r, w1,n,2.0);

    double sum = 0.0;
    u_long counter=0;
    
    double start = clock();
    
   
#pragma omp parallel shared(r,theta,phi,n,w1,w2,w3,sum,counter)
    {
        u_long l_counter = 0;
        double l_sum = 0;
        //int i,j,k,l,m,p;
#pragma omp for 
    for(int i=0;i<n;i++){
        //cout<<i<<" out of "<<n<<" iterations"<<endl;
        for(int j=0;j<n;j++){
            for(int k=0;k<n;k++){
                for(int l=0;l<n;l++){
                    for(int m=0;m<n;m++){
                        for(int p=0;p<n;p++){
                            l_sum+=w1(i+1)*w1(j+1)*w2(k)*w2(l)*w3(m)*w3(p)*\
                                    f_sub(r(i+1),r(j+1),theta(k),theta(l),phi(m),phi(p));
                            l_counter++;
                        }
                    }
                }
            }
        }
    }
#pragma omp critical
        {
         sum += l_sum; 
         counter += l_counter;
        }
    }
    sum *=1.0/1024;
    double stop = clock();
    double diff = timediff(start,stop);
    u_long n6 = n*n*n*n*n*n;
    cout<<"sum = "<<sum<<" in "<<diff<<" ms, steps: "<<counter<<" out of "<<n6<<endl;
    cout<<"correct sum  = "<<((5*PI*PI)/(16*16))<<endl;
#endif

    return 0;
}        

