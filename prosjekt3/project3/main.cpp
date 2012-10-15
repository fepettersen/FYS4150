/* 
 * File:   main.cpp
 * Author: fredrik
 *
 * Created on 15. oktober 2012, 09:58
 */

#include <cstdlib>

/*Integration program. Solves a 6 dimensional integral by Gauss-Legendre quadrature
  and by Monte Carlo integration*/
#include <armadillo>
#include <cmath>
#include <iostream>
using namespace std;
using namespace arma;

void gauleg(double x1, double x2, vec &x, vec &w, int n)
{
   int         m,j,i,s;
   double      z1,z,xm,xl,pp,p3,p2,p1;
   double      const  pi = 3.14159265359; 
   s=n-1;
   m  = (n + 1)/2;
   xm = 0.5 * (x2 + x1);
   xl = 0.5 * (x2 - x1);

   for(i = 1; i <= m; i++) {                             // loops over desired roots
      z = cos(pi * (i - 0.25)/(n + 0.5));

           /*
	   ** Starting with the above approximation to the i'th root
           ** we enter the main loop of refinement by Newtons method.
           */

      do {
         p1 =1.0;
	 p2 =0.0;

   	   /*
	   ** loop up recurrence relation to get the
           ** Legendre polynomial evaluated at x
           */

	 for(j = 1; j <= n; j++) {
	    p3 = p2;
	    p2 = p1;
	    p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3)/j;
	 }

	   /*
	   ** p1 is now the desired Legrendre polynomial. Next compute
           ** ppp its derivative by standard relation involving also p2,
           ** polynomial of one lower order.
           */
 
	 pp = n * (z * p1 - p2)/(z * z - 1.0);
	 z1 = z;
	 z  = z1 - p1/pp;                   // Newton's method
      } while(fabs(z - z1) > 1e-10);

          /* 
	  ** Scale the root to the desired interval and put in its symmetric
          ** counterpart. Compute the weight and its symmetric counterpart
          */
      x(i-1) = xm - xl * z;
      x(s) = xm + xl * z;
      w(i-1) = 2.0 * xl/((1.0 - z * z) * pp * pp);
      w(s) = w(i-1);
      s--;
      
   }
} // End_ function gauleg()
double f(double x1,double x2,double y1,double y2, double z1, double z2){
    double factor = 1.0/(fabs(sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2))));
    return factor*exp(-2*(x1+x2+y1+y2+z1+z2));
}
int main(int argc, char** argv){
    int n = 3;//atoi(argv[1]);
    vec x1(n),x2(n),y1(n),y2(n),z1(n),z2(n),w1(n),w2(n),w3(n),w4(n),w5(n),w6(n);
    x1.zeros();x2.zeros();y1.zeros();y2.zeros();z1.zeros();z2.zeros();
    w1.zeros();w2.zeros();w3.zeros();w4.zeros();w5.zeros();w6.zeros();
    
    double lambda = 1.0;
    gauleg(-lambda,lambda,x1,w1,n);
    gauleg(-lambda,lambda,x2,w2,n);
    gauleg(-lambda,lambda,y1,w3,n);
    gauleg(-lambda,lambda,y2,w4,n);
    gauleg(-lambda,lambda,z1,w5,n);
    gauleg(-lambda,lambda,z2,w6,n);
    w2.print("w2:");
    double sum = 0.0;
    int counter=0;
    
    for (int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            for(int k=0;k<n;k++){
                for(int l=0;l<n;l++){
                    for(int m=0;m<n;m++){
                        for(int p=0;p<n;p++){
                            sum+=w1(i)*w2(j)*w3(k)*w4(l)*w5(m)*w6(p)*f(x1(i),x2(j),y1(k),y2(l),z1(m),z2(p));
                            cout<<"f(x,y,z) = "<<f(x1(i),x2(j),y1(k),y2(l),z1(m),z2(p))<<endl;
                            counter++;
                        }
                    }
                }
            }
        }
    }
    cout<<"sum = "<<sum<<" after "<<counter<<" steps"<<endl;
    return 0;
}        

