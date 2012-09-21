/* 
 * File:   main.cpp
 * Author: Fredrik E Pettersen
 * Program description: This program is for project 1 in FYS4150. 
 * The aim is to solve a second order differential equation by linear algebra.
 * Created on 30. august 2012, 08:31
 */

#include <cstdlib>
#include <iostream>
#include <fstream>
#include "armadillo"
#include <time.h>
// #include <functions.h>


using namespace std;
using namespace arma;
/*
 Program for solving a second order differential equation -u''(x) = f(x) 
 */
double frhs(double x){
    //This function defines the right hand side of the diff.eq.
    return 100*exp(-10*x);
    //return 1.8;
}
double timediff(double time1, double time2){
    // This function returns the elapsed time in milliseconds
    return ((time2 - time1)*1000)/CLOCKS_PER_SEC;
}
double uexact(double x){
    //This function defines the exact solution, should it exist
    return 1.0-(1.0-exp(-10.0))*x-exp(-10.0*x);
}
int main(int argc, char** argv) {
    
//------ Initialization--------
    
    int i,n;
    double h,stop,stop_lu, start, start_lu;
    double a = -1.0;
    double b = 2.0;
    double c = -1.0;
    double *avec,*bvec,*cvec,*f,*v,*u, *temp; 
    
    cout << "Enter the size of the NxN matrix" << endl;
    cin >> n;
    mat A(n,n);
    h = 1./(n+1);
    //avec = new double[n];
    //bvec = new double[n];
    //cvec = new double[n];
    u = new double[n];
    f = new double[n];
    v = new double[n];
    double H = h*h; //saves n-1 flops while creating f
    
    for(i=0;i<n;i++){   
        //This loop will initialize arrays
        //avec[i] = a;
        //bvec[i] = b;
        //cvec[i] = c;
        f[i] = frhs((i+1)*h)*H; //Scale f by h**2 to fit the equation
        v[i] = H;               //setting v[i] to h**2 makes no difference
        u[i] = uexact((i+1)*h); //makes the exact solution
    }
    
    for(i=0;i<n;i++){
        //This loop will make the tridiagonal matrix A
        for(int j=0;j<n;j++){
            if(i==j){
                A(i,j) = 2.0;
            }
            else if(abs(i-j)==1){
                A(i,j) = -1.0;
            }
            else{
                A(i,j) = 0;
            }
        }
    }
    
//------ Solving the equations ----------
    
    double btemp;
    temp = new double[n];
    for(i=0;i<n;i++){temp[i]=0;} //set all elements in temp to 0 to avoid mess
    //btemp = bvec[1];  //For the "vetorized" solution
    btemp = 2;  //could use vectors also
    v[0] = f[0]/btemp;
    
    start = clock();    //start timing the computation
    /*
    for(i=1;i<n;i++){
        //forward substitution with vectors
        temp[i] = cvec[i-1]/btemp;
        btemp = bvec[i] - avec[i]*temp[i];
        v[i] = (f[i] - avec[i]*v[i-1])/btemp;
    }
    */
    for(i=1;i<n;i++){
        //forward substitution without vectors
        temp[i] = c/btemp;
        btemp = b +temp[i]; //I've done the multiplication with a allready
        v[i] = (f[i] +v[i-1])/btemp;
    }
    
    for(i=n-2;i>=0;i--){
        //Backward substitution 
        v[i] -= temp[i+1]*v[i+1];
    }
    
    stop = clock();
    cout<<"Computing time = "<< timediff(start,stop) <<" ms"<<endl;
    
//    ofstream myfile;
//    myfile.open("proj1plot.txt");
//    myfile<<0.0<<" "<<0.0;
//    for(i=0;i<=n;i++){
//        myfile << endl;
//        myfile << v[i]<<" "<< u[i];
//    }
//    myfile.close();
   
    //From here some predefined transformations are used
        
//----- LU decomposition ------
    
    mat L;
    mat U;
    mat P;
    mat y(1,n);
    mat z(1,n);
    double LUtemp;
    start_lu = clock();
    lu(L, U, P, A);   //This step is the LU decomposition
    stop_lu = clock();
    
    //substitution
    
    y[0] = f[0];
    for(i=1;i<n;i++){
        /* This loop calculates y from L*y = f. It does not consider calculations
        with 0. Results are only valid for L matrix in LU decomposition. */
        LUtemp=0;
        for(int j=0;j<i;j++){
            LUtemp += L(i,j)*y[j];
        }
        y[i] = f[i] -LUtemp;
    }
    
    double factor = 0;
    z[n-1] = y[n-1]/U(n-1,n-1);
   
    for(i=n-2;i>=0;i--){
        /*This loop continues the calculations from the previous one.
         Calculates U*z = y, where z = v from earlier in the program.
         Resluts are only valid for U matrix in LU decomposition. */
        LUtemp=0;
        for(int j=(n-1);j>=0;j--){
            if(j==i){
                factor = 1.0/U(i,j);
                break;
            }
            LUtemp += U(i,j)*z[j];
        }
        z[i] =(y[i] -LUtemp)*factor;

        
    }
    cout<<"Time elapsed LU decomp. = "<< timediff(start_lu,stop_lu) <<" ms"<<endl;
    z[0] = 0.5*z[1]; //For some reason this element is left out of the above loop
   
    ofstream myfile;
    myfile.open("proj1plot.txt");
    myfile<<0.0<<" "<<0.0<<" "<<0.0;
    for(i=0;i<=n;i++){
        myfile << endl;
        myfile << v[i]<<" "<< u[i]<< " "<<z[i];
    }
    myfile.close();
   
    return 0;
}