/* 
 * File:   diffusion.cpp
 * Author: candidate 55
 *
 * Created on 28. november 2012, 09:58
 */

#include "diffusion.h"


int main(int argc, char** argv){
    //--------Common variables--------------
    //int N = atoi(argv[1]);
    ofstream outfile;
    int tofile = atoi(argv[1]);
    int spacing = atoi(argv[2]);
    int FE1D = atoi(argv[3]);
    int BE1D = atoi(argv[4]);
    int CN1D = atoi(argv[5]);
    int FE2D = atoi(argv[6]);
    int LeapFrog = atoi(argv[7]);
    int nx = atoi(argv[8]);
    int n_t = atoi(argv[9]);

    vec u_new = zeros<vec>(nx+1);
    vec u_prev = u_new;
    double dx = 1.0/(nx);
    double dt = dx*dx/4.0;      //Stability criterion dt <= dx*dx/2
    double dtdx2 = dt/(dx*dx);
    mat U = zeros<mat>(nx+1,nx+1);
    mat U_p = U;
    mat U_pp = U;

    if(FE1D){
    //########################################
    //--------Forward Euler scheme----------##
    //########################################
        if(dtdx2>0.5){dtdx2 = 0.5;} //make sure the stability criterion is fulfilled
        u_prev = linspace<vec>(-1.0,0.0,nx+1);
        cout<<"alpha = "<<dtdx2<<endl;
        for(int n = 0;n <= n_t; n++){
            for(int i = 1; i <nx; i++){
                u_new(i) = dtdx2*(u_prev(i+1)-2*u_prev(i) + u_prev(i-1)) + u_prev(i);
            }
            u_new(0) = 0; u_new(nx) = 0;
            u_prev = u_new;
            /*write to file for plotting*/
            if(tofile && (n%spacing)==0){output(&outfile,u_prev,n,0,nx);cout<<n<<endl;}

        }
        cout<<"Explicit scheme finished. "<<endl;
        }
    if(BE1D){
    //#########################################
    //-------Backward Euler scheme-----------##
    //#########################################
        //dtdx2 = 1/(dx*dx*200);
        u_prev = linspace<vec>(-1.0,0.0,nx+1);
        double a = -dtdx2;
        double c = a;
        double b = 1+2*dtdx2;
        u_new.zeros();
        for(int n = 1;n<=n_t;n++){
            tridiag(a,b,c, u_new, u_prev,nx);
            u_new(0)=0; u_new(nx) = 0;
            //u_prev = u_new;
            for(int k=0;k<=nx;k++){
                u_prev(k) = u_new(k);
            }
            /*Write to file for plotting*/
            if(tofile && (n%spacing)==0){output(&outfile,u_prev,n,1,nx);}

        }
        cout<<"Backward Euler scheme finished"<<endl;
    }
    if(CN1D){
    //##########################################
    //------Crank Nicolson scheme-------------##
    //##########################################
    
        u_prev = linspace<vec>(-1.0,0.0,nx+1);
        double a1 = -dtdx2;
        double c1 = a1;
        double b1 = 2+2*dtdx2;
        double a2 = dtdx2;
        double c2 = a2;
        double b2 = 2-2*dtdx2;
        u_new = u_prev;
        for(int n=1; n<=n_t; n++){
            make_uprev(u_prev,u_new,a2,c2,b2,nx); 
            u_prev(0)=0;u_prev(nx)=0;
            tridiag(a1,b1,c1, u_new, u_prev,nx);
            //u_prev = u_new;
            for(int k=0;k<=nx;k++){
                u_prev(k) = u_new(k);
            }
            u_prev(0)=0;u_prev(nx)=0;
            /*Write to file for plotting*/
            if(tofile && (n%spacing)==0){output(&outfile,u_prev,n,2,nx);} 
        }
        cout<<"Crank Nicolson scheme finished"<<endl;
    }
    //########################################
    //---------------2D solvers-------------##
    //--------------Forward Euler-----------##
    //########################################
    if(0){
        double C = dtdx2;
        if(C >= 0.25)
        {//Insert for stability criterion!
            C = 0.2;
        }
        dx = 1.0/nx;
        dt = C*dx*dx;
        cout<<"dt = "<<dt<<" C = "<<C<<endl;
        initial_condition(U_p,dx,nx);
        for(int t=0; t<n_t; t++){
            for(int i=1; i<nx; i++){
                for(int j=1; j<nx; j++){

                    U(i,j) = U_p(i,j) + C*(U_p(i+1,j)-2*U_p(i,j)+U_p(i-1,j)) \
                    + C*(U_p(i,j+1)-2*U_p(i,j)+U_p(i,j-1));
                }
            }
            update_boundarys(U,t*dt,dx,nx);
            for(int k=0;k<=nx;k++){
                for(int l=0;l<=nx;l++){
                    U_p(k,l) = U(k,l);
                }
            }
             //Write to file for plotting
            if(tofile && (t%spacing)==0){output2D(&outfile,U_p,t,3,nx);}
        }
        cout<<"Forward Euler done!"<<endl;
    }
    if(LeapFrog){
        double C = dtdx2;
        if (C>(1.0/8.0))
        {
            C = 1.0/10.0;
        }
        dt = C*dx*dx;
        initial_condition(U_p,dx,nx);
        for(int i=1; i<nx; i++){
            for(int j=1; j<nx; j++){
                U(i,j) = U_p(i,j) + C*(U_p(i+1,j)-2*U_p(i,j)+U_p(i-1,j)) \
                + C*(U_p(i,j+1)-2*U_p(i,j)+U_p(i,j-1));
            }
        }
        update_boundarys(U,dt,dx,nx);
        U_pp = U_p;
        U_p = U;
        for(int t=1; t<n_t; t++){
            for(int i=1; i<nx; i++){
                for(int j=1; j<nx; j++){
                    U(i,j) = C*(U_p(i+1,j)-2*U_p(i,j)+U_p(i-1,j)) \
                    + C*(U_p(i,j+1)-2*U_p(i,j)+U_p(i,j-1))+0.5*U_pp(i,j);
                }
            }
            update_boundarys(U,t*dt,dx,nx);
            for(int k=0;k<=nx;k++){
                for(int l=0;l<=nx;l++){
                    U_pp(k,l) = U_p(k,l);
                }
            }
            for(int k=0;k<=nx;k++){
                for(int l=0;l<=nx;l++){
                    U_p(k,l) = U(k,l);
                }
            }
            //Write to file for plotting        
            if(tofile && (t%spacing)==0){output2D(&outfile,U_p,t,4,nx);}
    }
    cout<<"Leap Frog scheme finished"<<endl;
    }
    if(FE2D){
    //################################################
    //----------Euler Chromer scheme----------------##
    //################################################
        double C = dtdx2;
        if(C >= 0.25)
        {//Insert for stability criterion!
            C = 0.2;
        }
        dx = 1.0/nx;
        dt = C*dx*dx;

        cout<<"dt = "<<dt<<" C = "<<C<<endl;
        initial_condition(U_p,dx,nx);
        for(int t=1; t<=n_t; t++){
            for(int i=1; i<nx; i++){
                for(int j=1; j<nx; j++){
                    U(i,j) = U_p(i,j) + C*(U_p(i+1,j)-2*U_p(i,j)+U(i-1,j)) \
                    + C*(U_p(i,j+1)-2*U_p(i,j)+U(i,j-1));
                }
            }
            update_boundarys(U,t*dt,dx,nx);
            for(int k=0;k<=nx;k++){
                for(int l=0;l<=nx;l++){
                    U_p(k,l) = U(k,l);
                }
            }
             //Write to file for plotting
            if(tofile && (t%spacing)==0){output2D(&outfile,U_p,t,3,nx);}
        }
        cout<<"Euler Chromer done!"<<endl;
    }
    return 0;
}        
