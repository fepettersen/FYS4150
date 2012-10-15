/* 
 * File:   main.cpp
 * Author: Fredrik E Pettersen
 * Program description: This program is for project 1 in FYS4150. 
 * The aim is to solve a second order differential equation by linear algebra.
 * Created on 30. august 2012, 08:31
 */


#include "jacobi.h"

int main(int argc, char** argv) {
    
//------ Initialization--------
    
    int n,choose;
    double h,rho_0,rho_n, start, stop;
    rho_0 = 0.0;
//    cout << "Enter the size of the NxN matrix" << endl;
//    cin >> n;
//    cout<<"Enter rho_max"<<endl;
//    cin >> rho_n;
    n = atoi(argv[1]);
    rho_n = atof(argv[2]);
    choose = 1;//choose which potential to use 0 is without repulsive coloumbforce
    double omega = atof(argv[3]); 
    h = (rho_n - rho_0)/(n+1);
    vec rho= linspace<vec>(rho_0+h,rho_n-h,n);
    mat A = make_A(n, rho_0, rho_n, rho, choose,omega);
    mat R(n,n); R.eye();
    double eps = 1e-9;//atof(argv[3]); 
    int k=0;
    int l=0;
    int rotation_counter = 0;
    int max_it = 10*n*n;
    double max_offdiag = maxoffdiag(A, &k, &l, n); 
//------ Solving the equations ----------
    cout<<"Init ok"<<endl;
    
    vec eigval(n); eigval(0) = A(0,0);
    vec nondiag(n);
    mat eigvec(n,n); eigvec.eye();
    for(int i=1;i<n;i++){
        eigval(i) = A(i,i);
        nondiag(i) = A(i-1,i);
    }
    start = clock();
    eigvec = tqli(eigval,nondiag,n,eigvec);
    stop = clock();
    cout<<"tqli function done. "<<timediff(start,stop)<<" ms. "<<endl;
    /*
    start = clock();    //start timing the computation
    while(max_offdiag>eps && rotation_counter<max_it){
        max_offdiag = maxoffdiag(A, &k, &l,n);
        rotate(A,R,k,l,n);
        rotation_counter ++;
    }
    stop = clock();
    */
    //-------- Print some messages to screen-------------
    /*
    cout<<"Yarrr! We be done! We been runnin' "<<rotation_counter<<" rounds in circle. This be takin' ";
    cout<< timediff(start,stop)<<" ms" <<endl;
    cout<<n<<" by "<<n<<" matrix, rho_max= "<<rho_n<<", eps= "<<eps<<", largest element: "<<max_offdiag<<endl;
    */
    cout<<"Done! "<<n<<" by "<<n<<" matrix. rho_max = "<<rho_n<<" eps = "<<eps<<endl;
    /*
    double eigval1 = 0;
    double eigval2 = 0;
    double eigval3 = 0;
    */
    //print_diag(A,3,n, &eigval1, &eigval2, &eigval3);
    vec min  =sort_eigenvalues(eigval, 3, n);
    vec eigenvec0=eigvec.col(min(0));
    vec eigenvec1=eigvec.col(min(1));
    vec eigenvec2=eigvec.col(min(2));
    ofstream myfile;
    myfile.open(make_filename(n,rho_n,omega,0));
    for(int i=0;i<n;i++){
        myfile <<eigenvec0(i)<<"  "<<rho(i)<<" "<<eigenvec1(i)<<" "<<eigenvec2(i)<<endl;
    }
    myfile.close();
   
    return 0;
}