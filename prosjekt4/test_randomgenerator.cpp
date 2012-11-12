#include "ising.h"
int main(int argc, char** argv) {
	long idum = -1;
	for(int k =0;k<10000000;k++){
        cout<<int (ran3(&idum))<<", ";
    }
    cout<<endl;
    return 0;
}