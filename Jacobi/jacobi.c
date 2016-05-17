#include <stdio.h>

#define DIMENSION 101
#define XYZ_MAX 1

double f(double i, double j, double k); // solution function
void init_jacobi(double A[DIMENSION][DIMENSION][DIMENSION]); // function to initialize boundary
//double init_sol(double A[DIMENSION][DIMENSION][DIMENSION],i_max,j_max,k_max); // function to initialize solution



int main(int argc, char** argv)
{
    printf("made it here!");
    // approximation and real solution arrays
    double jacobi_A[DIMENSION][DIMENSION][DIMENSION];
    double jacobi_B[DIMENSION][DIMENSION][DIMENSION];
    double real_sol[DIMENSION][DIMENSION][DIMENSION];
    
    // initialize boundaries and real solution
    init_jacobi(jacobi_A);
    
    
    return(0);
}


// solution function
double f(double i, double j, double k)
{
    int x = i*XYZ_MAX/(DIMENSION-1);
    int y = j*XYZ_MAX/(DIMENSION-1);
    int z = k*XYZ_MAX/(DIMENSION-1);
    return(x+y+z);
}

void init_jacobi(double A[DIMENSION][DIMENSION][DIMENSION])
{
    int i,j,k;
    printf("made it here!");
    
    // init inside to 0
    for(i=1;i<DIMENSION;i++) {
        for(j=1;j<DIMENSION;j++) {
            for(k=1;k<DIMENSION;k++) {
                A[i][j][k]=0;
            }
        }
    }
    
    printf("made it here!");
    
    //init boundaries to match f(x,y,z)
    for(i=0;i<DIMENSION;i++) {
        for(j=0;j<DIMENSION;j++) {
            A[i][j][0]=f(i,j,0);
            A[i][j][DIMENSION-1]=f(i,j,DIMENSION-1);
        }
    }
    for(i=0;i<DIMENSION;i++) {
        for(k=0;k<DIMENSION;k++) {
            A[i][0][k]=f(i,0,k);
            A[i][DIMENSION-1][k]=f(i,DIMENSION-1,k);
        }
    }
    for(j=0;j<DIMENSION;j++) {
        for(k=0;k<DIMENSION;k++) {
            A[0][j][k]=f(0,j,k);
            A[DIMENSION-1][j][k]=f(DIMENSION-1,j,k);
        }
    }
}


