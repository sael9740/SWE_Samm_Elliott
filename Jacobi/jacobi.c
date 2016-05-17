#include <stdio.h>
#include <math.h>

#define DIMENSION 101
#define XYZ_MAX 1

double f(int i, int j, int k); // solution function
void init_jacobi(double A[DIMENSION][DIMENSION][DIMENSION]); // function to initialize boundary
void init_sol(double A[DIMENSION][DIMENSION][DIMENSION]); // function to initialize solution
double max_diff(double A[DIMENSION][DIMENSION][DIMENSION],
                double B[DIMENSION][DIMENSION][DIMENSION]); // calculates max difference of values in two domains



int main(int argc, char** argv)
{
    //printf("made it here!");
    
    // approximation and real solution arrays
    double jacobi_A[DIMENSION][DIMENSION][DIMENSION];
    double jacobi_B[DIMENSION][DIMENSION][DIMENSION];
    double real_sol[DIMENSION][DIMENSION][DIMENSION];
    
    // initialize boundaries and real solution
    init_jacobi(jacobi_A);
    init_sol(real_sol);
    
    max_diff(jacobi_A,real_sol);
    
    return(0);
}


// solution function
double f(int i, int j, int k)
{
    double x,y,z;
    
    x = ((double) i)*XYZ_MAX/(DIMENSION-1);
    y = ((double) j)*XYZ_MAX/(DIMENSION-1);
    z = ((double) k)*XYZ_MAX/(DIMENSION-1);
    
    return(x+y+z);
}

void init_jacobi(double A[DIMENSION][DIMENSION][DIMENSION])
{
    int i,j,k;
    //printf("made it here!");
    
    // init inside to 0
    for(i=1;i<DIMENSION-1;i++) {
        for(j=1;j<DIMENSION-1;j++) {
            for(k=1;k<DIMENSION-1;k++) {
                A[i][j][k]=0;
            }
        }
    }
    
    // init boundaries to match f(x,y,z)
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

// initializes real solution
void init_sol(double A[DIMENSION][DIMENSION][DIMENSION])
{
    int i,j,k;
    
    // init inside to 0
    for(i=0;i<DIMENSION;i++) {
        for(j=0;j<DIMENSION;j++) {
            for(k=0;k<DIMENSION;k++) {
                A[i][j][k]=f(i,j,k);
            }
        }
    }
}

// calculates max difference of the values in two domains
double max_diff(double A[DIMENSION][DIMENSION][DIMENSION],
                double B[DIMENSION][DIMENSION][DIMENSION])
{
    int i,j,k;
    int imax, jmax, kmax = 1000;
    double max = 0;
    
    // init inside to 0
    for(i=0;i<DIMENSION;i++) {
        for(j=0;j<DIMENSION;j++) {
            for(k=0;k<DIMENSION;k++) {
                if(labs(A[i][j][k]-B[i][j][k]) > max) {
                    max = labs(A[i][j][k]-B[i][j][k]);
                    imax=i;jmax=j;kmax=k;
                    printf("maximum difference: %f at (%d,%d,%d)\n f(i,j,k)=%f\nA[i][j][k]=%f B[i][j][k]=%f\n",max,imax,jmax,kmax,f(imax,jmax,kmax),A[i][j][k],B[i][j][k]);
                }
                //if(i==j&&j==k)
                //    printf("at (%d,%d,%d), the value of A_max is %f and B_max is %f and f(i,j,k)=%f\n",i,j,k,A[i][j][k],B[i][j][k],f(i,j,k));
            }
        }
    }
    
    printf("maximum difference: %f at (%d,%d,%d)\nf(i,j,k)=%f\n",max,imax,jmax,kmax,f(imax,jmax,kmax));
    printf("value of A_max is %f and B_max is %f\n",A[imax][jmax][kmax],B[imax][jmax][kmax]);
    return max;
}

