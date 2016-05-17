#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define DIMENSION 101
#define XYZ_MAX 1
#define TOLERANCE .1

typedef double Domain_t[DIMENSION][DIMENSION][DIMENSION];

double f(int i, int j, int k); // solution function
void init_jacobi(Domain_t A); // function to initialize boundary
void init_sol(Domain_t A); // function to initialize solution
double max_diff(Domain_t A, Domain_t B); // calculates max difference of values in two domains
void do_jacobi(Domain_t A, Domain_t B);



int main(int argc, char** argv)
{
    double t_start, t_total;
    //printf("made it here!");
    double err = TOLERANCE +1;
    
    Domain_t* jacobi_A = (Domain_t *) malloc(sizeof(Domain_t));
    Domain_t* jacobi_B = (Domain_t *) malloc(sizeof(Domain_t));
    Domain_t* real_sol = (Domain_t *) malloc(sizeof(Domain_t));
    Domain_t* dummy;
    
    // approximation and real solution arrays
    //double jacobi_A[DIMENSION][DIMENSION][DIMENSION];
    //double jacobi_B[DIMENSION][DIMENSION][DIMENSION];
    //double real_sol[DIMENSION][DIMENSION][DIMENSION];
    
    // initialize boundaries and real solution
    init_jacobi(*jacobi_A);
    init_jacobi(*jacobi_B);
    init_sol(*real_sol);
    
    t_start = omp_get_wtime();
    
    while(err > TOLERANCE) {
        do_jacobi(*jacobi_A,*jacobi_B);
        dummy = jacobi_A;
        jacobi_A = jacobi_B;
        jacobi_B = dummy;
        err = max_diff(*jacobi_A,*real_sol);
        printf("Error: %f\n",err);
    }
    
    t_total = omp_get_wtime() - t_start;
    printf("Converged in %f seconds",t_total);
    
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

void init_jacobi(Domain_t A)
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
void init_sol(Domain_t A)
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
double max_diff(Domain_t A, Domain_t B)
{
    int i,j,k;
    double max_val = 0;

    #pragma omp parallel for private(j,k) reduction(max:max_val)
    for(i=0;i<DIMENSION;i++) {
        for(j=0;j<DIMENSION;j++) {
            for(k=0;k<DIMENSION;k++) {
                if(fabs(A[i][j][k]-B[i][j][k]) > max_val) {
                    max_val = fabs(A[i][j][k]-B[i][j][k]);
                }
            }
        }
    }
    
    return max;
}

void do_jacobi(Domain_t A, Domain_t B)
{
    int i,j,k;
    
#pragma omp parallel for private(i,j,k)
    for(i=1;i<DIMENSION-1;i++) {
        for(j=1;j<DIMENSION-1;j++) {
            for(k=1;k<DIMENSION-1;k++) {
                B[i][j][k]=(A[i-1][j][k]+A[i+1][j][k]+A[i][j-1][k]+A[i][j+1][k]+A[i][j][k-1]+A[i][j][k+1])/6;
            }
        }
    }
}