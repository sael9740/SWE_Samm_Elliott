#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define DIMENSION 201
#define XYZ_MAX 1
#define TOLERANCE 1

typedef float Domain_t[DIMENSION][DIMENSION][DIMENSION];

float f(int i, int j, int k); // solution function
void init_jacobi(Domain_t A); // function to initialize boundary
void init_sol(Domain_t A); // function to initialize solution
float max_diff(Domain_t A, Domain_t B); // calculates max difference of values in two domains
void do_jacobi(Domain_t A, Domain_t B);



int main(int argc, char** argv)
{
    int iter = 0;
    float t_start, t_total;
    //printf("made it here!");
    float err = TOLERANCE +1;

    Domain_t jacobi_A, jacobi_B, real_sol;
    
    // initialize boundaries and real solution
    init_jacobi(jacobi_A);
    init_jacobi(jacobi_B);
    init_sol(real_sol);
    #ifdef _OPENACC
        #pragma acc data copyin(jacobi_A,jacobi_B,real_sol)
    #endif
    {
    t_start = omp_get_wtime();
    
    while(err > TOLERANCE) {
        iter += 2;
        
        do_jacobi(jacobi_A,jacobi_B);
        err = max_diff(jacobi_A,real_sol);
        
        if(iter%100 ==0)
            printf("Error: %f\n",err);
    }
    }
    t_total = omp_get_wtime() - t_start;
    printf("Converged to a maximum error of %f in %f seconds and %d iterations with %d total threads\n",TOLERANCE,t_total,iter,omp_get_max_threads());
    
    
    return(0);
}


// solution function
float f(int i, int j, int k)
{
    float x,y,z;
    
    x = ((float) i)*XYZ_MAX/(DIMENSION-1);
    y = ((float) j)*XYZ_MAX/(DIMENSION-1);
    z = ((float) k)*XYZ_MAX/(DIMENSION-1);
    
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
float max_diff(Domain_t A, Domain_t B)
{
    int i,j,k;
    float err, max_val = 0;
    
    #ifdef _OPENACC
        #pragma acc parallel loop private(j,k,err) reduction(max:max_val) independent collapse(3) present_or_copyin(A[:DIMENSION][:DIMENSION][:DIMENSION],B[:DIMENSION][:DIMENSION][:DIMENSION])
    #else
        #pragma omp parallel for private(j,k,err) reduction(max:max_val) collapse(3)
    #endif
    for(i=0;i<DIMENSION;i++) {
        for(j=0;j<DIMENSION;j++) {
            for(k=0;k<DIMENSION;k++) {
                err = fabs(A[i][j][k]-B[i][j][k]);
                if(err > max_val) {
                    max_val = err;
                }
            }
        }
    }
    
    printf("%f\n",max_val);
    
    return max_val;
}

void do_jacobi(Domain_t A, Domain_t B)
{
    int i,j,k;
    #ifdef _OPENACC
        #pragma acc parallel private(i,j,k) present_or_copyin(A[:DIMENSION][:DIMENSION][:DIMENSION],B[:DIMENSION][:DIMENSION][:DIMENSION])
    {
    #endif
    
    #ifdef _OPENACC
        #pragma acc loop  independent collapse(3)
    #else
        #pragma omp parallel for collapse(3)
    #endif
    for(i=1;i<DIMENSION-1;i++) {
        for(j=1;j<DIMENSION-1;j++) {
            for(k=1;k<DIMENSION-1;k++) {
                B[i][j][k]=(A[i-1][j][k]+A[i+1][j][k]+A[i][j-1][k]+A[i][j+1][k]+A[i][j][k-1]+A[i][j][k+1])/6;
            }
        }
    }
    #ifdef _OPENACC
        #pragma acc wait
    #endif
        
    #ifdef _OPENACC
        #pragma acc loop  independent collapse(3)
    #else
        #pragma omp parallel for collapse(3)
    #endif
    for(i=1;i<DIMENSION-1;i++) {
        for(j=1;j<DIMENSION-1;j++) {
            for(k=1;k<DIMENSION-1;k++) {
                A[i][j][k]=(B[i-1][j][k]+B[i+1][j][k]+B[i][j-1][k]+B[i][j+1][k]+B[i][j][k-1]+B[i][j][k+1])/6;
            }
        }
    }
    #ifdef _OPENACC
        #pragma acc wait
    }
    #endif
}








