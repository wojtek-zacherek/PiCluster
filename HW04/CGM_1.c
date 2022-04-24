#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>
#include <math.h>

#define DEBUG 1

void printVector(double *v1, int N){
    int i;
    for(i = 0; i < N; i++){
        printf("%f ",v1[i]);
    }
    printf("\n");
}

double** mallocMatrix(int numRows, int numCols){
    double **matrixPtr = (double**)malloc(numRows * sizeof(double*));
    int i;
    for (i = 0; i < numRows; i++){
        matrixPtr[i] = (double*)malloc(numCols * sizeof(double));
    }
    return matrixPtr;
}

void freeMatrix(double **A, int numRows, int numCols){
    int i;
    for (i = 0; i < numRows; i++){
        free(A[i]);
    }
    free(A);
}

double* mallocVector(int N){
    return (double *)malloc(N * sizeof(double));
}

void multVector(double **A,double *x,double *Ax,int N){
    int i,j;
    for(i = 0; i < N; i++){
        Ax[i] = 0;
        for(j = 0; j < N; j++){
            Ax[i] += A[i][j] * x[j];
        }
    }
}

void addVector(double *a, double *b, double *sum, int N){
    int i;
    for(i = 0; i < N; i++){
        sum[i] = a[i] + b[i];
    }
}

double scalarInnerProduct(double *v1, double *v2, int N){
    int i;
    double result = 0.0;
    for(i = 0; i < N; i++){
        result += v1[i]*v2[i];
    }
    return result;
}

void multVectScalar(double scalar, double *v1, double *product, int N){
    int i;
    for(i = 0; i < N; i++){
        product[i] = scalar*v1[i];
    }
}

double vectorNorm(double *v1, int N){
    int i;
    double result = 0;
    for(i = 0; i < N; i++){
        result += v1[i];
    }
    result = sqrt(result);
    return result;
}

void axpy(double **A, double *x, double *b, double* out, int N){
    double *temp = mallocVector(N);
    multVector(A,x,temp,N);
    addVector(temp,b,out,N);
    free(temp);
}

void saxpy(double scalar, double *x, double *b, double *out, int N){
    double *temp = mallocVector(N);
    multVectScalar(scalar, x, temp, N);
    addVector(temp, b, out, N);
    free(temp);
}

void copyArray(double *src, double *dest, int N){
    int i;
    for(i = 0; i < N; i++){
        dest[i] = src[i];
    }
}

void generateA(double **A, int N){
    int i,j;
    for(i = 0; i < N; i++){
        for(j = 0; j < N; j++){
            // A[i][j] = 1.0/(i + j + 1.0);
            if(i==j)
                A[i][j] = 1.0;
            else
                A[i][j] = 0.0;
        }
    }
}

void generateVector(double *b, int N, int offset){
    int i;
    for(i = 0; i < N; i++){
        if(offset > 0){
            b[i] = 1.0 / (i + 1.0 + (rand()/(double)RAND_MAX));
            b[i] = 0.10;
        }else{
            // b[i] = 1.0 / (i + 1.0);
            b[i] = 1;
        }
    }
}


int main(int argc, char *argv[]){
    int my_rank, comm_sz;

    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    srand(time(NULL) + my_rank);

    int N = 2;
    int k = 0;
    int kMax = 4*N;

    double *b = mallocVector(N);
    double *x0 = mallocVector(N);
    double *r0 = mallocVector(N);
    double *p0 = mallocVector(N);
    double **A = mallocMatrix(N,N);

    generateA(A,N);
    generateVector(b, N, 0);
    generateVector(x0, N, 1);
    axpy(A, x0, b, r0, N);
    copyArray(r0, p0, N);

    double *p_current = mallocVector(N);
    double *v_current = mallocVector(N);
    double *r_current = mallocVector(N);
    double *x_current = mallocVector(N);
    double *x_next = mallocVector(N);
    double *r_next = mallocVector(N);
    double *p_next = mallocVector(N);
    double *finalResult = mallocVector(N);
    double alpha_current;
    double tolerance = 0.0001;
    double residual;
    double beta_current;

    copyArray(p0, p_current, N);
    copyArray(r0, r_current, N);
    copyArray(x0, x_current, N);

    while(k < kMax){
        if(k >= kMax){break;}
        printf("-------k = %d--------\n",k);
        multVector(A, p_current, v_current, N);     // Step 8
        printf("v%d: ",k);
        printVector(v_current,N);
        
        alpha_current = scalarInnerProduct(r_current, r_current, N) / scalarInnerProduct(p_current, v_current, N);        // Step 9
        printf("alpha_%d: %f\n",k, alpha_current);

        saxpy(alpha_current, p_current, x_current, x_next, N);      // Step 10
        printf("x%d: ",k+1);
        printVector(x_next,N);

        alpha_current *= -1;        // Step 11
        saxpy(-1*alpha_current, v_current, r_current, r_next, N);
        alpha_current *= -1;
        printf("r%d: ",k+1);
        printVector(r_next,N);
        
        residual = vectorNorm(r_next, N);        // Step 12
        if(residual < tolerance){break;}
        printf("res%d: %f\n",k, residual);

        beta_current = scalarInnerProduct(r_next, r_next, N) / scalarInnerProduct(r_current, r_current, N);     // Step 13
        printf("beta%d: %f\n",k, beta_current);

        saxpy(beta_current, p_current, r_next, p_next, N);      // Step 14
        printf("p%d: ",k+1);
        printVector(p_next,N);

        k = k + 1;      // Step 15

        copyArray(x_next, x_current, N);
        copyArray(r_next, r_current, N);
        copyArray(p_next, p_current, N);
    }
    printf("------- end --------\n");

    copyArray(x_next, finalResult, N);

    printf("Target: \n");
    printVector(b, N);

    printf("Calculated X: \n");
    printVector(finalResult, N);

    printf("Calculated value: \n");
    double *debugRes = mallocVector(N);
    multVector(A, finalResult, debugRes, N);
    printVector(debugRes, N);
    free(debugRes);

    free(b);
    free(x0);
    free(r0);
    free(p0);
    free(p_current);
    free(v_current);
    free(r_current);
    free(x_current);
    free(x_next);
    free(r_next);
    free(p_next);
    free(finalResult);
    freeMatrix(A,N,N);

    MPI_Finalize();
    return 0;
}