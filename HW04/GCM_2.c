#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>
#include <math.h>

#define DEBUG 1

void printVector(double *v1, int N){
    int i;
    for(i = 0; i < N; i++){
        printf("%.3f ",v1[i]);
    }
    printf("\n");
}
void printMatrix(double **A, int N){
    int i,j;
    for(i = 0; i < N; i++){
        for(j = 0; j < N; j++){
            printf("%.3f ",A[i][j]);
        }
        printf("\n");
    }
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
        result += v1[i]*v1[i];
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
                A[i][j] = 1.0/(N*(j+1));
        }
    }
}

void generateVector(double *b, int N, int offset){
    int i;
    for(i = 0; i < N; i++){
        if(offset > 0){
            b[i] = 1.0 / (i + 1.0 + 10*(rand()/(double)RAND_MAX));
            b[i] = 0.0;
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

    int N = 30;
    int k = 1;
    int kMax = N;
    double residuals[kMax];

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
    double *v = mallocVector(N);
    double *r_current = mallocVector(N);
    double *x_current = mallocVector(N);
    double *finalResult = mallocVector(N);
    double alpha_current;
    double tolerance = 0.00000001;
    double residual;
    double beta_current;
    double alpha_next;
    double *r_previous = mallocVector(N);
    double *p_previous = mallocVector(N);
    double *x_previous = mallocVector(N);
    double alpha_previous;
    double *x_diff = mallocVector(N);

    copyArray(p0, p_previous, N);
    copyArray(r0, r_previous, N);
    copyArray(x0, x_previous, N);

    printf("b%d: ",0);
    printVector(b,N);
    printf("A: \n");
    printMatrix(A,N);

    while(k < kMax){
        if(k >= kMax){break;}
        printf("-------k = %d--------\n",k);
        
        multVector(A, p_previous, v, N);
        alpha_previous = pow(vectorNorm(r_previous,N),2) / scalarInnerProduct(p_previous,v,N);

        saxpy(alpha_previous,p_previous,x_previous,x_current,N);
        saxpy(-1*alpha_previous,v,r_previous,r_current,N);

        multVectScalar(-1, x_current, x_current, N);
        addVector(x_previous, x_current,x_diff,N);
        multVectScalar(-1, x_current, x_current, N);
        residual = vectorNorm(x_diff,N);
        residuals[k-1] = residual;
        if(residual < tolerance){
            printf("Tolerance reached\n");
            break;
        }
        

        beta_current = pow(vectorNorm(r_current,N),2) / pow(vectorNorm(r_previous,N),2);
        saxpy(beta_current,p_previous,r_current,p_current,N);

        copyArray(r_current,r_previous,N);
        copyArray(p_current,p_previous,N);
        copyArray(x_current,x_previous,N);

        k++;
    }
    printf("------- end --------\n");

    printf("Residuals: ");
    printVector(residuals,kMax);

    copyArray(x_current, finalResult, N);

    printf("Target: ");
    printVector(b, N);

    printf("Calculated X: ");
    printVector(finalResult, N);

    printf("Calculated value: ");
    double *debugRes = mallocVector(N);
    multVector(A, finalResult, debugRes, N);
    printVector(debugRes, N);
    free(debugRes);

    free(b);
    free(x0);
    free(r0);
    free(p0);
    free(p_current);
    free(v);
    free(r_current);
    free(x_current);
    free(finalResult);
    free(r_previous);
    free(p_previous);
    free(x_previous);
    free(x_diff);
    freeMatrix(A,N,N);

    MPI_Finalize();
    return 0;
}