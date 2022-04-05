#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>

#define DEBUG 1
#define DEBUG_KERNEL 0

void printVector(double *vector, int N);

void multiple(double *col, double *vector, double *out, int N, int N2, int my_rank){
    // N = row #.
    // N2 = col #

    int i, j;

    for(i = 0; i < N; i++){  // loop row
        out[i] = 0;
        for(j = 0; j < N2; j++){  // loop col
            out[i] = out[i] + col[i + j*N] * vector[j];
            if(DEBUG_KERNEL && my_rank==1){
                printf("%.2f  x  %.2f  =  %.3f\n", col[j*N + i], vector[j], col[j*N + i]*vector[j]);
            }
        }
        if(DEBUG_KERNEL && my_rank==1){
            printf("%d - %d    %.3f\n", my_rank, i, out[i]);
        }
    }
    if(DEBUG_KERNEL && my_rank == 1){
        printVector(out,N);
    }
}

int check(double *matrix, double *vector, double *vectorResult, int N){
    int i, j, status;
    double *tempOut = malloc(N*sizeof(double));
    for(i = 0; i < N; i++){
        tempOut[i] = 0;
        for(j = 0; j < N; j++){
            tempOut[i] = tempOut[i] + matrix[i*N + j] * vector[j];
        }
    }
    status = 0;
    for(i = 0; i < N; i++){
        if(vectorResult[i] != tempOut[i]){
            printf("Error at row %d: %.2f != %.2f\n", i, vectorResult[i], tempOut[i]);
            status = 1;
        }
    }
    free(tempOut);
    return status;
}

void printVector(double *vector, int N){
    int i;
    for(i = 0; i < N; i++){
        printf("  %.3f",vector[i]);
    }
    printf("\n");
    printf("-- end printout --\n\n");
}

void printMatrix(double *matrix, int N){
    int i, j;
    for(i = 0; i < N; i++){
        for(j = 0; j < N; j++){
            printf("  %.3f",matrix[N*i + j]);
        }
        printf("\n");
    }
    printf("-- end printout --\n\n");
}

void printMatrixColWise(double *matrix, int N){
    int i, j;
    for(i = 0; i < N; i++){
        for(j = 0; j < N; j++){
            printf("  %.3f",matrix[i + N*j]);
        }
        printf("\n");
    }
    printf("-- end printout --\n\n");
}

void transposeSquareMatrix(double *matrix, double *out, int N){
    int i, j;
    for(i = 0; i < N; i++){
        for(j = 0; j < N; j++){
            out[i*N + j] = matrix[i + j*N];
        }
    }
}

int main(int argc, char *argv[]){
    int my_rank, comm_sz;
    int i;

    MPI_Init(NULL,NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    srand(time(NULL));

    double *vector, *matrix, *col, *matrixTransposed;
    double *singleResult, *vectorResult;

    // Assumptions:
    //      N % comm_sz = 0
    //      Block partitioning (1 matrix col per proc)
    int N;
    int alpha = 2;
    N = alpha * comm_sz;

    vector = malloc(N * sizeof(double));
    col = malloc(alpha * N * sizeof(double));
    singleResult = malloc(alpha * N * sizeof(double));

    // Send vector
    // Scatter matrix
    if(my_rank == 0){
        if(DEBUG){
            printf("Number nodes = %d\n",comm_sz);
            printf("N = %d\n",N);
        }

        matrix = malloc(N*N * sizeof(double));
        matrixTransposed = malloc(N*N * sizeof(double));
        vectorResult = malloc(N * sizeof(double));
        for(i = 0; i < N*N; i++){
            matrix[i] = rand() % 20;
        }
        for(i = 0; i < N; i++){
            vector[i] = rand() % 5;
        }
        transposeSquareMatrix(matrix,matrixTransposed,N);
        
        if(DEBUG){
            printf("--Vector--\n");
            printVector(vector,N);
            printf("--Matrix--\n");
            printMatrix(matrix,N);
            printf("--Transposed--\n");
            printMatrix(matrixTransposed,N);
        }
    }else{

    }
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Bcast(vector,N,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    if(my_rank == 0 && DEBUG){
        printf("Vector transfer done\n");
    }

    MPI_Scatter(matrixTransposed, N*alpha, MPI_DOUBLE, col, N*alpha, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    if(my_rank == 0 && DEBUG){
        printf("Matrix transfer done\n");
    }

    multiple(col, vector + my_rank*alpha, singleResult, N, alpha, my_rank);
    MPI_Barrier(MPI_COMM_WORLD);
    if(my_rank == 0 && DEBUG){
        printf("Multiplication done.\n");
    }

    MPI_Reduce(singleResult, vectorResult, N, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    if(my_rank == 0 && DEBUG){
        printf("Reduction summation done.\n");
    }

    if(my_rank == 0){
        printf("--Multiplcation Result--\n");
        printVector(vectorResult, N);
        int status = check(matrix, vector, vectorResult, N);
        if(status){
            printf("Errors occured. Good luck finding them :) .\n");
        }else{
            printf("No errors, computation is correct.\n");
        }
    }

    

    if(my_rank == 0){
        free(matrix);
        free(matrixTransposed);
        free(vectorResult);
    }
    free(vector);
    free(col);
    free(singleResult);

    MPI_Finalize();
    return 0;
}
