#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>
#include <math.h>

#define DEBUG 1
#define DEBUG_KERNEL 0

void printVector(double *vector, int N);

void multiple(double *matrix, double *vector, double *out, int alpha, int block_row, int my_rank){
    // N = row #.
    // N2 = col #

    int i, j;
    for(i = 0; i < alpha; i++){
        out[i + block_row*alpha] = 0;
        for(j = 0; j < alpha; j++){
                // printf("----%d Testing Matrix val: %.2f x %.2f\n",my_rank, matrix[i*alpha + j], vector[j]);
            out[i + block_row*alpha] = out[i + block_row*alpha] + matrix[i*alpha + j] * vector[j];
        }
        // printf("%d: %.2f\n",my_rank,out[i + block_row*alpha]);
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
        printf("  %.2f",vector[i]);
    }
    printf("\n");
    printf("-- end printout --\n\n");
}

void printMatrix(double *matrix, int N){
    int i, j;
    for(i = 0; i < N; i++){
        for(j = 0; j < N; j++){
            printf("  %.2f",matrix[N*i + j]);
        }
        printf("\n");
    }
    printf("-- end printout --\n\n");
}

void printMatrixColWise(double *matrix, int N){
    int i, j;
    for(i = 0; i < N; i++){
        for(j = 0; j < N; j++){
            printf("  %.2f",matrix[i + N*j]);
        }
        printf("\n");
    }
    printf("-- end printout --\n\n");
}

void printSubMatricies(double *matrix, int blocks, int elementsPerBlock){
    int i, j;
    for(i = 0; i < blocks; i++){
        for(j = 0; j < elementsPerBlock; j++){
            printf("  %.2f",matrix[i*elementsPerBlock + j]);
        }
        printf("\n");
    }
    printf("-- end printout --\n\n");
}

void organizeSubMatricies(double *matrix, double *out, int N, int comm_sz, int alpha, int comm_sqr){ 
    // 0 1 2 9 10 11 18 19 20
    int i, j, k;
    for(i = 0; i < comm_sz; i++){
        for(j = 0; j < alpha*alpha; j++){
            int block = comm_sz;
            int internal_col = (j % alpha);
            int internal_row = (j / alpha);
            // offset
            // block
            int block_col = i % comm_sqr;   // block col
            int block_row = i / comm_sqr;   // block row

            k = (internal_row * N) + block_row * N * alpha + 
                (internal_col) + block_col * alpha;
            printf("%d %d     %d %d     %d\n",internal_col,internal_row,block_col,block_row,k);

            out[i*alpha*alpha + j] = matrix[k];
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

    double *vector, *matrix, *matrixReorganized, *subMatrix;
    double *singleResult, *vectorResult;

    // Assumptions:
    //      N % comm_sz = 0.
    //      comm_sz is a perfect square.
    //      Block partitioning 
    //      alpha is block size. blocks are squares
    int N;
    int alpha = 7;
    int numBlocks = comm_sz;
    int numBlocks1Dim = sqrt(numBlocks);
    int elementsPerBlock = alpha*alpha;
    N = alpha * numBlocks1Dim;

    vector = malloc(N * sizeof(double));
    subMatrix = malloc(elementsPerBlock * sizeof(double));
    singleResult = malloc(alpha * N * sizeof(double));
    for(i = 0; i < alpha*N; i++){
        singleResult[i] = 0;
    }

    // Send vector
    // Scatter matrix
    if(my_rank == 0){
        if(DEBUG){
            printf("Number nodes = %d\n",comm_sz);
            printf("N = %d\n",N);
        }

        matrix = malloc(N*N * sizeof(double));
        matrixReorganized = malloc(N*N * sizeof(double));
        vectorResult = malloc(N * sizeof(double));
        for(i = 0; i < N*N; i++){
            // matrix[i] = rand() % 20;
            matrix[i] = i;
        }
        for(i = 0; i < N; i++){
            // vector[i] = rand() % 5;
            vector[i] = i;
        }
        // transposeSquareMatrix(matrix,matrixReorganized,N);
        organizeSubMatricies(matrix,matrixReorganized,N,comm_sz,alpha, numBlocks1Dim);
        
        if(DEBUG){
            printf("--Vector--\n");
            printVector(vector,N);
            printf("--Matrix--\n");
            printMatrix(matrix,N);
            printf("--Transposed--\n");
            printSubMatricies(matrixReorganized,numBlocks,elementsPerBlock);
        }
    }else{

    }
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Bcast(vector,N,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    if(my_rank == 0 && DEBUG){
        printf("Vector transfer done\n");
    }

    MPI_Scatter(matrixReorganized, elementsPerBlock, MPI_DOUBLE, subMatrix, elementsPerBlock, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    printf("   %d blockRow: %d\n",my_rank, my_rank/numBlocks1Dim);
    printf("   %d vectvalue: %d\n",my_rank, (alpha * (my_rank % numBlocks1Dim)));


    if(my_rank == 0 && DEBUG){
        printf("Matrix transfer done\n");
    }

    multiple(subMatrix, vector + (alpha * (my_rank % numBlocks1Dim)), singleResult, alpha, my_rank / numBlocks1Dim, my_rank);
    
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
        free(matrixReorganized);
        free(vectorResult);
    }
    free(vector);
    free(subMatrix);
    free(singleResult);

    MPI_Finalize();
    return 0;
}
