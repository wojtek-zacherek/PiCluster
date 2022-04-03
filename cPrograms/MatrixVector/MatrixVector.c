#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>

#define DEBUG 1

void multiple(double *row, double *vector, double *out, int N, int N2){
    double sum = 0;
    int i, j;
    for(j = 0; j < N2; j++){        // Row select
        sum = 0;
        for(i = 0; i < N; i++){     // Col select
            sum+= row[j*N + i] * vector[i];
        }
        *(out + j) = sum;
    }
}

void printVector(double *vector, int N){
    int i;
    for(i = 0; i < N; i++){
        printf("%.3f  ",vector[i]);
    }
    printf("\n");
    printf("-- end printout --\n");
}

void printMatrix(double *matrix, int N){
    int i, j;
    for(i = 0; i < N; i++){
        for(j = 0; j < N; j++){
            printf("%.3f  ",matrix[N*i + j]);
        }
        printf("\n");
    }
    printf("-- end printout --\n");
}

int main(int argc, char *argv[]){
    int my_rank, comm_sz;
    int i;
    
    MPI_Init(NULL,NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    srand(time(NULL));

    double *vector, *matrix, *row;
    double *singleResult, *vectorResult;
    
    // Assumptions: 
    //      N % comm_sz = 0
    //      Block partitioning (1 matrix row per proc)
    int N;
    int alpha = 2;
    N = alpha * comm_sz;

    vector = malloc(N * sizeof(double));
    row = malloc(alpha * N * sizeof(double));
    singleResult = malloc(alpha * sizeof(double));

    // Send vector
    // Scatter matrix
    if(my_rank == 0){
        if(DEBUG){
            printf("Number nodes = %d\n",comm_sz);
            printf("N = %d\n",N);
        }
        matrix = malloc(N*N * sizeof(double));
        vectorResult = malloc(N * sizeof(double));  
        for(i = 0; i < N*N; i++){
            matrix[i] = rand() % 20;
        }
        for(i = 0; i < N; i++){
            vector[i] = rand() % 5;
        }
    }else{
    
    }

    printf("Proc %d eached first barrier\n",my_rank);
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Bcast(vector,N,MPI_DOUBLE,0,MPI_COMM_WORLD);
    if(my_rank == 0 && DEBUG){
        printf("Vector transfer done\n");
    }
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Scatter(matrix, N*alpha, MPI_DOUBLE,
                row, N*alpha, MPI_DOUBLE,
                0, MPI_COMM_WORLD);
    if(my_rank == 0 && DEBUG){
        printf("Matrix transfer done\n");
    }
    MPI_Barrier(MPI_COMM_WORLD);

    multiple(row, vector, singleResult, N, alpha);
    if(my_rank == 0 && DEBUG){
        printf("Multiplication done\n");
    }   
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Gather(singleResult, alpha, MPI_DOUBLE,
                vectorResult, alpha, MPI_DOUBLE,
                0,MPI_COMM_WORLD);
    if(my_rank == 0 && DEBUG){
        printf("Gather done\n");
    }
    MPI_Barrier(MPI_COMM_WORLD);

    if(my_rank == 0){
        printVector(vector,N);
        printMatrix(matrix,N);
        printVector(vectorResult, N);
    }

    if(my_rank == 0){
        free(matrix);
        free(vectorResult);
    }
    free(vector);
    free(row);
    free(singleResult);

    MPI_Finalize();
    return 0;
}
