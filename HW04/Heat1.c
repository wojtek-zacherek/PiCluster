#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>

#define DEBUG 0
#define NUM_PER_BLOCK_AXIS 4

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

void generateVector(double *b, int N, int offset){
    int i;
    for(i = 0; i < N; i++){
        if(offset > 0){
            b[i] = 1.0 / (i + 1.0 + 10*(rand()/(double)RAND_MAX));
        }else{
            b[i] = 1;
        }
    }
}

int main(int argc, char *argv[]){
    int my_rank, comm_sz;
    
    MPI_Init(NULL,NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    srand(time(NULL)+ my_rank);

    // Ensure N >= 4

    int N = 50;
    int NN = NUM_PER_BLOCK_AXIS * N;
    int mu = 0.1;

    double *u = mallocVector(NN);
    double *q = mallocVector(NN);

    generateVector(u,NN,1);

    int i, j;
    int indicies[5];
    int internalIndex = -1;
    int blockNum = -1;
    int blockNum2 = -1;
    int blockNum3 = -1;

    printf("Vector u: ");
    printVector(u, NN);

    for(i = 0; i < NN; i++){
        internalIndex = i % 4;
        blockNum = i / 4;
        blockNum2 = (N + blockNum + 1) % N;
        blockNum3 = (N + blockNum - 1) % N;
        indicies[0] = 4*blockNum  + ((4 + internalIndex) % 4);
        indicies[1] = 4*blockNum  + ((4 + internalIndex + 1) % 4);
        indicies[2] = 4*blockNum  + ((4 + internalIndex - 1) % 4);
        indicies[3] = 4*blockNum2 + ((4 + internalIndex) % 4);
        indicies[4] = 4*blockNum3 + ((4 + internalIndex) % 4);
        if (DEBUG){
            printf("%d %d %d %d - %d  %d  %d  %d  %d\n",internalIndex, blockNum, blockNum2, blockNum3, indicies[0], indicies[1], indicies[2], indicies[3], indicies[4]);
        }

        q[i] = 0;
        q[i] += 4*mu;
        for(j = 0; j < 5; j++){
            q[i] += -1*u[indicies[j]];
        }
    }

    printf("Output q: ");
    printVector(q, NN);

    free(q);
    free(u);
    MPI_Finalize();
    return 0;
}