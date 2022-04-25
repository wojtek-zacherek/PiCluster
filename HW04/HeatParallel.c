#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>
#include <math.h>

#define DEBUG 1
#define NUM_PER_BLOCK_AXIS 4

void printVector(double *v1, int N){
    int i;
    for(i = 0; i < N; i++){
        printf("%.3f ",v1[i]);
    }
    printf("\n");
}
void printVectorInt(int *v1, int N){
    int i;
    for(i = 0; i < N; i++){
        printf("%d ",v1[i]);
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
    printf("Hello\n");
    int my_rank, comm_sz;
    
    MPI_Init(NULL,NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    srand(time(NULL)+ my_rank);
    // Ensure N >= 4

    int N;
    int NN;
    double mu = 0.1;

    double *u = mallocVector(NN);
    double *q = mallocVector(NN);
    double *qFinal = mallocVector(NN);

    int i, j;
    int indicies[5];
    int internalIndex = -1;
    int blockNum = -1;
    int blockNum2 = -1;
    int blockNum3 = -1;

    // // {numBlocksN, startingBlock#(inclusive), endingBlock#(exclusive)}
    int _a = 3;
    int runData[_a * comm_sz];
    int *myData;
    myData = malloc(_a * sizeof(int));
    if(my_rank == 0){
        printf("Hello2\n");
        N = 4;
        NN = NUM_PER_BLOCK_AXIS * N;
        generateVector(u, NN, 1);
        
        int blocksPerProc = N / comm_sz;
        int blocksExtra = N % comm_sz;
        int currentCounter = 0;
        for(i = 0; i < comm_sz; i++){
            runData[_a*i + 0] = N;
            runData[_a*i + 1] = currentCounter;
            if(blocksExtra > 0){
                currentCounter += blocksPerProc + 1;
                blocksExtra--;
            }else{
                currentCounter += blocksPerProc;
            }
            runData[_a*i + 2] = currentCounter;
        }
        
        if(DEBUG){
            printf("Number nodes = %d\n",comm_sz);
        }
        if(DEBUG){
            printf("N = %d\n",N);
            printf("NN = %d\n",NN);
            printf("mu = %f\n",mu);
            printf("Labor division: ");
            printVectorInt(runData, _a * comm_sz);
            printf("Vector u: ");
            printVector(u, NN);
        }
        printf("done 0\n");
    }
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Scatter(runData, _a, MPI_INT, myData, _a, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    if(my_rank == 0 && DEBUG){
        printf("Division of labor done\n");
    }

    MPI_Bcast(u, NN, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    if(my_rank == 0 && DEBUG){
        printf("Vector transfer done\n");
    }


    N = myData[0];
    NN = NUM_PER_BLOCK_AXIS * N;
    for(i = myData[1] * NUM_PER_BLOCK_AXIS; i < myData[2] * NUM_PER_BLOCK_AXIS; i++){
        internalIndex = i % 4;
        blockNum = i / 4;
        blockNum2 = (N + blockNum + 1) % N;
        blockNum3 = (N + blockNum - 1) % N;
        indicies[0] = 4*blockNum  + ((4 + internalIndex) % 4);
        indicies[1] = 4*blockNum  + ((4 + internalIndex + 1) % 4);
        indicies[2] = 4*blockNum  + ((4 + internalIndex - 1) % 4);
        indicies[3] = 4*blockNum2 + ((4 + internalIndex) % 4);
        indicies[4] = 4*blockNum3 + ((4 + internalIndex) % 4);
        // if (DEBUG){
        //     printf("%d %d %d %d - %d  %d  %d  %d  %d\n",internalIndex, blockNum, blockNum2, blockNum3, indicies[0], indicies[1], indicies[2], indicies[3], indicies[4]);
        // }

        q[i] = 0;
        q[i] += 4*mu;
        for(j = 0; j < 5; j++){
            q[i] += -1*u[indicies[j]];
        }
    }

    MPI_Reduce(q, qFinal, NN, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    if(my_rank == 0 && DEBUG){
        printf("Vector transfer done\n");
    }

    printf("Output q: ");
    printVector(q, NN);
    if(my_rank == 0){
        printf("qFinal: ");
        printVector(qFinal, NN);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    // free(q);
    // free(u);
    // free(myData);
    // if(my_rank == 0){
    //     free(qFinal);
    // }
    
    MPI_Finalize();
    return 0;
}