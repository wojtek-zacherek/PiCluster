// Gaussian Elimination with Partial Pivoting Draft
// Wojciech Zacherek
// Used for running in VSCode with MPI plugin: mpirun -np 4 --mca btl_vader_single_copy_mechanism none GEPP.bin

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>
#include <math.h>
#include "UsefulFunction.h"

#define DEBUG 1
char response = 'n';
int N_user = 3;
double matInputs[] = {1,1,1,1,1,2,1,1,2};
double bInputs[] = {1,3,-1};

void doSomething(int my_rank, int comm_sz, int N, int *rowsPerProc, int *N_rows, int *maxRowNum, int *procStillInPlay){
    // Calculate # of rows for each processors
    *N_rows = (int)(N * ((double)1/comm_sz));
    int N_extra = N - comm_sz*(*N_rows);
    maxRowNum = 0;
    int i;
    for(i = 0; i < comm_sz; i++){
        rowsPerProc[i] = *N_rows;
        if(i < N_extra){
            rowsPerProc[i]++;
        }
        if(i <= my_rank){
            maxRowNum += rowsPerProc[i];
        }
        if(rowsPerProc[i] > 0){
            procStillInPlay[i] = 1;
        }else{
            procStillInPlay[i] = 0;
        }
    }
    *N_rows = rowsPerProc[my_rank];
    
}

void kernel(int my_rank, int comm_sz){
    
    // Necesary data
    double *submatrix;
    double *b;
    double *subB;
    int *globalRowOrderList;      // N_rows
    int *localRowList;      // N        <-- only row
    double *incomingRow;   // N + 1    <-- row and b vector entry
    int N, N_rows, N_extra;
    int i, j, k;
    double alpha;
    int rowsPerProc[comm_sz];
    int maxRowNum;
    

    double *matrix;
    int internalRow = 0, internalCol = 0;
    double compareListValues[comm_sz];
    int compareListIndex[comm_sz];
    int compareListProc[comm_sz];
    double colMinValue;
    int rowMinIndex;
    int procStillInPlay[comm_sz]; // 1 = yes, 0 = no
    int numRowsCompleted = 0;
    // Input N, matrix 
    if(my_rank == 0){

        
        //// User input matrix dimension N OR use default
        if(response == 'y' || response == 'Y'){
            // TODO
            N = N_user;
        }else{
            N = 3;
        }
        ////////////////////////////////
        MPI_Bcast(&N,1,MPI_INT,0,MPI_COMM_WORLD);
        ////////////////////////////////
        doSomething(my_rank, comm_sz, N, rowsPerProc, &N_rows, &maxRowNum, procStillInPlay); // Calculate # of rows for each processors

        //// User input matrix OR generate
        matrix = (double *)malloc(N * N * sizeof(double));
        if(response == 'y' || response == 'Y'){
            copyArray(matInputs,matrix,N*N);
        }else{
            for(i = 0; i < N*N; i++){
                matrix[i] = i + 1;
            }
            matrix[0] = 2;  // Make matrix nonsingular
        }
        
        
        //// User input b OR generate
        b = malloc(N * sizeof(double));
        if(response == 'y' || response == 'Y'){
            copyArray(bInputs,b,N);
        }else{
            
            for(i = 0; i < N; i++){
                b[i] = i + 1;
            }
        }

        printf("Matrix = \n");
        printSubMatrix(matrix,N,N);
        printf("Vector b = ");
        printVectorDouble(b,N);
    }else{
        ////////////////////////////////
        MPI_Bcast(&N,1,MPI_INT,0,MPI_COMM_WORLD);
        ////////////////////////////////
        doSomething(my_rank, comm_sz, N, rowsPerProc, &N_rows, &maxRowNum, procStillInPlay); // Calculate # of rows for each processors
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if(N_rows >= 1){
        globalRowOrderList = (int *)malloc(N * sizeof(int));
        incomingRow = (double *)malloc((N+1) * sizeof(double));
        submatrix = (double *)malloc(N_rows*N*sizeof(double));
        subB = (double *)malloc(N_rows*sizeof(double));
        if(my_rank == 0){
            copyArray(matrix,submatrix,N*N_rows);
            copyArray(b,subB,N_rows);
            
            /////// PROC 0 DISTRIBUTION BEGIN
            int sendOffset = N_rows;
            for (j = 1; j < comm_sz; j++){
                MPI_Send(matrix + N*sendOffset, N*rowsPerProc[j], MPI_DOUBLE, j, 0, MPI_COMM_WORLD);
                MPI_Send(b + sendOffset, rowsPerProc[j], MPI_DOUBLE, j, 0, MPI_COMM_WORLD);
                sendOffset += rowsPerProc[j];
            }
            /////// PROC 0 DISTRIBUTION EMD

        }else{
            // b = (double *)malloc(N_rows*sizeof(double));

            /////// SLAVE DISTRIBUTION BEGIN
            MPI_Recv(submatrix, N*N_rows, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(subB, N_rows, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            /////// SLAVE DISTRIBUTION END
        }

        localRowList = malloc(N_rows * sizeof(int));
        localRowList[0] = 0;
        for(i = 0; i < my_rank; i++){
            localRowList[0] += rowsPerProc[i];
        }
        for(i = 1; i < N_rows; i++){
            localRowList[i] = i + localRowList[0];
        }
        for(i = 0; i < N; i++){
            globalRowOrderList[i] = i;
        }

        internalRow = 0;
        while(internalCol < N){
            //////////////
            int firstSet = 0;
            double tempVal;
            
            for(i = 0; i < N_rows; i++){
                if(localRowList[i] != -1){
                    // if(my_rank==1){printSubMatrix(submatrix,N,N_rows);}
                    tempVal = fabs(submatrix[i * N + internalCol]);
                    // printf("   inter %d    proc %d    value %.4f\n",internalCol,my_rank,tempVal);
                    if(firstSet != 0){
                        if(colMinValue > tempVal){
                            colMinValue = tempVal;
                            rowMinIndex = localRowList[i];
                        }
                    }else{
                        colMinValue = fabs(submatrix[i * N + internalCol]);
                        rowMinIndex = localRowList[i];
                        firstSet = 1;
                    }
                }
            }


            double minVal;
            int minIndex;
            int minProc = -1;
            if(my_rank == 0){
                int numProcToCount = 0;
                if(firstSet == 1){
                    compareListValues[0] = colMinValue;
                    compareListIndex[0] = rowMinIndex;
                    compareListProc[0] = 0;
                    numProcToCount = 1;
                }


                /////// PROC 0 BAREBONE BEGIN
                for(j = 1; j < comm_sz; j++){
                    if(procStillInPlay[j] == 1){
                        MPI_Recv(compareListValues + numProcToCount, 1, MPI_DOUBLE, j, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        MPI_Recv(compareListIndex + numProcToCount, 1, MPI_INT, j, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        MPI_Recv(compareListProc + numProcToCount, 1, MPI_INT, j, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        numProcToCount++;
                    }
                }
                
                firstSet = 0;
                for(j = 0; j < numProcToCount; j++){
                    if(firstSet != 0){
                        tempVal = compareListValues[j];
                        if(minVal > tempVal){
                            minVal = compareListValues[j];
                            minIndex = compareListIndex[j];    
                            minProc = compareListProc[j];
                        }
                    }else{
                        minVal = compareListValues[j];
                        minIndex = compareListIndex[j];
                        minProc = compareListProc[j];
                        firstSet = 1;
                    }
                }
                if(minVal == 0){
                    // SKIP Iteration somehow
                }
                printf("..........Col %d Global Min of %.3f at row %d in proc %d\n",internalCol, minVal,minIndex,minProc);
                for(j = 1; j < comm_sz; j++){
                    if(procStillInPlay[j] == 1){
                        MPI_Send(&minProc, 1, MPI_INT, j, 0, MPI_COMM_WORLD);
                        MPI_Send(&minIndex, 1, MPI_INT, j, 0, MPI_COMM_WORLD);
                        MPI_Send(&minVal, 1, MPI_DOUBLE, j, 0, MPI_COMM_WORLD);
                    }
                }
                /////// PROC 0 BAREBONE END
            }else{
                if(firstSet == 1){
                    MPI_Send(&colMinValue, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
                    MPI_Send(&rowMinIndex, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
                    MPI_Send(&my_rank, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);

                    MPI_Recv(&minProc, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Recv(&minIndex, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Recv(&minVal, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
            
            if(minVal > 0){
            if(globalRowOrderList[internalRow] != minIndex){
                int temp = globalRowOrderList[internalRow];
                globalRowOrderList[internalRow] = minIndex;
                globalRowOrderList[minIndex] = temp;
            }
            internalRow++;
            if(minProc == my_rank){
                numRowsCompleted++;
                int tempJForCopy;
                for(j = 0; j < N_rows; j++){
                    if(localRowList[j] == rowMinIndex){
                        localRowList[j] = -1;
                        tempJForCopy = j;
                        break;
                    }
                }


                int inPlay = 1;
                if(numRowsCompleted >= N_rows){
                    inPlay = 0;
                }
                for(j = 0; j < comm_sz; j++){
                    if(j != my_rank){
                        MPI_Send(&inPlay, 1, MPI_INT, j, 0, MPI_COMM_WORLD);
                    }
                }
                copyArray(submatrix + N*tempJForCopy, incomingRow,N);
                incomingRow[N] = subB[tempJForCopy];
                for(j = 0; j < comm_sz; j++){
                    if(j != my_rank){
                        MPI_Send(incomingRow, N+1, MPI_DOUBLE, j, 0, MPI_COMM_WORLD);
                    }
                }
            }else{
                int tempPlay;
                MPI_Recv(&tempPlay, 1, MPI_INT, minProc, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                MPI_Recv(incomingRow, N+1, MPI_DOUBLE, minProc, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                procStillInPlay[minProc] = tempPlay;
            }

            for(k = 0; k < N_rows; k++){
                if(localRowList[k] != -1){
                    alpha = submatrix[k * N + internalCol] / incomingRow[internalCol] ;
                    submatrix[k * N + internalCol] = 0;
                    for(j = internalCol + 1; j < N; j++){
                        submatrix[k * N + j] = submatrix[k * N  + j] - alpha*incomingRow[j];
                    }
                    subB[k] = subB[k] - alpha*incomingRow[N];
                }
            }
            }

            // ITERATIVE DEUBGGING OF RESULTS
            if(my_rank == 0){
                int offset = N_rows;
                copyArray(submatrix,matrix,N*N_rows);
                copyArray(subB,b,N_rows);
                for(j = 1; j < comm_sz; j++){
                    if(j < N){
                        MPI_Recv(matrix + N*offset, N*rowsPerProc[j], MPI_DOUBLE, j, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        MPI_Recv(b + offset, rowsPerProc[j], MPI_DOUBLE, j, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        offset += rowsPerProc[j];
                    }
                }
                printf("Matrix reduced %d = \n", internalCol);
                printSubMatrix(matrix,N,N);
                printf("Vector b reduced %d = \n",internalCol);
                printVectorDouble(b,N);
                printf("--------------\n");
            }else{
                MPI_Send(submatrix, N*N_rows, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
                MPI_Send(subB, N_rows, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            }

            internalCol++;
        }

        if(my_rank == 0){
            int offset = N_rows;
            copyArray(submatrix,matrix,N*N_rows);
            for(j = 1; j < comm_sz; j++){
                if(j < N){
                    MPI_Recv(matrix + N*offset, N*rowsPerProc[j], MPI_DOUBLE, j, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Recv(b + offset, rowsPerProc[j], MPI_DOUBLE, j, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    offset += rowsPerProc[j];
                }
            }

            double *bNew = malloc(N*sizeof(double));
            for(j = 0; j < N; j++){
                // bNew[j] = b[globalRowOrderList[j]];
            }
            printf("Matrix reduced = \n");
            printSubMatrix(matrix,N,N);
            printf("Vector b reduced = \n");
            printVectorDouble(b,N);

            free(matrix);
            free(bNew);
            free(b);
        }else{
            MPI_Send(submatrix, N*N_rows, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            MPI_Send(subB, N_rows, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            
        }

        free(submatrix);
        free(globalRowOrderList);
        free(localRowList);
        free(incomingRow);
        free(subB);
    }
}


int main(int argc, char *argv[]){
    int my_rank, comm_sz;
    
    MPI_Init(NULL,NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    srand(time(NULL)+ my_rank);

    kernel(my_rank,comm_sz);

    MPI_Finalize();
    return 0;
}

