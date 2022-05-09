#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>
#include <math.h>

#define DEBUG 1
// mpirun -np 4 --mca btl_vader_single_copy_mechanism none GEPP.bin

void printVectorDouble(double *v1, int N){
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
void printSubMatrix(double *A, int nCols, int nRows){
    int i,j;
    for(i = 0; i < nRows; i++){
        for(j = 0; j < nCols; j++){
            printf("%.3f ",A[i * nCols + j]);
        }
        printf("\n");
    }
}



void populateArray(double *target, int N){
    int i;
    int min = 0;
    int max = 10;
    for(i = 0; i < N; i++){
        target[i] = (double)((rand() % (-1 + (max-min))) + (double)rand()/(double)RAND_MAX);
    }
}
void populateMatrix(double *target, int N){
    int i, j;
    int min = 0;
    int max = 10;
    for(i = 0; i < N; i++){
        for(j = 0; j < N; j++){
            target[i*N + j] = (double)((rand() % (-1 + (max-min))) + (double)rand()/(double)RAND_MAX);
        }
    }
}

void mallocMatrixInt(int **target, int N){
    int i;
    target = (int **)malloc(N * sizeof(int*));
    for(i = 0; i < N; i++){
        target[i] = (int *)malloc(N * sizeof(int));
    }
}
void mallocMatrixDouble(double **target, int N){
    int i;
    target = (double **)malloc(N * sizeof(double*));
    for(i = 0; i < N; i++){
        target[i] = (double *)malloc(N * sizeof(double));
    }
}
void freeMatrixInt(int **target, int N){
    int i;
    for(i = 0; i < N; i++){
        free(target[i]);
    }
    free(target);
}
void freeMatrixDouble(double **target, int N){
    int i;
    for(i = 0; i < N; i++){
        free(target[i]);
    }
    free(target);
}

double getMinDouble(double a, double b){
    return (a > b) ? b : a;
}
double getMaxDouble(double a, double b){
    return (a > b) ? a : b;
}

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
    int *globalRowOrderList;      // N_rows
    int *localRowList;      // N        <-- only row
    double **incomingRow;   // N + 1    <-- row and b vector entry
    int N, N_rows, N_extra;
    int i, j, k;
    double alpha;
    int rowsPerProc[comm_sz];
    int maxRowNum;
    

    // rowsPerProc = malloc(comm_sz * sizeof(int));
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
        N = 4;
        double *matrix;
        

        // TODO: mallocMatrix doesnt work for some odd reason, it worked in all my previous programs
        // mallocMatrixDouble(matrix,N);
        // matrix = (double **)malloc(N * sizeof(double*));
        // for(i = 0; i < N; i++){
        //     matrix[i] = (double *)malloc(N * sizeof(double));
        // }
        matrix = (double *)malloc(N * N * sizeof(double));
        globalRowOrderList = (int *)malloc(N * sizeof(int));

        b = malloc(N * sizeof(double));
        populateMatrix(matrix,N);
        populateArray(b,N);
        // printVectorDouble(b,N);
        printSubMatrix(matrix,N,N);
        printf("matrix %.4f\n",matrix[1]);
        // printVectorDouble(matrix,N*N);

        // Send N;
        MPI_Bcast(&N,1,MPI_INT,0,MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        doSomething(my_rank, comm_sz, N, rowsPerProc, &N_rows, &maxRowNum, procStillInPlay); // Calculate # of rows for each processors
        printf("Proc %d N = %d\nProc %d #rows = %d\n",my_rank,N,my_rank,N_rows);
        printVectorInt(procStillInPlay,comm_sz);
        printf("%d %d %d %d\n",rowsPerProc[0],rowsPerProc[1],rowsPerProc[2],rowsPerProc[3]);

        int sendOffset = N*N_rows;
        for (j = 1; j < comm_sz; j++){
            // MPI_Send(matrix + sendOffset, N*rowsPerProc[j], MPI_DOUBLE, j, 0, MPI_COMM_WORLD);
            MPI_Send(matrix + sendOffset, N*rowsPerProc[j], MPI_DOUBLE, j, 0, MPI_COMM_WORLD);
            sendOffset += N*rowsPerProc[j];
        }
        localRowList = malloc(N_rows * sizeof(int));
        localRowList[0] = 0;
        for(i = 1; i < N_rows; i++){
            localRowList[i] = i + localRowList[0];
        }
        printVectorInt(localRowList,N_rows);
        

        submatrix = matrix;
        
        internalRow = 0;
        // printSubMatrix(submatrix,N,N_rows);
        while(internalCol < N){
            
            // for(internalRow = 0; internalRow < ; internalRow++){
            //     for(i = internalRow; i < N; i++){
            //         submatrix[i * N + internalCol];
            //     }
            // }
            int firstSet = 0;
            double tempVal;
            
            for(i = 0; i < N_rows; i++){
                if(localRowList[i] != -1){
                    tempVal = fabs(submatrix[i * N + internalCol]);
                    printf("%.4f\n",tempVal);
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
            // printf(" <><><> Proc %d pass %d Minimum Value %.4f at row %d\n",my_rank,internalCol,colMinValue,rowMinIndex);

            int numProcToCount = 0;
            if(firstSet == 1){
                compareListValues[0] = colMinValue;
                compareListIndex[0] = rowMinIndex;
                compareListProc[0] = 0;
                numProcToCount = 1;
            }
            for(j = 1; j < comm_sz; j++){
                if(procStillInPlay[j] == 1){
                    MPI_Recv(compareListValues + numProcToCount, 1, MPI_DOUBLE, j, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Recv(compareListIndex + numProcToCount, 1, MPI_INT, j, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Recv(compareListProc + numProcToCount, 1, MPI_INT, j, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    numProcToCount++;
                }
            }
            // printVectorDouble(compareListValues,numProcToCount);
            
            double minVal;
            int minIndex;
            int minProc = -1;
            
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
            printf("..........Col %d Global Min of %.3f at row %d in proc %d\n",internalCol, minVal,minIndex,minProc);
            for(j = 1; j < comm_sz; j++){
                if(procStillInPlay[j] == 1){
                    MPI_Send(&minProc, 1, MPI_INT, j, 0, MPI_COMM_WORLD);
                }
            }
            if(minProc == my_rank){
                numRowsCompleted++;
                for(j = 0; j < N_rows; j++){
                    if(localRowList[j] == rowMinIndex){
                        localRowList[j] = -1;
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
            }else{
                int tempPlay;
                MPI_Recv(&tempPlay, 1, MPI_INT, minProc, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                procStillInPlay[minProc] = tempPlay;
            }
            printVectorInt(procStillInPlay,comm_sz);

            

            internalCol++;
            // MPI_Barrier(MPI_COMM_WORLD);
        }

        
        free(globalRowOrderList);
        free(localRowList);
        free(matrix);
        free(b);
        printf("--->Proc %d done\n",my_rank);
    }else{

        MPI_Bcast(&N,1,MPI_INT,0,MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        doSomething(my_rank, comm_sz, N, rowsPerProc, &N_rows, &maxRowNum, procStillInPlay); // Calculate # of rows for each processors
        // printf("Proc %d N = %d\nProc %d #rows = %d\n",my_rank,N,my_rank,N_rows);

        if(N_rows >= 1){
            submatrix = (double *)malloc(N_rows*N*sizeof(double));
            globalRowOrderList = (int *)malloc(N * sizeof(int));

            
            MPI_Recv(submatrix, N*N_rows, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            localRowList = malloc(N_rows * sizeof(int));
            localRowList[0] = 0;
            for(i = 0; i < my_rank; i++){
                localRowList[0] += rowsPerProc[i];
            }
            for(i = 1; i < N_rows; i++){
                localRowList[i] = i + localRowList[0];
            }
            
            // printVectorInt(localRowList,N_rows);
            // printSubMatrix(submatrix,N,N_rows);
            


            internalCol = 0;
            while(internalCol < N){
                int firstSet = 0;
                double tempVal;
                for(i = 0; i < N_rows; i++){
                    if(localRowList[i] != -1){
                        tempVal = fabs(submatrix[i * N + internalCol]);
                        
                        if(firstSet != 0){
                            if(colMinValue > tempVal){
                                colMinValue = tempVal;
                                rowMinIndex = localRowList[i];
                            }
                        }else{
                            // printf("%d %.4f\n", my_rank, tempVal);
                            colMinValue = fabs(submatrix[i * N + internalCol]);
                            rowMinIndex = localRowList[i];
                            firstSet = 1;
                        }
                    }
                }
                // printf(" <><><> Proc %d Minimum Value %.4f at row %d\n",my_rank,colMinValue,rowMinIndex);

                if(firstSet == 1){
                MPI_Send(&colMinValue, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
                MPI_Send(&rowMinIndex, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
                MPI_Send(&my_rank, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
                }

                int minProc;
                MPI_Recv(&minProc, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                // if(minProc == my_rank){
                //     for(j = 0; j < N_rows; j++){
                //         if(localRowList[j] == rowMinIndex){
                //             localRowList[j] = -1;
                //             break;
                //         }
                //     }
                // }
                printf("---minProc %d\n",minProc);
                if(minProc == my_rank){
                    numRowsCompleted++;
                    printf("Proc 1 match, num rowcomp %d, numrows %d\n",numRowsCompleted,N_rows);
                    for(j = 0; j < N_rows; j++){
                        if(localRowList[j] == rowMinIndex){
                            localRowList[j] = -1;
                            printf("Proc 1 found row %d\n",j);
                            break;
                        }
                    }

                    int inPlay = 1;
                    if(numRowsCompleted >= N_rows){
                        inPlay = 0;
                        printf("Proc 1 out of play\n");
                    }
                    for(j = 0; j < comm_sz; j++){
                        if(j != my_rank){
                            MPI_Send(&inPlay, 1, MPI_INT, j, 0, MPI_COMM_WORLD);
                        }
                    }
                }else{
                    int tempPlay;
                    MPI_Recv(&tempPlay, 1, MPI_INT, minProc, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                    procStillInPlay[minProc] = tempPlay;
                }

                internalCol++;
            }









            free(submatrix);
            
            free(globalRowOrderList);
            free(localRowList);
        }
    }
    
    
    
    

    
    // Calculate # of rows for each processors

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