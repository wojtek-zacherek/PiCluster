#include <stdio.h>
#include <stdlib.h>

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
void printSubMatrix(int **A, int nCols, int nRows){
    int i,j;
    for(i = 0; i < nRows; i++){
        for(j = 0; j < nCols; j++){
            printf("%u ",A[i][j]);
        }
        printf("\n");
    }
}
void printSubMatrixFlip(int **A, int nCols, int nRows){
    int i,j;
    for(i = 0; i < nCols; i++){
        for(j = 0; j < nRows; j++){
            printf("%u ",A[j][i]);
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


void copyArray(int *src, int *dest, int N){
    int i;
    for(i = 0; i < N; i++){
        dest[i] = src[i];
    }
}
