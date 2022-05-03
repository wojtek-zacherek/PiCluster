#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>
#include <float.h>

#define DEBUG 1

enum programType{Float,Double,Int};

void printVector(float *v1, int N){
    int i;
    for(i = 0; i < N; i++){
        printf("%.3f ",v1[i]);
    }
    printf("\n");
}
int populateArray(float *target, int N){
    int i;
    for(i = 0; i < N; i++){
        target[i] = (float)(rand());
    }
}

// https://linuxhint.com/min-function-c/
int getMin(int a, int b){
    return (a > b) ? b : a;
}


void findMinMax(float *target, int N, float *min, float *max, int searchNum){
    *min = *max = target[0];
    int start = rand() % (N - searchNum);
    int i;

    printf("Starting search at i=%d for %d values\n",start,searchNum);
    for(i = 0; i < searchNum; i++){
        if(target[i + start] < *min){
            *min = target[i + start];
        }else if(target[i + start] > *max){
            *max = target[i + start];
        }
    }
}


int main(int argc, char *argv[]){
    int my_rank, comm_sz;
    
    MPI_Init(NULL,NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    srand(time(NULL)+ my_rank);

    int numBins = 10;
    int numElements = 30;
    float *arrayToSort;
    void **bins;
    float min, max;
    int searchNum = getMin(10, 0.1*numElements);
    enum programType pType = Float;
    

    printf("Populating array with %d floats\n",numElements);
    arrayToSort = (float *)malloc(numElements * sizeof(float));
    populateArray(arrayToSort, numElements);
    printVector(arrayToSort,numElements);
    findMinMax(arrayToSort,numElements,&min,&max,searchNum);
    printf("Min = %f \nMax = %f \n",min,max);


    free(arrayToSort);
    MPI_Finalize();
    return 0;
}