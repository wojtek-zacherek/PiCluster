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
int populateArray(float *target, int N, int bins, int Nperbin){
    int i, j;
    int min = 0;
    int max = 1000000;
    int tracker[bins];
    for(i = 0; i < bins; i++){
        tracker[i] = 0;
    }
    for(i = 0; i < bins; i++){
        for(j = 0; j < Nperbin; j++){
            target[i * Nperbin + j] = (rand() % (-1 + max/bins)) + (i * ((max-min) / bins)) + (float)rand()/RAND_MAX;
        }
        // target[i] = (rand() % 1000) + (rand() / RAND_MAX);
        
    }
}

// https://linuxhint.com/min-function-c/
int getMin(int a, int b){
    return (a > b) ? b : a;
}
int getMax(int a, int b){
    return (a > b) ? a : b;
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

void copyArray(double *src, double *dest, int N){
    int i;
    for(i = 0; i < N; i++){
        dest[i] = src[i];
    }
}

int cmpfunc (const void * a, const void * b) {
   return ( *(int*)a - *(int*)b );
}

int checkSortiness(float *target, int N){
    int i;
    for(i = 0; i < N-1; i++){
        if(target[i] > target[i+1]){
            printf("Sorting check error between indicies %d & %d\n",i,i+1);
            return 1;
        }
    }
    printf("Sort verified!\n");
    return 0;
}

int main(int argc, char *argv[]){
    int my_rank, comm_sz;
    
    MPI_Init(NULL,NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    srand(time(NULL)+ my_rank);

    if(my_rank==0){
        printf("Important State Variables:\n\tnumber of proccessors: %d\n\tcurrent processor rank: %d\n",comm_sz,my_rank);
    }
    // Assumptions: we know the number of elements. Buckets are same size. Data is is divided into buckets when created.

    int i;
    int numBins = comm_sz;
    int numElements = numBins * 10000000;
    float *arrayToSort;
    void **bins;
    float min, max;
    int searchNum = getMax(10, 0.1*numElements);
    float elementsPer = numElements / numBins;
    int elementsExtra = numElements % numBins;
    float *localBucket;
    if(my_rank==0){
        printf("Program Variables:\n\tnumBins: %d\n\tnumElements: %d\n\telementsPer:%.4f\n",numBins,numElements,elementsPer);
    }
    
    if(my_rank == 0){
        printf("Populating array with %d floats\n",numElements);
        arrayToSort = (float *)malloc(numElements * sizeof(float));
        populateArray(arrayToSort, numElements, numBins, numElements/numBins);
        // printVector(arrayToSort,numElements);
        findMinMax(arrayToSort,numElements,&min,&max,searchNum);
        printf("Local Min = %f \nLocal Max = %f \n",min,max);    
    }
    
    localBucket = (float *)malloc(elementsPer * sizeof(float));
    MPI_Scatter(arrayToSort,elementsPer,MPI_FLOAT,localBucket,elementsPer,MPI_FLOAT,0,MPI_COMM_WORLD);
    qsort(localBucket, elementsPer, sizeof(float), cmpfunc);
    MPI_Gather(localBucket,elementsPer,MPI_FLOAT,arrayToSort,elementsPer,MPI_FLOAT,0,MPI_COMM_WORLD);
    

    if(my_rank == 0){
        // printVector(arrayToSort,numElements);
        checkSortiness(arrayToSort,numElements);
        free(arrayToSort);
    }
    free(localBucket);
    MPI_Finalize();
    return 0;
}