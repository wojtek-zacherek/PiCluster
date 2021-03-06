#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>
#include <float.h>

#define DEBUG 1

// -np 4 --mca btl_vader_single_copy_mechanism none

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
    int max = 10;
    int tracker[bins];
    // for(i = 0; i < bins; i++){
    //     tracker[i] = 0;
    // }
    for(i = 0; i < N; i++){
        target[i] = (rand()) + ((float)rand() / RAND_MAX);
        target[i] = (float)((rand() % (-1 + max)) + (float)rand()/(float)RAND_MAX);
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

void copyArray(float *src, float *dest, int N){
    int i;
    for(i = 0; i < N; i++){
        dest[i] = src[i];
    }
}

int cmpfunc (const void * a, const void * b) {
   return ( *(float*)a - *(float*)b );
}
int compare (const void * a, const void * b)
{
  float fa = *(const float*) a;
  float fb = *(const float*) b;
  return (fa > fb) - (fa < fb);
}

int checkSortiness(float *target, int N){
    int i;
    for(i = 0; i < N-1; i++){
        if(target[i] > target[i+1]){
            printf("!!!!Sorting check error between indicies %d & %d!!!!\n\t%.2f , %.2f, diff=%.2f\n",i,i+1,target[i],target[i+1],target[i]-target[i+1]);
            return 1;
        }
    }
    printf("!!!!Sort verified!!!!!\n");
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
    int numElements = numBins * 100000000;
    float *arrayToSort;
    void **bins;
    float min, max;
    int searchNum = getMax(10, 0.1*numElements);
    float elementsPer = numElements / numBins;
    int elementsExtra = numElements % numBins;
    float *localBucket;
    float padding = 0.05;
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

    float **arraysBinned;
    int elementsProc[comm_sz];
    if(my_rank == 0){
        arraysBinned = (float **)malloc(numBins * sizeof(float*));
        for(i = 0; i < numBins; i++){
            arraysBinned[i] = (float *)malloc(elementsPer * sizeof(float));
        }


        float borders[comm_sz - 1];
        for(i = 1; i < comm_sz; i++){
            borders[i - 1] = i * ((float)(max-min)/(numBins)) + min;
            printf("  -  %f\n",borders[i - 1]);
        }

        // printVector(borders,comm_sz-1);
        
        int currentSizes[comm_sz];
        for(i = 0; i < comm_sz; i++){
            elementsProc[i] = 0;
            currentSizes[i] = elementsPer;
        }
        int j;
        float currEL = 0;
        for(i = 0; i < numElements; i++){
            currEL = arrayToSort[i];
            
            for(j = 0; j < comm_sz; j++){
                if((j == comm_sz - 1) || (currEL < borders[j])){
                    if(elementsProc[j] + 1 > currentSizes[j]){
                        currentSizes[j] = currentSizes[j] + getMax(10,(int)(padding*elementsPer));
                        arraysBinned[j] = realloc(arraysBinned[j],(currentSizes[j] * sizeof(float)));
                        printf("Increasing bucket %d to %d. Curr num of elements: %d\n",j,currentSizes[j],elementsProc[j]);
                    }
                    arraysBinned[j][elementsProc[j]] = currEL;
                    elementsProc[j]++;
                    break;
                }
            }
            // if(j == comm_sz - 2){
            //     arraysBinned[j+1][elementsProc[j+1]] = currEL;
            //     elementsProc[j+1]++;
            // }

        }
        printf("Done Bucket sorting\n");
        int sum = 0;
        for(i = 0; i < numBins; i++){
            // printVector(arraysBinned[i],elementsProc[i]);
            printf("%d elePro = %d\n",i,elementsProc[i]);
            sum+=elementsProc[i];
        }
        printf("Sum:%d\n",sum);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    int myNum = -1;
    MPI_Scatter(elementsProc,1,MPI_INT,&myNum,1,MPI_INT,0,MPI_COMM_WORLD);
    printf("  My num=%d\n",myNum);

    if(my_rank==0){
        localBucket = arraysBinned[0];
        // MPI_Request requests[comm_sz - 1];
        for(i = 1; i < comm_sz; i++){
            // MPI_Isend(arraysBinned[i],elementsProc[i],MPI_FLOAT,i,0,MPI_COMM_WORLD,requests[i-1]);
            MPI_Send(arraysBinned[i],elementsProc[i],MPI_FLOAT,i,0,MPI_COMM_WORLD);
        }

        // printVector(localBucket,myNum);
        qsort(localBucket, myNum, sizeof(float), compare);
        // printVector(localBucket,myNum);
        // for(i = 1; i < comm_sz; i++){
        //     MPI_Wait(requests[i-1], MPI_STATUS_IGNORE);
        // }
        copyArray(localBucket,arrayToSort,myNum);
        int offset = myNum;
        MPI_Barrier(MPI_COMM_WORLD);
        for(i = 1; i < comm_sz; i++){
            printf("        Received From proc%d\n",i);
            MPI_Recv(arrayToSort + offset, elementsProc[i], MPI_FLOAT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            offset += elementsProc[i];
        }
        printf("        Done recieing\n");
        

    }else{
        
        localBucket = (float *)malloc(myNum * sizeof(float));
        MPI_Recv(localBucket,myNum,MPI_FLOAT,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        printf("Received %d\n",my_rank);
        
        qsort(localBucket, myNum, sizeof(float), compare);
        
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Send(localBucket,myNum,MPI_FLOAT,0,0,MPI_COMM_WORLD);
        printf("Sent %d\n",my_rank);
    }

    
    
    // MPI_Scatter(arrayToSort,elementsPer,MPI_FLOAT,localBucket,elementsPer,MPI_FLOAT,0,MPI_COMM_WORLD);
    // qsort(localBucket, elementsPer, sizeof(float), cmpfunc);
    // MPI_Gather(localBucket,elementsPer,MPI_FLOAT,arrayToSort,elementsPer,MPI_FLOAT,0,MPI_COMM_WORLD);
    

    if(my_rank == 0){
        // printVector(arrayToSort,numElements);
        // printf("DONEDONE\n");
        checkSortiness(arrayToSort,numElements);
        // printf("Freeing arrayToSort\n");
        free(arrayToSort);
        for(i = 0; i < numBins; i++){
            printf("Freeing %d\n",i);
            free(arraysBinned[i]);
        }
        // printf("DONEDONE\n");
        free(arraysBinned);
        // printf("11DONEDONE\n");
        // free(localBucket);
        // printf("22DONEDONE\n");
    }else{
        // printf("LB %d\n",my_rank);
        free(localBucket);
        // printf("Freed %d\n",my_rank);
    }
    
    MPI_Finalize();
    return 0;
}