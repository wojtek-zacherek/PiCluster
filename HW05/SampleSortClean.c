#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>
#include <float.h>
#include <unistd.h>
#include <sys/time.h>

#define DEBUG 0

// -np 4 --mca btl_vader_single_copy_mechanism none

void printVector(float *v1, int N){
    int i;
    for(i = 0; i < N; i++){
        printf("%.3f ",v1[i]);
    }
    printf("\n");
}
void populateArray(float *target, int N){
    int i;
    int min = 0;
    int max = 10;
    for(i = 0; i < N; i++){
        target[i] = (float)((rand() % (-1 + (max-min))) + (float)rand()/(float)RAND_MAX);
    }
}

// https://linuxhint.com/min-function-c/
int getMinInt(int a, int b){
    return (a > b) ? b : a;
}
int getMaxInt(int a, int b){
    return (a > b) ? a : b;
}


void findMinMax(float *target, int N, float *min, float *max, int searchNum){
    *min = *max = target[0];
    int start = rand() % (N - searchNum);
    int i;

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
    struct timeval programStart, programEnd, variableSetupStart, variableSetupEnd, arraySetupStart, arraySetupEnd, stage1Start, stage1End, distributionStart, distributionEnd, stage2Start, stage2End;
    gettimeofday(&programStart,NULL);

    int my_rank, comm_sz;
    int i, j;
    
    MPI_Init(NULL,NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    srand(time(NULL)+ my_rank);
    

    if(my_rank==0){
        printf("Important State Variables:\n\tnumber of proccessors: %d\n\troot processor rank: %d\n",comm_sz,my_rank);
    }
    // Assumptions: we know the number of elements. Buckets are same size. Data is is divided into buckets when created.

    gettimeofday(&variableSetupStart,NULL);
    int NBins = comm_sz;
    int N = 6938503;
    float *targetArray;
    float min, max;
    float padding = 0.05;
    int searchNum = getMaxInt(10, 0.1*N);
    
    float **bucketArray;
    float *localBucket;//
    int runningBucketSizes[comm_sz];//
    int averageBucketSize = N / NBins + 1;
    int myBucketSize = -1;//
    gettimeofday(&variableSetupEnd,NULL);

    
    
    if(my_rank==0){
        printf("Important Program Variables:\n\tNBins: %d\n\tN: %d\n\taverageBucketSize: %d\n\tpadding: %.4f\n",NBins,N,averageBucketSize,padding);
    }
    
    if(my_rank == 0){
        gettimeofday(&arraySetupStart,NULL);
        targetArray = (float *)malloc(N * sizeof(float));
        bucketArray = (float **)malloc(NBins * sizeof(float*));
        for(i = 0; i < NBins; i++){
            bucketArray[i] = (float *)malloc(averageBucketSize * sizeof(float));
        }
        if(DEBUG) printf("Populating array with %d floats\n",N);
        populateArray(targetArray, N);
        findMinMax(targetArray,N,&min,&max,searchNum);
        if(DEBUG) printf("\tApproximate Min = %f\n\tApproximate Max = %f\n",min,max);
        gettimeofday(&arraySetupEnd,NULL);
    }

    
    if(my_rank == 0){
        gettimeofday(&stage1Start,NULL);
        float borders[comm_sz - 1];
        for(i = 1; i < comm_sz; i++){
            borders[i - 1] = i * ((float)(max-min)/(NBins)) + min;
            if(DEBUG) printf("\tborder %d: %f\n",i-1,borders[i - 1]);
        }
        
        int currentSizes[comm_sz];
        for(i = 0; i < comm_sz; i++){
            runningBucketSizes[i] = 0;
            currentSizes[i] = averageBucketSize;
        }
        
        float currEL = 0;
        if(DEBUG) printf("Beginning stage 1 bucket sort\n");
        for(i = 0; i < N; i++){
            currEL = targetArray[i];
            
            for(j = 0; j < comm_sz; j++){
                if((j == comm_sz - 1) || (currEL < borders[j])){
                    if(runningBucketSizes[j] + 1 > currentSizes[j]){
                        currentSizes[j] = currentSizes[j] + getMaxInt(10,(int)(padding*averageBucketSize));
                        bucketArray[j] = realloc(bucketArray[j],(currentSizes[j] * sizeof(float)));
                        if(DEBUG) printf("\tIncreasing bucket %d to %d. Curr num of elements: %d\n",j,currentSizes[j],runningBucketSizes[j]);
                    }
                    bucketArray[j][runningBucketSizes[j]] = currEL;
                    runningBucketSizes[j]++;
                    break;
                }
            }
        }
        if(DEBUG) printf("Complete stage 1 bucket sort\n");
        gettimeofday(&stage1End,NULL);

        gettimeofday(&distributionStart,NULL);
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Scatter(runningBucketSizes,1,MPI_INT,&myBucketSize,1,MPI_INT,0,MPI_COMM_WORLD);
        localBucket = bucketArray[0];
        MPI_Request requests[comm_sz];
        for(i = 1; i < comm_sz; i++){
            MPI_Isend(bucketArray[i],runningBucketSizes[i],MPI_FLOAT,i,0,MPI_COMM_WORLD,&requests[i]);
            // MPI_Send(bucketArray[i],runningBucketSizes[i],MPI_FLOAT,i,0,MPI_COMM_WORLD);
        }
        gettimeofday(&distributionEnd,NULL);

        gettimeofday(&stage2Start,NULL);
        if(DEBUG) printf("Beginning stage 2 bucket sort\n");
        qsort(localBucket, myBucketSize, sizeof(float), compare);
        copyArray(localBucket,targetArray,myBucketSize);
        for(i = 1; i < comm_sz; i++){
            MPI_Wait(&requests[i], MPI_STATUS_IGNORE);
        }

        int recvOffset = myBucketSize;
        for(i = 1; i < comm_sz; i++){
            MPI_Recv(targetArray + recvOffset, runningBucketSizes[i], MPI_FLOAT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            recvOffset += runningBucketSizes[i];
        }

        if(DEBUG) printf("Complete stage 2 bucket sort\n");

    }else{
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Scatter(runningBucketSizes,1,MPI_INT,&myBucketSize,1,MPI_INT,0,MPI_COMM_WORLD);

        localBucket = (float *)malloc(myBucketSize * sizeof(float));
        MPI_Recv(localBucket,myBucketSize,MPI_FLOAT,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        qsort(localBucket, myBucketSize, sizeof(float), compare);
        MPI_Send(localBucket,myBucketSize,MPI_FLOAT,0,0,MPI_COMM_WORLD);

    }
    

    
    

    
    // if(my_rank==0){
    //     printf("Beginning stage 2 bucket sort\n");
    //     localBucket = bucketArray[0];
    //     MPI_Request requests[comm_sz];
    //     for(i = 1; i < comm_sz; i++){
    //         MPI_Isend(bucketArray[i],runningBucketSizes[i],MPI_FLOAT,i,0,MPI_COMM_WORLD,&requests[i]);
    //         // MPI_Send(bucketArray[i],runningBucketSizes[i],MPI_FLOAT,i,0,MPI_COMM_WORLD);
    //     }


    //     qsort(localBucket, myBucketSize, sizeof(float), compare);
    //     copyArray(localBucket,targetArray,myBucketSize);
    //     for(i = 1; i < comm_sz; i++){
    //         MPI_Wait(&requests[i], MPI_STATUS_IGNORE);
    //     }

    //     int recvOffset = myBucketSize;
    //     for(i = 1; i < comm_sz; i++){
    //         MPI_Recv(targetArray + recvOffset, runningBucketSizes[i], MPI_FLOAT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //         recvOffset += runningBucketSizes[i];
    //     }

    //     printf("Complete stage 2 bucket sort\n");
    // }else{
    //     localBucket = (float *)malloc(myBucketSize * sizeof(float));
    //     MPI_Recv(localBucket,myBucketSize,MPI_FLOAT,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    //     qsort(localBucket, myBucketSize, sizeof(float), compare);
    //     MPI_Send(localBucket,myBucketSize,MPI_FLOAT,0,0,MPI_COMM_WORLD);
    // }
    

    if(my_rank == 0){
        checkSortiness(targetArray,N);
        free(targetArray);
        for(i = 0; i < NBins; i++){
            free(bucketArray[i]);
        }
        free(bucketArray);
    }else{
        free(localBucket);
    }
    
    MPI_Finalize();
    gettimeofday(&programEnd,NULL);

    return 0;
}