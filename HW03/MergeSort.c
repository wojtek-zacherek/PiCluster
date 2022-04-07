#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>

#define DEBUG 1

// https://stackoverflow.com/questions/1787996/c-library-function-to-perform-sort

int compare_function(const void *a,const void *b) {
double *x = (double *) a;
double *y = (double *) b;
// return *x - *y; // this is WRONG...
if (*x < *y) return -1;
else if (*x > *y) return 1; return 0;
}

int main(int argc, char *argv[]){
    int my_rank, comm_sz;
    
    MPI_Init(NULL,NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    srand(time(NULL) + my_rank);

    // Assuming n % comm_sz = 0
    int n;
    int alpha = 4;
    if(argc == 2){
        n = atoi(argv[1]);
    }else{
        n = alpha * comm_sz;
    }

    MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);

    double *localList, *beginningList;;
    int local_n;

    local_n = n / comm_sz;
    localList = malloc(local_n * sizeof(double));
    int i, j;
    for(i = 0; i < local_n; i++){
        localList[i] = rand() % 100;
    }

    if(my_rank == 0){
        printf("n = %d,   local_n = %d\n",n,local_n);
    }

    
    if(my_rank == 0){
        beginningList = malloc(n * sizeof(double));
    }

    qsort(localList, local_n, sizeof(double), compare_function);
    for(j = 0; j < local_n; j++){
        printf("%.1f  ",localList[j]);
    }
    printf("\n");

    MPI_Gather(localList, local_n, MPI_DOUBLE,
                beginningList, local_n, MPI_DOUBLE,
                0,MPI_COMM_WORLD);
    
    if(my_rank == 0){
        for(i = 0; i < comm_sz; i++){
            for(j = 0; j < local_n; j++){
                printf("%.1f  ",beginningList[i * local_n + j]);
            }
        }
        printf("\n");
    }

    

    printf("Beginning Merge\n");
    MPI_Barrier(MPI_COMM_WORLD);
    int divisor = 1;
    int amount = 1;
    double *localListNew;
    double *ptrTemp;


    while(divisor <= comm_sz){
        divisor = 2 * divisor;

        if(my_rank % divisor != 0){
            MPI_Send(&local_n, 1, MPI_INT, my_rank - (divisor/2), 0, MPI_COMM_WORLD);
            MPI_Send(localList, local_n, MPI_DOUBLE, my_rank - (divisor/2), 0, MPI_COMM_WORLD);
            printf("Goodbye %d\n",my_rank);
            break;
        }else{
            // If there is a tree-neighbor to pull from, do it. Otherwise, do nothing and iterate.
            if(my_rank==0){
                printf("Hello\n");
            }
            if(my_rank + divisor/2 < comm_sz){
                int incomingSize = 0;
                MPI_Recv(&incomingSize, 1, MPI_INT, my_rank + (divisor/2), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                localListNew = malloc((local_n + incomingSize)  * sizeof(double));
                for(i = 0; i < local_n; i++){
                    localListNew[i] = localList[i];
                }
                MPI_Recv(localListNew + local_n, incomingSize, MPI_DOUBLE, my_rank + (divisor/2), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                free(localList);        // relieve old data
                localList = localListNew;
                local_n = local_n + incomingSize;
                qsort(localList, local_n, sizeof(double), compare_function);
            }else{
                
            }
        }
        
    }



    MPI_Barrier(MPI_COMM_WORLD);
    if(my_rank == 0){
        printf("\n--Begin Results--\n");
        printf("\n");
        // for(i = 0; i < comm_sz; i++){
            for(j = 0; j < local_n; j++){
                printf("%.1f  ",localList[j]);
            }
            printf("\n");
        // }
        printf("\n--End Results--\n");
    
    }
    
    if(my_rank == 0){
        free(beginningList);
    }
    free(localList);

    MPI_Finalize();
    return 0;
}
