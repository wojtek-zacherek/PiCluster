#include <stdio.h>
#include <mpi.h>

#define DEBUG 1

int countNumPrimes(int start, int end){
    int i, j;
    int totalPrimes = 0;
    int flag;
    for(i = start; i <= end; i++){  
        if(i != 0 && i != 1){
            flag = 0;
            for(j = 2; j <= i/2; j++){
                if(i % j == 0){
                    flag = 1;
                    break;
                }
            }
            if(flag == 0){
                totalPrimes++;
            }
        }else{

        }
    }
    return totalPrimes;
}

int main(char argc, char *argv[]){
    int my_rank, comm_sz;
    int maxNum, totalPrime;

    if(argc <= 1){
        printf("Missing argument: max number\n");
        return 0;
    }

    maxNum = atoi(argv[1]);

    MPI_Init(NULL,NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);

    double startIdx = 0 + ((double)my_rank/comm_sz);
    double endIdx = (((double)my_rank + 1)/comm_sz);
    int delT = maxNum - 0;

    int start = startIdx * delT;
    int end = endIdx * delT - 1;

    if(my_rank == comm_sz - 1){
        end = maxNum;
    }
    if(DEBUG)
        printf("Processor %d running on domain %d to %d\n",my_rank,start,end);
    MPI_Barrier(MPI_COMM_WORLD);
    int result = countNumPrimes(start,end);
    if(DEBUG)
        printf("Processor %d found %d primes\n",my_rank,result);

    MPI_Reduce(&result , &totalPrime , 1, MPI_INT , MPI_SUM , 0 , MPI_COMM_WORLD );
    MPI_Barrier(MPI_COMM_WORLD);
    if(my_rank == 0){
        printf("Total number of prime numbers = %d\n",totalPrime);
    }

    MPI_Finalize();
    return 0;
}
