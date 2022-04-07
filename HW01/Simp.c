// qsub -I -l nodes=1:ppn=10
// mpirun Simp.bin

#include <math.h>
#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>

#define DEBUG 0

enum functionSelect{
    LINEAR = 0,
    QUAD = 1,
    CUBIC = 2,
    SIN = 3,
    COS = 4
};
enum functionSelect FS;

double f(double x){
    // return x*x;
    switch(FS){
        case LINEAR:
            return x;
        case QUAD:
            return x*x;
        case CUBIC:
            return x*x*x;
        case SIN:
            return sin(x);
        case COS:
            return cos(x);
        default:
            return 1;
    };
    
}

double Simpson(double begin, double end, int numPairs, int rank, int comm_sz){
    double my_sum = 0.0;
    double delX = 0;
    double rate = 0;


    rate = (end - begin)/numPairs;
    delX = rate/2.0;
    int i;
    
    if(rank < comm_sz){

    }

    if(numPairs == 0){
        double intermediateValue;
        if(rank == 0){
            intermediateValue = f(begin);
            MPI_Send(&intermediateValue, 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
        }else{
            MPI_Recv(&intermediateValue, 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Send(&intermediateValue, 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
        }
    }else{

        

        /* Handle last value computation
            All cores will compute their own last value.
            All cores (except first) will wait for computed value from previous core.
        */
        double values[2*(int)numPairs + 1];
        
        // Handle first value computation
        values[2 * numPairs] = f(begin + (2 * numPairs)*rate/2.0);
        double toSend = values[2 * numPairs];
        double tempVar;

        if(rank == 0){
            values[0] = f(begin);
            MPI_Send(&toSend, 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
        }else{
            MPI_Recv(&tempVar, 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            values[0] = tempVar;
            if(rank < comm_sz - 1){
                MPI_Send(&toSend, 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
            }
        }
        // Handle intermediate values
        for(i = 1; i < 2*numPairs; i++){
            values[i] = f(begin + i*rate/2.0);
        }
        
        for(i = 0; i < 2*numPairs + 1; i++){
            if(i == 0){
                my_sum += values[i];
            }else if(i == 2*numPairs){
                my_sum += values[i];
            }else{
                if(i % 2 ==1){
                    my_sum += 4*values[i];
                }else{
                    my_sum += 2*values[i];
                }
            }
        }
        my_sum = delX * my_sum / 3.0;
    }
    
    if(DEBUG){
        printf("mysum = %f\n",my_sum);
    }
    return my_sum;
}

int main(int argc, char *argv[]){
    int my_rank, comm_sz;
    double value, total_sum = 0;
    double local_a, local_b;

    // Setting default parameters
    int n_segmentPairs = 10;
    double xBegin = 0;
    double xEnd = 3;
    FS = LINEAR;

    // Simple argument parser
    // mpirun Simp.bin beginning_X_value end_X_value number_regions
    // mpirun Simp.bin 2 8.24 143
    if(argc >= 4){
        xBegin = atof(argv[1]);
        xEnd = atof(argv[2]);
        n_segmentPairs = atoi(argv[3]);
        if(argc >= 5){
            FS = atoi(argv[4]);
        }
    }

    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    if(my_rank == 0){
        printf("xBegin = %f xEnd = %f SegPairs = %d FS = %d\n",xBegin, xEnd, n_segmentPairs, FS);
    }
    
    local_a = (my_rank/(0.0 + comm_sz))*n_segmentPairs;          // Seg DP
    local_b = ((my_rank + 1)/(0.0 + comm_sz))*n_segmentPairs;    // Seg DP


    // Cover edge case inaccuricies that arise from rounding.
    // Transalte the estimated division into discrete segments that align with the number of segments that user wants.
    if(my_rank == 0){
        local_a = 0;
    }else{
        local_a = (int)local_a;
    }
    if(my_rank == comm_sz - 1){
        local_b = n_segmentPairs;
    }else{
        local_b = (int)local_b;
    }
    int numPairs = local_b-local_a;

    // Convert from segment boundaries to corresponding x-values.
    local_a = (local_a / n_segmentPairs) * (xEnd - xBegin) + xBegin;
    local_b = (local_b / n_segmentPairs) * (xEnd - xBegin) + xBegin;

    if(DEBUG){
        printf("Proc %d of %d: Simpson calculation from %.3f to %.3f, with %d segment pairs\n",my_rank + 1, comm_sz,local_a, local_b, numPairs);
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    value = Simpson(local_a, local_b, numPairs, my_rank, comm_sz);    

    if( my_rank != 0){
        MPI_Send(&value, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }else{
        int i;
        total_sum = value;
        for(i = 1; i < comm_sz;i++){
            MPI_Recv(&value, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            total_sum += value;
        }
        
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if( my_rank == 0){
        printf("Area Summation = %0.3f\n",total_sum);
    }
    MPI_Finalize();

    return 0;
}
