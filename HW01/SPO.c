#include <stdio.h>
#include <mpi.h>

int main(int argc, char *argv[]){
    int my_rank, comm_sz, ready = 1;

    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);


    if(my_rank > 0){
        MPI_Recv(&ready, 1, MPI_INT, my_rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("Proc %d of %d > Does anyone have a toothpick?\n",my_rank + 1, comm_sz);
    }else{
        printf("Proc %d of %d > Does anyone have a toothpick?\n",my_rank + 1, comm_sz);
    }
    if(my_rank < comm_sz - 1){
        MPI_Send(&ready, 1, MPI_INT, my_rank + 1, 0, MPI_COMM_WORLD);
    }
    
    MPI_Finalize();
    return 0;

}
