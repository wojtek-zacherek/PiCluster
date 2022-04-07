#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define DEBUG 0

double f(double);

int type = 0;

double Trap(double left_endpt, double right_endpt,int trap_count,
    double base_len){
    double estimate, x;
    int i;

    estimate = (f(left_endpt) + f(right_endpt))/2.0;
    for(i = 1; i <= trap_count-1; i++){
        x = left_endpt + i*base_len;
        estimate += f(x);
    }
    estimate = estimate*base_len;

    return estimate;
}

double f(double x){
    //return pow(x,2);
    if(type == 0){
        return x*x;
    }
    if(type == 1){
        return x;
    }
    return 1;
}

int main(int argc, char *argv[]){
    int my_rank, comm_sz, n = 1024;
    double a = 0.0, b = 3.0, h, local_a, local_b;
    double local_int, total_int;
    double local_n;
    int source;

    if(argc > 1){
        n = atoi(argv[1]);
        // printf("n = %d\n",n);
        if(argc > 2){
            type = atoi(argv[2]);
            // printf("type = %d\n",type);
        }
    }


    MPI_Init(NULL,NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);

    h = (b-a)/n;
    local_n = ((double)n)/comm_sz;

    // local_a = a + my_rank * local_n * h;
    // local_b = local_a + local_n * h;

    local_a = (my_rank/(0.0 + comm_sz))*n;
    local_b = ((my_rank + 1)/(0.0 + comm_sz))*n;

    if(my_rank == 0){
        local_a = 0;
    }else{
        local_a = (int)local_a;
    }
    if(my_rank == comm_sz - 1){
        local_b = n;
    }else{
        local_b = (int)local_b;
    }
    int numPairs = local_b-local_a;

    local_a = (local_a / n) * (b - a) + a;
    local_b = (local_b / n) * (b - a) + a;

    // local_int = Trap(local_a,local_b,local_n,h);
    local_int = Trap(local_a,local_b,numPairs,h);
    if(DEBUG){
        printf("Proc %d of %d: Trap calculation from %.3f to %.3f, with %d segment pairs and %f baselen\n",my_rank + 1, comm_sz,local_a,local_b,numPairs,h);
    }

    if(my_rank != 0){
        MPI_Send(&local_int, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }else{
        total_int = local_int;
        for(source = 1; source < comm_sz; source++){
            MPI_Recv(&local_int, 1, MPI_DOUBLE, source, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            total_int += local_int;
        }
    }

    if(my_rank == 0){
        printf ( "With n = %d trapezoids , our estimate \n" , n );
        printf ( "of the integral from %.3f to %.3f = %.3f\n" , a , b , total_int );
    }

    MPI_Finalize();
    return 0;

}
