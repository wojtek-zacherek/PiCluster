#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <sys/time.h>

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
    int my_rank, comm_sz, n = 8096;
    int nn[] = {1024,2048,4096,8096,16384,32768,65536,131072,262144,524288,1048576,2097152,4194304,8388608};
    int nn_results1[14], nn_results2[14];
    double a = 0.0, b = 3.0, h, local_a, local_b;
    double local_int, total_int;
    double local_n;
    int source;
    struct timeval startProgram, startKernel, endKernel, endCollection, endProgram;
    int repeats = 100;
    int i, j;

    gettimeofday(&startProgram, NULL);

    if(argc > 1){
        n = atoi(argv[1]);
        // printf("n = %d\n",n);
        if(argc > 2){
            repeats = atoi(argv[2]);
            // printf("type = %d\n",type);
        }
    }


    MPI_Init(NULL,NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);

    for(i = 0; i < repeats; i++){
    
        for(j = 0; j < 14; j++){
            h = (b-a)/n;
            local_n = ((double)n)/comm_sz;

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

            
            MPI_Barrier(MPI_COMM_WORLD);
            gettimeofday(&startKernel, NULL);
            local_int = Trap(local_a,local_b,numPairs,h);
            MPI_Barrier(MPI_COMM_WORLD);
            gettimeofday(&endKernel, NULL);

            if(DEBUG){
                printf("Proc %d of %d: Trap calculation from %.3f to %.3f, with %d segment pairs and %f baselen\n",my_rank + 1, comm_sz,local_a,local_b,numPairs,h);
            }


            MPI_Reduce(&local_int , &total_int , 1, MPI_DOUBLE , MPI_SUM , 0 , MPI_COMM_WORLD ) ;
            MPI_Barrier(MPI_COMM_WORLD);
            gettimeofday(&endProgram, NULL);

            unsigned long kernalTime = (endKernel.tv_sec - startKernel.tv_sec) * 1000000 + endKernel.tv_usec - startKernel.tv_usec;
            unsigned long transferTime = (endProgram.tv_sec - endKernel.tv_sec) * 1000000 + endProgram.tv_usec - endKernel.tv_usec;
            nn_results1[j] = kernalTime;
            nn_results2[j] = transferTime;
        }
        
        
        

        if(my_rank == 0){

            // Save Results to CSV file
            FILE *file;
            char name[100];
            snprintf(name, sizeof(name), "Data_TTRT_%d_%d.csv", comm_sz, n);

            if (file = fopen(name, "r")) {
                fclose(file); 
                file = fopen(name, "a");
            } else {
                file = fopen(name, "a");
                for(j = 0; j < 14; j++){
                    fprintf(file,"Kernel Time %d (us), Transfer Time %d (us), ", nn[j], nn[j]);
                }
                fprintf(file,"\n");
            }
            for(j = 0; j < 14; j++){
                fprintf(file,"%lu, %lu,",nn_results1[j], nn_results2[j]);
            }
            fprintf(file,"\n");
            fclose(file);
        }
    }


    MPI_Finalize();
    return 0;

}
