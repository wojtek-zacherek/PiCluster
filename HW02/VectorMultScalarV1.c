#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <time.h>

#define INPUT 0
#define FIXED_N 1

double dotProductScale(double* a, double* b, double* c, double alpha, int n){
    int i;
    for(i = 0; i < n; i++){
        c[i] = alpha*a[i]*b[i];
    }

}

int main(int argc, char *argv[]){
    int my_rank, comm_sz, n = 256;
    double *a = NULL, *b = NULL, *c = NULL;

    double alpha = -1;
    int numPerProc = -1;
    int i;
    srand(time(NULL));

    MPI_Init(NULL,NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);

    if(FIXED_N){
        n = comm_sz * 100;
    }
    numPerProc = n / comm_sz;

    if ( my_rank == 0) {
        a = malloc(n*sizeof(double));
        b = malloc(n*sizeof(double));
        c = malloc(n*sizeof(double));
        if(INPUT){
            printf("Enter the vector A\n");
            for (i = 0; i < n; i ++){
                scanf("%lf", &a[i]);
            }
            printf("Enter the vector B\n");
            for(i = 0; i < n; i++){
                scanf("%lf", &a[i]);
            }
            printf("Enter the scalar\n");
            scanf("%lf", &alpha);

        }else{
            for(i = 0; i < n; i++){
                a[i] = rand() % 100;
                b[i] = rand() % 50;
            }
            alpha = rand() % 15;

        }

        for(i = 1; i < comm_sz; i++){
            MPI_Send(a + i * numPerProc, numPerProc, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
            MPI_Send(b + i * numPerProc, numPerProc, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
            MPI_Send(&alpha, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
        }

    }else{
        a = malloc(numPerProc*sizeof(double));
        b = malloc(numPerProc*sizeof(double));
        c = malloc(numPerProc*sizeof(double));

        MPI_Recv(a, numPerProc, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(b, numPerProc, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&alpha, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    dotProductScale(a, b, c, alpha, numPerProc);

    if(my_rank == 0){
        for(i = 1; i < comm_sz; i++){
            MPI_Recv(c + i * numPerProc, numPerProc, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }else{
        MPI_Send(c, numPerProc, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    if(my_rank == 0){
        double* c_test = malloc(n*sizeof(double));
        dotProductScale(a, b, c_test, alpha, n);
        int error = 0;
        for(i = 0; i < n; i++){
            if(c[i] != c_test[i]){
                printf("\ti=%d, alpha = %lf, a=%lf, b=%lf, c=%lf, c_test=%lf\n",i,alpha,a[i],b[i],c[i],c_test[i]);
                error = 1;
            }
        }
        if(error == 0){
            printf("Successfully verified.\n");
            printf("Sample Comparrison:\n");
            for(i = 0; i < 10; i++){
                printf("\ti=%d, alpha = %lf, a=%lf, b=%lf, c=%lf, c_test=%lf\n",i,alpha,a[i],b[i],c[i],c_test[i]);
            }

        }
        free(c_test);
    }    

    free(a);
    free(b);
    free(c);
    return 0;
}

// MPI_Scatter ( a , local_n , MPI_DOUBLE , local_a , local_n , MPI_DOUBLE , 0, comm ) ;
// MPI_Scatter ( a , local_n , MPI_DOUBLE , local_a , local_n , MPI_DOUBLE , 0, comm ) ;
