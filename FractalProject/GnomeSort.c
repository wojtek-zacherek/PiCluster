/**
 * @file GnomeSort.c
 * @author Wojciech Zacherek (zacherw@rose-hulman.com)
 * @brief 
 * @version 0.1
 * @date 2022-05-21
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>
#include <math.h>
#include "UsefulFunction.h"
#include <sys/time.h>
#include <string.h>

#define DEBUG 1
// make && mpirun -np 4 --mca btl_vader_single_copy_mechanism none FractalMain.bin

struct timePair
{
    struct timeval start;
    struct timeval end;
    char name[30];
};
struct timePair program, dataInitMain, individualSort, combineSort;
#define sizeTimeArray 4
struct timePair timeArray[sizeTimeArray];
double calculateTimes(struct timePair x)
{
    double time_taken = x.end.tv_sec + x.end.tv_usec / 1e6 - x.start.tv_sec - x.start.tv_usec / 1e6;
    return time_taken;
}
// gettimeofday(&timeArray[1].start,NULL); strcpy(timeArray[1].name,"dataInitKernel");
// gettimeofday(&timeArray[1].end,NULL);

int cmpfunc(const void *a, const void *b)
{
    return (*(int *)a - *(int *)b);
}

void swap(int *a, int *b)
{
    int t = *a;
    *a = *b;
    *b = t;
}

void gnomeSort(int arr[], int n)
{
    int index = 0;

    while (index < n)
    {
        if (index == 0)
            index++;
        if (arr[index] >= arr[index - 1])
            index++;
        else
        {
            swap(&arr[index], &arr[index - 1]);
            index--;
        }
    }
    return;
}

void popArr(int a[], int N)
{
    int i;
    for (i = 0; i < N; i++)
    {
        a[i] = (int)rand() % 1000;
    }
}

void divideWork(int my_rank, int comm_sz, int N, int workDivision[])
{
    int nExtra = N % comm_sz;
    int rate = N / comm_sz;
    int i;

    for (i = 0; i < comm_sz; i++)
    {
        workDivision[i] = rate;
        if (nExtra > 0)
        {
            workDivision[i] += 1;
            nExtra--;
        }
    }
}

void printArray(int arr[], int n)
{
    printf("Sorted sequence after Gnome sort:\n");
    for (int i = 0; i < n; i++)
        printf("%d ", arr[i]);
    printf("\n");
}

void verifySort(int arr[], int N)
{
    int last = arr[0];
    int i;
    for (i = 0; i < N; i++)
    {
        if (arr[i] < last)
        {
            printf("Error at %d: %d !< %d\n", i, last, arr[i]);
        }
        else
        {
            last = arr[i];
        }
    }
    if (i == N)
    {
        printf("Sorting verified!\n");
    }
}

int main(int argc, char *argv[])
{
    int my_rank, comm_sz;

    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    srand(time(NULL) + my_rank);
    int i;

    if (DEBUG == 1 && my_rank == 0)
    {
        gettimeofday(&timeArray[0].start, NULL);
        strcpy(timeArray[0].name, "Program");
        gettimeofday(&timeArray[1].start, NULL);
        strcpy(timeArray[1].name, "Init Data");
    }
    int N = 10;

    //////////////////////////////////////////////////
    int a[N];
    int workDivision[comm_sz];

    divideWork(my_rank, comm_sz, N, workDivision);
    int aLocal[workDivision[my_rank]];
    int *currentWorkingArr;
    currentWorkingArr = malloc(workDivision[my_rank] * sizeof(int));
    if (my_rank == 0)
    {
        popArr(a, N);
        // printVectorInt(a, N);
        copyArray(a, currentWorkingArr, workDivision[my_rank]);

        if (DEBUG == 1 && my_rank == 0)
        {
            gettimeofday(&timeArray[1].end, NULL);
        }

        int offset = workDivision[my_rank];
        for (i = 1; i < comm_sz; i++)
        {
            MPI_Send(a + offset, workDivision[i], MPI_INT, i, 0, MPI_COMM_WORLD);
            offset += workDivision[i];
        }
    }
    else
    {
        MPI_Recv(currentWorkingArr, workDivision[my_rank], MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    if (DEBUG == 1 && my_rank == 0)
    {
        gettimeofday(&timeArray[2].start, NULL);
        strcpy(timeArray[2].name, "Gnome Sort");
    }
    gnomeSort(currentWorkingArr, workDivision[my_rank]);
    if (DEBUG == 1 && my_rank == 0)
    {
        gettimeofday(&timeArray[2].end, NULL);
    }


    if (DEBUG == 1 && my_rank == 0)
    {
        gettimeofday(&timeArray[3].start, NULL);
        strcpy(timeArray[3].name, "Combine Sort");
    }
    int index = 2;
    while (index/2 <= comm_sz)
    {
        if (my_rank % index == 0)
        {

            // Check if other half exists
            int targetIndex = my_rank + (index / 2);
            if (targetIndex < comm_sz)
            {
                int incomingSize = 0;
                int localSize = workDivision[my_rank];

                // Get incoming size, make array to fit both incoming and local arrays, and recv incoming data
                MPI_Recv(&incomingSize, 1, MPI_INT, targetIndex, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                int totalSize = localSize + incomingSize;
                int *arrIncoming = malloc(totalSize * sizeof(int));
                // printf("Proc %d: localSize = %d, incoming = %d, total = %d\n",my_rank,localSize,incomingSize,totalSize);

                // printVectorInt(currentWorkingArr,localSize);
                copyArray(currentWorkingArr, arrIncoming, localSize);
                MPI_Recv(arrIncoming + localSize, incomingSize, MPI_INT, targetIndex, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                // Combine sort
                qsort(arrIncoming, totalSize, sizeof(int), cmpfunc);

                // Update local values
                workDivision[my_rank] = totalSize;
                free(currentWorkingArr);
                currentWorkingArr = arrIncoming;
            }
            else
            {
                // Otherwise do nothing.
            }
        }
        else
        {

            // Get proc to which data will be sent to.
            int targetIndex = my_rank - (index / 2);
            int outgoingSize = workDivision[my_rank];
            MPI_Send(&outgoingSize, 1, MPI_INT, targetIndex, 0, MPI_COMM_WORLD);
            MPI_Send(currentWorkingArr, outgoingSize, MPI_INT, targetIndex, 0, MPI_COMM_WORLD);
            // Once data is sent, this proc is done.
            break;
        }
        index *= 2;
    }
    if (DEBUG == 1 && my_rank == 0)
    {
        gettimeofday(&timeArray[3].end, NULL);
    }

    if (my_rank == 0)
    {
        // printArray(currentWorkingArr, N);
        verifySort(currentWorkingArr, N);
    }

    if (DEBUG == 1 && my_rank == 0)
    {
        gettimeofday(&timeArray[0].end, NULL);
    }

    if (DEBUG && my_rank == 0)
    {
        for (i = 0; i < sizeTimeArray; i++)
        {
            printf("Time for \"%s\" = %.3f s\n", timeArray[i].name, calculateTimes(timeArray[i]));
        }
    }
    free(currentWorkingArr);
    MPI_Finalize();
    return 0;
}
