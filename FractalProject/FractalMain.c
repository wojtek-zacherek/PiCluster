#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>
#include <math.h>
#include "UsefulFunction.h"


#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#define DEBUG 1
// make && mpirun -np 4 --mca btl_vader_single_copy_mechanism none FractalMain.bin

struct complexDouble{
    double R;
    double I;
};
typedef struct rgb{
    char r;
    char g;
    char b;
} rgb;
typedef struct colorSelection {
    struct rgb *selection;
    uint numEntries;
} colorSelection;
uint maxValue = 0;

void generateFractal(int my_rank, int comm_sz, uint xRes, uint yRes, uint iter, double xMin, double xMax, double yMin, double yMax, double rMax);
uint Mandelbrot(double C_R, double C_I, uint maxInt, double rMax);
void complexPow(struct complexDouble *operand, struct complexDouble *resultant, double power);
void complexSum(struct complexDouble *operand1, struct complexDouble *operand2, struct complexDouble *resultant);
void scalarSum(struct complexDouble *operand1, double R, struct complexDouble *resultant);
void complexADiffB(struct complexDouble *operand1, struct complexDouble *operand2, struct complexDouble *resultant);
void complexBDiffA(struct complexDouble *operand1, struct complexDouble *operand2, struct complexDouble *resultant);
void scalarADiffB(struct complexDouble *operand1, double R, struct complexDouble *resultant);
void scalarBDiffA(struct complexDouble *operand1, double R, struct complexDouble *resultant);
void complexMagnitude(struct complexDouble *operand, double *resultant);
void complexProduct(struct complexDouble *operand1, struct complexDouble *operand2, struct complexDouble *resultant);
void scalarProduct(struct complexDouble *operand1, double R, struct complexDouble *resultant);
void saveFile(uint xResolution, uint yResolution, uint thresh, uint **matrix, uint iter);
void makeColourfull(char **target, uint **matrix, uint xResolution, uint yResolution, uint iter);

// void checkLogic(){
//     struct complexDouble someValue = {1,2};
//     struct complexDouble result;

//     printf("%.2f + i(%.2f)\n",someValue.R,someValue.I);
//     complexPow(&someValue, &result, 2);
//     printf("%.2f + i(%.2f)\n",result.R,result.I);
// }

int main(int argc, char *argv[]){
    int my_rank, comm_sz;
    
    MPI_Init(NULL,NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    srand(time(NULL)+ my_rank);

    // checkLogic();
    MPI_Barrier(MPI_COMM_WORLD);

    double xMin = -2;
    double xMax = 1;
    double yMin = -1.2;
    double yMax = 1.2;
    // double xMin = -1.2;
    // double xMax = 1.2;
    // double yMin = -1.2;
    // double yMax = 1.2;

    // double xMin = -1.78;
    // double xMax = 1.78;
    // double yMin = -1;
    // double yMax = 1;

    uint xRes = 1920*1;    //889
    uint yRes = xRes / ((xMax - xMin) / (yMax - yMin));    //500
    // yRes = 1080;
    uint iter = 127;    // 127
    double rMax = 2;

    

    generateFractal(my_rank, comm_sz, xRes, yRes, iter, xMin, xMax, yMin, yMax, rMax);

    MPI_Finalize();
    return 0;
}

void generateFractal(int my_rank, int comm_sz, uint xRes, uint yRes, uint iter, double xMin, double xMax, double yMin, double yMax, double rMax){
    int i, j;
    uint resOffset = 1;
    uint **matrix;
    uint **submatrix;

    if(my_rank == 0){
        matrix = malloc((xRes + resOffset) * sizeof(uint*));
        for(i = 0; i < (xRes + resOffset); i++){
            matrix[i] = malloc((yRes + resOffset) * sizeof(uint));
        }
    }

    double xSlope = (double)(xMax - xMin) / xRes;
    double ySlope = (double)(yMax - yMin) / yRes;

    uint xDivs = comm_sz;
    uint yDivs = 1; 
    int xNumList[comm_sz];
    int yNumList[comm_sz];

    double xDivStep = ((double)xRes + resOffset) / xDivs;
    double yDivStep = ((double)yRes + resOffset) / yDivs;
    uint xStart, xEnd, yStart, yEnd;
    for(i = 0; i < comm_sz; i++){
        xStart = (i) * xDivStep;
        xEnd =  (i) * xDivStep + xDivStep;
        yStart = (yDivs - (0 + 1)) * yDivStep;
        yEnd = (yDivs - (0 + 1)) * yDivStep + yDivStep;
        xNumList[i] = xEnd - xStart;
        yNumList[i] = yEnd - yStart;
        
    }
    if(my_rank == 0){
        printVectorInt(xNumList,comm_sz);
    }
    xStart = ((my_rank + xDivs)- xDivs) * xDivStep;
    xEnd =  ((my_rank + xDivs) - xDivs) * xDivStep + xDivStep;
    yStart = (yDivs - (0 + 1)) * yDivStep;
    yEnd = (yDivs - (0 + 1)) * yDivStep + yDivStep;

    printf("Rank %d, xDivs %d, xStep %.2f, yDivs %d, yStep %.2f || xStart %d, xEnd %d, yStart %d, yEnd %d\n", my_rank, xDivs, xDivStep, yDivs, yDivStep, xStart, xEnd, yStart, yEnd);
    submatrix = malloc((xEnd - xStart) * sizeof(uint*));
    for(i = 0; i < ((xEnd - xStart)); i++){
        submatrix[i] = malloc((yEnd - yStart) * sizeof(uint));
        
    }
    
    uint val = 0;
    double tempX, tempY;
    // maxValue = 0;
    for(uint i = 0; i < (xEnd - xStart); i++){
        tempX = xSlope*((double)i - 0) + my_rank*(xMax - xMin) / comm_sz + xMin;
        for(uint j = 0; j < (yEnd - yStart); j++){
            tempY = ySlope*(j - 0) + yMin;
            // if(DEBUG == 1){
            //     printf("%d %d : %f %f\n",i,j,tempX,tempY);
            // }
            // uint val = Mandelbrot(tempX, tempY, thresh, iter);
            // if(my_rank == 0)
            //     printf("X,Y = %.2f,%.2f\n",tempX,tempY);
            val = Mandelbrot(tempX, tempY, iter, rMax);
            submatrix[i][j] = val;
            // sem_wait(&mutex);
            if(maxValue < submatrix[i][j]){
                maxValue = submatrix[i][j];
            }
            // sem_post(&mutex);
        }
    }
    printf("End %d\n",my_rank);
    
    if(my_rank == 0){
        // printSubMatrixFlip(submatrix, (yEnd - yStart), (xEnd - xStart));
    }

    
    
    if(my_rank == 0){
        for(i = 0; i < xNumList[0]; i++){
            copyArray(submatrix[i], matrix[i], yNumList[0]);
        }
        int offset = xNumList[0];
        for(i = 1; i < comm_sz; i++){
            for(j = 0; j < xNumList[i]; j++){
                MPI_Recv(matrix[offset], yNumList[i], MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                offset += 1;
            }
        }
    }else{
        for(j = 0; j < xNumList[my_rank]; j++){
            MPI_Send(submatrix[j], yNumList[my_rank], MPI_INT, 0, 0, MPI_COMM_WORLD);
        }
    }

    if(my_rank==0){
        // printSubMatrixFlip(matrix,yRes+1,xRes+1);
        saveFile(xRes, yRes, rMax, matrix, iter);
    }
}



uint Mandelbrot(double C_R, double C_I, uint maxInt, double rMax){
    
    // Z_(n+1) = (Z_n)^2 + C
    // Z is a complex point in the complex plain
    // C = Z_0
    uint i, maxIterations = maxInt, numIterations;
    double currentValue = 0;

    struct complexDouble Zn = {C_R, C_I};
    struct complexDouble C = {C_R, C_I};
    // struct complexDouble C = {0, -0.75};
    // struct complexDouble C = {0.28, 0.008};
    // struct complexDouble C = {0.3,-.01};
    // struct complexDouble C = {-0.835,-0.2321};
     

    struct complexDouble step1, step2, step3, step4, step5, step6, step7;
    
    
    for(i = 0; i < maxIterations; i++){
        complexMagnitude(&Zn, &currentValue);
        if(currentValue > rMax){
            break;
        }
        numIterations++;

        // MAIN COMPUTE BEGIN
        // complexPow(&Zn, &step1, 2);
        // complexSum(&step1, &C, &step2);
        // Zn.R = step2.R;
        // Zn.I = step2.I;

        complexPow(&Zn, &step1, 2);
        complexSum(&step1, &C, &step2);
        // scalarADiffB(&step2,cos(C.R),&step3);
        // scalarADiffB(&step2,cos(C.R)*sin(C.I),&step3);
        scalarProduct(&step2,cos(rand()),&step3);       // Basically a blurring function
        Zn.R = step3.R;
        Zn.I = step3.I;

        // MAIN COMPUTE END

        
        
    }
    numIterations--;
    return i;
}

void complexPow(struct complexDouble *operand, struct complexDouble *resultant, double power){
    double r = sqrt(pow(operand->R, 2) + pow(operand->I, 2));
    double theta;
    if(operand->R == 0){
        // theta = 2*acos(0.0);
        theta = 3.1415/2;
        // theta = 0;
    }else{
        theta = atan(operand->I / operand->R);
    }
    // theta = atan(operand->I / operand->R);
    double a = pow(r,power)*cos(power*theta);
    double b = pow(r,power)*sin(power*theta);
    resultant->R = a;
    resultant->I = b;
}

void complexSum(struct complexDouble *operand1, struct complexDouble *operand2, struct complexDouble *resultant){
    double a = operand1->R + operand2->R;
    double b = operand1->I + operand2->I;
    resultant->R = a;
    resultant->I = b;
}
void scalarSum(struct complexDouble *operand1, double R, struct complexDouble *resultant){
    double a = operand1->R + R;
    double b = operand1->I + 0;
    resultant->R = a;
    resultant->I = b;
}
void complexADiffB(struct complexDouble *operand1, struct complexDouble *operand2, struct complexDouble *resultant){
    double a = operand1->R - operand2->R;
    double b = operand1->I - operand2->I;
    resultant->R = a;
    resultant->I = b;
}
void complexBDiffA(struct complexDouble *operand1, struct complexDouble *operand2, struct complexDouble *resultant){
    double a = operand2->R - operand1->R;
    double b = operand2->I - operand1->I;
    resultant->R = a;
    resultant->I = b;
}
void scalarADiffB(struct complexDouble *operand1, double R, struct complexDouble *resultant){
    double a = operand1->R - R;
    double b = operand1->I - 0;
    resultant->R = a;
    resultant->I = b;
}
void scalarBDiffA(struct complexDouble *operand1, double R, struct complexDouble *resultant){
    double a = R - operand1->R;
    double b = 0 - operand1->I;
    resultant->R = a;
    resultant->I = b;
}


void complexMagnitude(struct complexDouble *operand, double *resultant){
    *resultant = sqrt(pow(operand->R, 2) + pow(operand->I, 2));
}

void scalarProduct(struct complexDouble *operand1, double R, struct complexDouble *resultant){
    double a = operand1->R * R ;
    double b = operand1->I * R + operand1->R * 0;
    resultant->R = a;
    resultant->I = b;
}

void complexProduct(struct complexDouble *operand1, struct complexDouble *operand2, struct complexDouble *resultant){
    double a = operand1->R * operand2->R - operand1->I * operand2->I;
    double b = operand1->I * operand2->R + operand1->R * operand2->I;
    resultant->R = a;
    resultant->I = b;
}

void saveFile(uint xResolution, uint yResolution, uint thresh, uint **matrix, uint iter){
    FILE* pgmimg;
    size_t fileTime = time(NULL);
    size_t fileSize = snprintf(NULL, 0, "%u_%u_%u_%lu.pgm", xResolution, yResolution, thresh, fileTime) + 1;
    size_t fileSize2 = snprintf(NULL, 0, "%u_%u_%u_%lu.jpeg", xResolution, yResolution, thresh, fileTime) + 1;
    size_t fileSize3 = snprintf(NULL, 0, "%u_%u_%u_%lu.png", xResolution, yResolution, thresh, fileTime) + 1;
    char* filename = malloc(fileSize);
    char* filename2 = malloc(fileSize2);
    char* filename3 = malloc(fileSize3);
    snprintf(filename, fileSize, "%u_%u_%u_%lu.pgm", xResolution, yResolution, thresh, fileTime);
    snprintf(filename2, fileSize2, "%u_%u_%u_%lu.jpeg", xResolution, yResolution, thresh, fileTime);
    snprintf(filename3, fileSize3, "%u_%u_%u_%lu.png", xResolution, yResolution, thresh, fileTime);

    char *data;
    maxValue = iter+1;
    makeColourfull(&data, matrix, xResolution, yResolution, iter);
    printf("Value: %d %d %d\n",data[0],data[xResolution*yResolution/2 + xResolution/2],data[xResolution*yResolution]);

    stbi_write_jpg(filename2,xResolution,yResolution,3,data,90);

}

void makeColourfull(char **target, uint **matrix, uint xResolution, uint yResolution, uint iter){
    int numColors = iter;
    char maxRGBValue = 255;
    struct rgb colors[numColors];
    struct rgb redwhite[512];
    struct rgb blue[256];
    struct rgb bluewhite[256];
    struct rgb blackredorangewhite[256+128+128];
    colorSelection colorSelected;

    colorSelected.selection = redwhite; colorSelected.numEntries = 511;
    colorSelected.selection = colors; colorSelected.numEntries = iter;
    colorSelected.selection = bluewhite; colorSelected.numEntries = 255;
    // colorSelected.selection = blue; colorSelected.numEntries = 255;
    colorSelected.selection = blackredorangewhite; colorSelected.numEntries = 511;

    for(int i = 0; i < 256; i++){
        blue[i].r = 0;
        blue[i].g = 0;
        blue[i].b = i;
    }

    for(int i = 0; i < 256; i++){
        blackredorangewhite[i].r = i;
        blackredorangewhite[i].g = 0;
        blackredorangewhite[i].b = 0;
    }
    for(int i = 0; i < 128; i++){
        blackredorangewhite[i + 256].r = 255;
        blackredorangewhite[i + 256].g = i;
        blackredorangewhite[i + 256].b = 0;
    }
    for(int i = 0; i < 128; i++){
        blackredorangewhite[i + 256 + 128].r = 255;
        blackredorangewhite[i + 256 + 128].g = 128 + i;
        blackredorangewhite[i + 256 + 128].b = 2*i;
    }

    for(int i = 0; i < 128; i++){
        bluewhite[i].r = 0;
        bluewhite[i].g = 0;
        bluewhite[i].b = 2*i;
    }
    for(int i = 0; i < 128; i++){
        bluewhite[i + 128].r = 2*i;
        bluewhite[i + 128].g = 2*i;
        bluewhite[i + 128].b = 255;
    }

    for(int i = 0; i < numColors; i++){
        colors[i].r = i;
        colors[i].g = 0;
        colors[i].b = 0;
    }


    for(int i = 0; i < 256; i++){
        redwhite[i].r = i;
        redwhite[i].g = 0;
        redwhite[i].b = 0;
    }
    for(int i = 0; i < 256; i++){
        redwhite[i + 256].r = 255;
        redwhite[i + 256].g = i;
        redwhite[i + 256].b = i;
    }

    // ratio = colorSelected.numEntries/maxValue
    printf("max = %d\n",maxValue);
    
    *target = malloc(xResolution*yResolution*sizeof(char)*3);
    char *data = *target;
    uint index = 0;
    uint index2 = 0;
    for(uint i = 0; i < xResolution; i++){
        for(uint j = 0; j < yResolution; j++){
            index = (int)(((long)matrix[i][j] * colorSelected.numEntries) / maxValue);
            index2 = 3*(i + j*xResolution);
            
            data[index2 + 0] = colorSelected.selection[index].r;
            data[index2 + 1] = colorSelected.selection[index].g;
            data[index2 + 2] = colorSelected.selection[index].b;
        }
    }
}