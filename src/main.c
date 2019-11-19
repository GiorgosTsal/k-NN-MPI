#include <stdio.h>
#include <stdlib.h>
#include "../inc/knn.h"




const int n = 3, d = 4, m = 3, k=2;


// Driver code 
int main() 
{ 
	
	//double *X = (double *)malloc(n * d * sizeof(double));
//	double *Y = (double *)malloc(m * d * sizeof(double));
	//populateArray(X, n, d);
	double X[3][4] = {{1000,1,1,1}, {3,3,3,3}, {6,6,6,6}};
	printf("Array X is: ");
	printf("\n");
	printArray(X, n, d);
	printf("\n");
	printf("Array Y is: ");
	printf("\n");
//	populateArray(Y, m, d);
	double Y[3][4] = {{0,1,2,3}, {8,8,4,8}, {11,9,9,9}};
	printArray(Y, m, d);

	int *idxX = (int *)malloc(n *sizeof(int));
	int *idxY = (int *)malloc(m *sizeof(int));
	idxX = indexArray(n);
	idxY = indexArray(m);
	
	printf("\n");


	printf("begin \n");
    knnresult result;
	result = kNN(X, Y, n, m, d, k);
	
	printf("end \n");
   

    return 0; 
} 


