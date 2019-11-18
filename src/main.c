#include <stdio.h>
#include <stdlib.h>
#include "../inc/knn.h"




const int n = 10, d = 2, m = 10, k=2;


// Driver code 
int main() 
{ 
	
	double *X = (double *)malloc(n * d * sizeof(double));
	double *Y = (double *)malloc(m * d * sizeof(double));
	populateArray(X, n, d);
	printArray(X, n, d);
	printf("\n");
	populateArray(Y, m, d);
	printArray(Y, m, d);

	int *idxX = (int *)malloc(n *sizeof(int));
	int *idxY = (int *)malloc(m *sizeof(int));
	idxX = indexArray(n);
	idxY = indexArray(m);
	printIndicesArray(idxX, n);
	printf("\n");
	printIndicesArray(idxY, m);

	printf("begin \n");
    knnresult result;
	result = kNN(X, Y, n, m, d, k);
	
	printf("end \n");
   

    return 0; 
} 


