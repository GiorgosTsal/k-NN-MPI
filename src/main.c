#include <stdio.h>
#include <stdlib.h>
#include "../inc/knn.h"




const int n = 3, d = 4, m = 3, k=3;


// Driver code 
int main() 
{ 

	
	double *X = (double *)malloc(n * d * sizeof(double));
	double *Y = (double *)malloc(m * d * sizeof(double));
	populateArray(X, n, d);
	printf("Array X is: ");
	printf("\n");
	printArray(X, n, d);
	printf("\n");
	printf("Array Y is: ");
	printf("\n");
	populateArray(Y, m, d);
	printArray(Y, m, d);

	int *idxX = (int *)malloc(n *sizeof(int));
	int *idxY = (int *)malloc(m *sizeof(int));
	idxX = indexArray(n);
	idxY = indexArray(m);
	
	printf("\n");


	printf("begin \n");
	double  total_time;
	struct timeval start, end;
	gettimeofday (&start, NULL);
    knnresult result;
	result = kNN(X, Y, n, m, d, k);
	gettimeofday (&end, NULL);
	total_time = (double)((end.tv_usec - start.tv_usec)/1.0e6 + end.tv_sec - start.tv_sec);
	printf("Time taken by kNN function: : %f sec\n",total_time);
	

	printf("end \n");
   
    return 0; 
} 


