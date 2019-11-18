
// C program to find groups of unknown 
// Points using K nearest neighbour algorithm. 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../inc/knn.h"



//Populate array/dataset with random double values
void populateArray(double *array, int n, int d){
	  for (int i = 0; i < n; i++) {
			for (int j = 0; j < d; j++) {
				*(array + i*d + j) = (double)rand()/RAND_MAX;
			}
	  }
}
////Helper to visualize results
void printArray(double *array, int n, int d)
{
    for (int i = 0; i < n; ++i)
    {
       // cout << i << ": ";
        for (int j = 0; j < d; ++j)
            printf("%f\n", *(array + i*d + j));
			putchar('\n');
			fflush(stdout);
    }
}
//Print Indices of an array
void printIndicesArray(int *array, int size)
{
    for (int i = 0; i < size; ++i)
    {
			printf("%d\n", *(array + i));
    }
}
//Indexes an array
int * indexArray(int size){
	int * idxOfarray = (int *)malloc(size *sizeof(int));
	for (int i = 0; i < size; i++) {
			*(idxOfarray + i) = i;
	}
	return idxOfarray;
}
 
double dist(){
	
}

//! Compute k nearest neighbors of each point in X [n-by-d]
/*!
	\param	X		Corpus data points			[n-by-d]
	\param	Y		Query data points			[m-by-d]
	\param	n		Number of corpus points		[scalar]
	\param	m		Number of query points		[scalar]
	\param	d		Number of dimensions		[scalar]
	\param	k		Number of neighbors			[scalar]

	\return The kNN result
*/
knnresult kNN(double* X, double* Y, int n, int m, int d, int k){
	//// Fill distances of all points from p 
	 printf("ds");
	double * distance = (double *) malloc((n*m)*sizeof(double));
	 
	double *xRow = (double *)malloc(n * sizeof(double));
	double *yRow = (double *)malloc(m * sizeof(double));
	
	for(int i=0; i<n; i++){
		for(int j=0; j<d; j++){
			xRow[i] += (*(X+i*d+j)) * (*(X+i*d+j));
		}
	}
	
	for(int i=0; i<m; i++){
		for(int j=0; j<d; j++){
			yRow[i] += (*(Y+i*d+j)) * (*(Y+i*d+j));
		}
	}
		
	for(int i=0; i<n; i++){
		for(int j=0; j<m; j++){
			*(distance + i*m + j) += xRow[i] + yRow[j];
			*(distance + i*m + j) = sqrt( *(distance + i*m + j) );
		}
	}

	printArray(distance, n, m);
   
}




