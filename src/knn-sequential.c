
// C program to find groups of unknown 
// Points using K nearest neighbour algorithm. 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
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
void printArray(double *array, int n, int d){
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
void printIndicesArray(int *array, int size){
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
//custom swap 
void swap(double *p, double *q){
    double tmp;

    tmp = *p;
    *p = *q;
    *q = tmp;
}
// Standard partition process of QuickSort().
// It considers the last element as pivot
// and moves all smaller element to left of
// it and greater elements to right
int partition(double arr[], int l, int r){
    double x = arr[r];
    int i = l;
    for (int j = l; j <= r - 1; j++) {
        if (*(arr +j) <= x) {
            swap(arr+i, arr +j); 
            i++;
        }
    }
    swap(arr + i, arr + r); 
    return i;
}
// This function returns k'th smallest
// element in arr[l..r] using QuickSort
// based method.  ASSUMPTION: ALL ELEMENTS
// IN ARR[] ARE DISTINCT
double kthSmallest(double arr[], int l, int r, int k){
    // If k is smaller than number of
    // elements in array
    if (k > 0 && k <= r - l + 1) {

        // Partition the array around last
        // element and get position of pivot
        // element in sorted array
        int index = partition(arr, l, r);

        // If position is same as k
        if (index - l == k - 1)
            return arr[index];

        // If position is more, recur
        // for left subarray
        if (index - l > k - 1)
            return kthSmallest(arr, l, index - 1, k);

        // Else recur for right subarray
        return kthSmallest(arr, index + 1, r,
            k - index + l - 1);
    }
    // If k is more than number of
    // elements in array
    return INT_MAX;
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

	double * distance = (double *) malloc(m*sizeof(double));
	 
	double *xRow = (double *)malloc(n * sizeof(double));
	double *yRow = (double *)malloc(m * sizeof(double));
	
	for(int i=0; i<n; i++){
		for(int j=0; j<d; j++){
			*(xRow +i) += (*(X+i*d+j)) * (*(X+i*d+j)); //sum(X.^2,2) 
		}
	}
	printf("xRow: \n");
	printArray(xRow, n, 1);
	for(int i=0; i<m; i++){
		for(int j=0; j<d; j++){
			*(yRow +i) += (*(Y+i*d+j)) * (*(Y+i*d+j)); //sum(Y.^2,2)
		}
	}
	printf("yRow: \n");
	printArray(yRow, m, 1);
	for(int i=0; i<n; i++){
		for(int j=0; j<m; j++){
			*(distance + i*m + j) += sqrt(*(xRow +i) + *(yRow +i));
		}
	}
	printf("Distances: \n");
	printArray(distance, m, 1);

	double *kNeighbours = (double *)malloc(k * sizeof(double));
	for(int i=1; i<=k; i++){
		printf("oi %d neibhor ine: %f ",i , kthSmallest(distance, 0, m-1, i));
		*(kNeighbours + i) = kthSmallest(distance, 0, m-1, i);
		printf("\n");
		//na to valo na gyrnaei k ta indices tous-kamia domi gia to dist??
	}
	
   
}




