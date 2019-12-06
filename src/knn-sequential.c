
// C program to find groups of unknown 
// Points using K nearest neighbour algorithm. 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include "../inc/knn.h"
#include <cblas.h>   //sudo apt-get install libatlas-base-dev 



//Populate array/dataset with random double values
void populateArray(double *array, int n, int d){
	  for (int i = 0; i < n; i++) {
			for (int j = 0; j < d; j++) {
				*(array + i*d + j) = (double)rand()/RAND_MAX;
			}
	  }
}
//Helper to visualize results
void printArray(double *array, int n, int d){
    for (int i = 0; i < n; ++i)
    {
       // cout << i << ": ";
        for (int j = 0; j < d; ++j){
            printf("%f ", *(array + i*d + j));
		}
		printf("\n");
	}
}
//custom swap for doubles
void swap_d(double *p, double *q){
    double tmp;

    tmp = *p;
    *p = *q;
    *q = tmp;
}
//custom swap for doubles
void swap_i(int *p, int *q){
    double tmp;
    tmp = *p;
    *p = *q;
    *q = tmp;
}
// Standard partition process of QuickSort().
// It considers the last element as pivot
// and moves all smaller element to left of
// it and greater elements to right
int partitionWithIndex(double *arr, int *idx,int l, int r){
	double x = arr[r];
	int i = l;
	for (int j = l; j <= r - 1; j++) {
		if (*(arr +j) <= x) {
			
			swap_d(arr+i, arr +j); 
			swap_i(idx+i, idx +j); 
			i++;
		}
	}
	swap_d(arr + i, arr + r); 
	swap_i(idx+i, idx +r);
	return i;
}
int partition(double *arr, int l, int r){
	double x = arr[r];
	int i = l;
	for (int j = l; j <= r - 1; j++) {
		if (*(arr +j) <= x) {
			swap_d(arr+i, arr +j); 
			
			i++;
		}
	}
	swap_d(arr + i, arr + r); 
	
	return i;
}
// This function returns k'th smallest with index
// element in arr[l..r] using QuickSort
// based method.  ASSUMPTION: ALL ELEMENTS
// IN ARR[] ARE DISTINCT
double kthSmallestWithIndex(double *arr, int *idx,int l, int r, int k){
    // If k is smaller than number of
    // elements in array
    if (k > 0 && k <= r - l + 1) {

        // Partition the array around last
        // element and get position of pivot
        // element in sorted array
        int index = partitionWithIndex(arr, idx, l, r);

        // If position is same as k
        if (index - l == k - 1)
            return *(arr + index);

        // If position is more, recur
        // for left subarray
        if (index - l > k - 1)
            return kthSmallestWithIndex(arr, idx,l, index - 1, k);

        // Else recur for right subarray
        return kthSmallestWithIndex(arr, idx,index + 1, r,
            k - index + l - 1);
    }
    // If k is more than number of
    // elements in array
    return INT_MAX;
} 
double kthSmallest(double *arr, int l, int r, int k){
    // If k is smaller than number of
    // elements in array
    if (k > 0 && k <= r - l + 1) {

        // Partition the array around last
        // element and get position of pivot
        // element in sorted array
        int index = partition(arr, l, r);

        // If position is same as k
        if (index - l == k - 1)
            return *(arr + index);

        // If position is more, recur
        // for left subarray
        if (index - l > k - 1)
            return kthSmallest(arr, l, index - 1, k);

        // Else recur for right subarray
        return kthSmallest(arr,index + 1, r,
            k - index + l - 1);
    }
    // If k is more than number of
    // elements in array
    return INT_MAX;
} 

//function to calc euclidean distances with cblas lib based on matlab given type: sqrt(sum(X.^2,2) -2* X*Y.' + sum(Y.^2,2).') <= (tautotita)
double* calcDistanceBlas(double * X, double * Y, int n, int m, int d, int k){
		double count = 0;
		double alpha= -2;
		double beta = 0;
		double * distance = (double *)malloc(m * n *sizeof(double));
		double * Xnrm = (double *)malloc(n *sizeof(double)); // L2 norm of X
		double * Ynrm = (double *)malloc(m *sizeof(double)); // L2 norm of Y
		double * bothMatrix = (double *)malloc(m * n *sizeof(double));
		
		// dgemm routine, calculates the product of double precision matrices (-2*X*Y')
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, n, m, d, alpha, X, d, Y, d, beta, bothMatrix, m);

		
		//	double cblas_dnrm2(const int __N, const double *__X, const int __incX);
		
		//	N		Length of vector X.
		//	X		Vector X.
		//	incX 	Stride within X. For example, if incX is 7, every 7th element is used.
		
		for(int i = 0; i < n; i++){
			count = cblas_dnrm2(d, X+i*d, 1); //cblas_dnrm2:	Computes the L2 norm (Euclidian length) of a vector (double precision).
			*(Xnrm +i) = count * count;
		}
		for(int i = 0; i < m; i++){
			count = cblas_dnrm2(d, Y+i*d, 1); //cblas_dnrm2:	Computes the L2 norm (Euclidian length) of a vector (double precision).
			*(Ynrm +i) = count * count; 
		}

		for (int i=0; i<n; i++){
			for(int j=0; j<m; j++){
				*(bothMatrix + i*m+j) += *(Xnrm +i) + *(Ynrm +j);
			}
		}

		for(int i = 0; i < n*m; i++){
			*(distance +i) = sqrt(*(bothMatrix +i));
		}

		// Allocate
		free(Xnrm);
		free(Ynrm);
		free(bothMatrix);

		return distance;
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

	knnresult knnres;
	double * distance = (double *)malloc(m * n  *sizeof(double));
	double * distancetmp = (double *)malloc(m * n  *sizeof(double));
	int *indexes = (int*)malloc(m * n  *sizeof(int));	
	double alpha= 1;
	double beta = 0;
	
	knnres.m=m;
	knnres.k=k;
	
	
	for(int i=0;i<m;i++){
		for(int j=0;j<n;j++){
			indexes[i*n+j]=j;
		}
	}
	
	double * unitMatrix = (double *)malloc(n * n  *sizeof(double)); //size n is the n × n square matrix
	//populate unit matrix
	for(int i = 0; i < n;i++) {
		for(int j = 0; j < n;j++) {
			if(i == j) {
				*(unitMatrix + i*n + i) = 1; //ones on the main diagonal
			} 
		}
	}
	
	double* tempDistance = calcDistanceBlas(X, Y, n, m, d, k);
	cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, m, n, n, alpha, tempDistance, m, unitMatrix, n, beta, distance, n);

	// Allocate
	free(unitMatrix);
	free(tempDistance);
	
	knnres.ndist=(double *)malloc(m*k * sizeof(double));
	knnres.nidx=(int *)malloc(m*k * sizeof(int));
	//Calculates the minimum distance of each point of y from X dataset
	for(int i=0;i<m;i++){
		printf("\n Gia i= %d : \n",i);
		for(int j=0;j<k;j++){
			printf("MPAINEI:	");
			knnres.ndist[i*k+j]=kthSmallestWithIndex(&distance[i*n], &indexes[i*n],0, n-1, j+1);
			
			knnres.nidx[i*k+j]=indexes[i*n+j];
			
			printf("O geitonas %d exei apostasi: %lf \n",knnres.nidx[i*k+j],knnres.ndist[i*k+j]);
			
		}
	}
	
	
	free(distance);
	free(indexes);
	return knnres;	
}





