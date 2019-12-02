// C program to find groups of unknown 
// Points using K nearest neighbour algorithm and MPI lib for communication. 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include "../inc/knn.h"
#include <mpi.h>
#include <string.h>




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
	int *indexes = (int*)malloc(m * n  *sizeof(int));
	
	//Declaration
	//void cblas_dgemm(const enum CBLAS_ORDER __Order, const enum CBLAS_TRANSPOSE __TransA, const enum CBLAS_TRANSPOSE __TransB, const int __M, const int __N, const int __K, const double __alpha, const double *__A, const int __lda, const double *__B, const int __ldb, const double __beta, double *__C, const int __ldc);
	
	//This function multiplies A * B and multiplies the resulting matrix by alpha. It then multiplies matrix C by beta. It stores the sum of these two products in matrix C.

	//Thus, it calculates either

	//C←αAB + βC

	//or

	//C←αBA + βC

	//with optional use of transposed forms of A, B, or both.
	
	
	//Parameters

	//Order

		//Specifies row-major (C) or column-major (Fortran) data ordering.
	//TransA

		//Specifies whether to transpose matrix A.
	//TransB

		//Specifies whether to transpose matrix B.
	//M

		//Number of rows in matrices A and C.
	//N

		//Number of columns in matrices B and C.
	//K

		//Number of columns in matrix A; number of rows in matrix B.
	//alpha

		//Scaling factor for the product of matrices A and B.
	//A

		//Matrix A.
	//lda

		//The size of the first dimention of matrix A; if you are passing a matrix A[m][n], the value should be m.
	//B

		//Matrix B.
	//ldb

		//The size of the first dimention of matrix B; if you are passing a matrix B[m][n], the value should be m.
	//beta

		//Scaling factor for matrix C.
	//C

		//Matrix C.
	//ldc

		//The size of the first dimention of matrix C; if you are passing a matrix C[m][n], the value should be m.
		
	
	knnres.m=m;
	knnres.k=k;
	
	// Calculates euclidean distance.
	for(int i=0;i<m;i++){
		for(int j=0;j<n;j++){
			double dist=0;
			for(int s=0;s<d;s++){
				dist+=(X[j*d+s]-Y[i*d+s])*(X[j*d+s]-Y[i*d+s]);		
			}
			indexes[i*n+j]=j;
			distance[i*n+j]=sqrt(dist);
		}
	}
	
	knnres.ndist=(double *)malloc(m*k * sizeof(double));
	knnres.nidx=(int *)malloc(m*k * sizeof(int));
	//Calculates the minimum distance of each point of y from X dataset
	for(int i=0;i<m;i++){
		for(int j=0;j<k;j++){
			knnres.ndist[i*k+j]=kthSmallestWithIndex(&distance[i*n], &indexes[i*n],0, n-1, j+1);
			knnres.nidx[i*k+j]=indexes[i*n+j];
		}
	}
	
	free(distance);
	free(indexes);
	return knnres;	
}

//! Compute distributed all-kNN of points in X
/*!
	\param	X	Data points				[n-by-d]
	\param	n	Number of data points	[scalar]
	\param	d	Number of dimensions	[scalar]
	\param	k	Number of neighbors		[scalar]
  
	\return	The kNN result
*/
knnresult distrAllkNN(double * X, int n, int d, int k){

	// Find out rank, size
	int process_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &process_rank); //the current ID of the running process
	int pr;
	MPI_Comm_size(MPI_COMM_WORLD, &pr); //total number of running MPI processes
	//Set MPI communication tag and status.
    int tag = 0; //suggested 0
    MPI_Status status;
	
	// which process to send?
	int from , to;
	if(process_rank == 0)
		from = pr - 1;
	else
		from = process_rank - 1;
	if(process_rank == pr - 1)
		to = 0;
	else
		to = process_rank + 1;


	int branks = process_rank;


	// Generate ranks matrix
	int* process_ranks = calloc(n, sizeof(int));
	for(int i = 0; i < n; i++){
		// 0 => last n points, 1 => [0, n], 2 => [n+1, 2n] ktlp 
		if(process_rank != 0)
			*(process_ranks +i) = (process_rank - 1) * n + i;
		else
			*(process_ranks +i) = (pr-1) * n + i;
	}
	// run the kNN code to get the first result
	knnresult res = kNN(X, X, n, n, d, k);


	// map ranks in order to send only the corpus
	for(int i = 0; i < n * k; i++){
		*(res.nidx + i) = *(process_ranks + res.nidx[i]);
	}

	knnresult tempRes;

	double * Y = (double *)malloc(n * d  *sizeof(double));
	memcpy(Y, X, n * d * sizeof(double));
	double * tempY = (double *)malloc(n * d  *sizeof(double));

	int idxRes = 0;
	int idxNres = 0;
	double* tempResDis = (double *)malloc(n * k  *sizeof(double));
	int* tempResIDs = (int *)malloc(n * k  *sizeof(int));

	for(int i = 0; i < pr-1; i++){

		// if even, first send and then receive
		if(process_rank % 2 == 0){
			
			MPI_Send(Y, n*d, MPI_DOUBLE, to, tag, MPI_COMM_WORLD);

			MPI_Recv(Y, n*d, MPI_DOUBLE, from, tag, MPI_COMM_WORLD, &status);
			
		}else{//first recv and then send

			MPI_Recv(tempY, n*d, MPI_DOUBLE, from, tag, MPI_COMM_WORLD, &status);

			MPI_Send(Y, n*d, MPI_DOUBLE, to, tag, MPI_COMM_WORLD);
			// Write tempY the originl
			memcpy(Y, tempY, n * d * sizeof(double));
		}

		branks = branks - 1;
		if(branks < 0)
			branks = pr-1;
		for(int j = 0; j < n; j++){
			if(branks != 0)
				*(process_ranks + j) = (branks - 1) * n + j;
			else
				*(process_ranks + j) = (pr-1) * n + j;
		}

		tempRes = kNN(Y, X, n, n, d, k);

		// map ranks back 
		for(int i = 0; i < n * k; i++){
			tempRes.nidx[i] = *(process_ranks + tempRes.nidx[i] );
		}	

		memcpy(tempResDis, res.ndist, n * k * sizeof(double));
		memcpy(tempResIDs, res.nidx, n * k * sizeof(int));

		// loop through query points for merging
		for(int s = 0; s < n; s++){
			idxRes = 0;
			idxNres = 0;

			// loop through k neibhors to find if there is a shorter dist
			for(int j = 0; j < k; j++){

				if(*(tempResDis + s*k + idxRes) < *(tempRes.ndist + s*k + idxNres)){
					*(res.ndist + s*k + j) = *(tempResDis + s*k + idxRes);
					res.nidx[s*k + j] = *(tempResIDs + s*k + idxRes);
					idxRes++;
				}else{
					*(res.ndist + s*k + j) = *(tempResDis + s*k + idxNres);
					res.nidx[s*k + j] = tempRes.nidx[s*k + idxNres];
					idxNres++;
				}
			}
		}
	}

	return res;
}

