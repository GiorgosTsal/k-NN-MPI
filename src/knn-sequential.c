
// C program to find groups of unknown 
// Points using K nearest neighbour algorithm. 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include "../inc/knn.h"
#include <cblas.h>


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
            printf("%f\n", *(array + i*d + j));
			putchar('\n');
			fflush(stdout);
		}
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
	
	
	//test calculation
	double * dist = (double *)malloc(m * n  *sizeof(double));
	for (int j = 0; j < m; j++ ){
		for (int i = 0; i < n; i++ ){
			double dista = 0;
			for (int l = 0; l < d; l++){
			dista += ( X[l*n+i] - Y[l*m+j] ) * ( X[l*n+i] - Y[l*m+j] );
			}
			*(dist+i*n+j) = sqrt(dista);
		}	
	}
	for (int j = 0; j < m; j++ ){
			for (int i = 0; i < n; i++ ){
			printf("to dianysmaa distances einai: %lf \n", *(dist +i*n+j));
		}
		printf("\n");
	}
	
	
	
	// Calculates euclidean distance.
	for(int i=0;i<m;i++){
		for(int j=0;j<n;j++){
			double dist=0;
			for(int k=0;k<d;k++){
				dist+=(X[j*d+k]-Y[i*d+k])*(X[j*d+k]-Y[i*d+k]);		
			}
			distance[i*n+j]=sqrt(dist);
		}
	}

	
	
		
	double * tempDis = (double *)malloc(m *sizeof(double));
	//Calculates the minimum distance of each point of Y from X dataset
	for(int i=0;i<m;i++){
		printf("makis is: %d \n", distance[i*n]); 
		*(tempDis + i)=kthSmallest(&distance[i*n], 0, n-1, 1);
		if (*(tempDis + i) == 0)
		{
			*(tempDis + i)=kthSmallest(&distance[i*n], 0, n-1, 2);
		}
		printf("tempDis[ %d ] =  %lf \n",i,tempDis[i]);
		indexes[i]=i;
	}
	
	free(dist);
	
	knnres.ndist=(double *)malloc(k * sizeof(double));
	knnres.nidx=(int *)malloc(k * sizeof(int));
	
	//Calculates k nearest neighbors
	for(int i=0; i<k; i++){
		
		knnres.ndist[i]=kthSmallestWithIndex(tempDis, indexes,0, m-1, i+1);
		knnres.nidx[i]=indexes[i];
		printf("oi %d neibhor ine: %lf ",i , knnres.ndist[i] );
		printf("To id ine: %d ", knnres.nidx[i] );
		printf("\n");
	}
	free(tempDis);
	free(indexes);
	return knnres;	
}





