
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
        for (int j = 0; j < d; ++j){
            printf("%f\n", *(array + i*d + j));
			putchar('\n');
			fflush(stdout);
		}
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
int partition(double *arr, int *idx,int l, int r){
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
// This function returns k'th smallest
// element in arr[l..r] using QuickSort
// based method.  ASSUMPTION: ALL ELEMENTS
// IN ARR[] ARE DISTINCT
double kthSmallest(double *arr, int *idx,int l, int r, int k){
    // If k is smaller than number of
    // elements in array
    if (k > 0 && k <= r - l + 1) {

        // Partition the array around last
        // element and get position of pivot
        // element in sorted array
        int index = partition(arr, idx, l, r);

        // If position is same as k
        if (index - l == k - 1)
            return *(arr + index);

        // If position is more, recur
        // for left subarray
        if (index - l > k - 1)
            return kthSmallest(arr, idx,l, index - 1, k);

        // Else recur for right subarray
        return kthSmallest(arr, idx,index + 1, r,
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
	double * Tdistance = (double *)malloc(m * n  *sizeof(double));
	double *xRow = (double *)malloc(m * n  *sizeof(double));
	double *yRow = (double *)malloc(m * sizeof(double));
	int *indexes = (int*)malloc(m * n  *sizeof(int));

	

	for(int i=0; i<m; i++){
		for(int j=0; j<n; j++){
			*(indexes+i*n+j)=j;
			printf("indexes1: %d \n", *(indexes+i*n+j));
		}
	}
	/*
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
			*(distance + i*m + j) += xRow[i] + yRow[j];
			*(distance + i*m + j) = sqrt( *(distance + i*m + j) );
			printf("Distances1: %d %f \n",i*m + j,*(distance + i*m + j));
		}
		
	}
	*/
	// Calculates euclidean distance.
	for(int i=0;i<n;i++){
		for(int j=0;j<m;j++){
			double dist=0;
			for(int k=0;k<d;k++){
				
				dist+=(X[i*d+k]-Y[j*d+k])*(X[i*d+k]-Y[j*d+k]);
				
				
			}
			distance[i*m+j]=sqrt(dist);
		
		}
	
	
	}
	
	
	//	After that, only keep the k shortest distances and the indices of the corresponding points (using k-select like in the
	//	previous assignment?).

	
	printf("\n=== Distances === \n" );
	for(int i=0; i<n; i++){
		for(int j=0; j<m; j++){
			printf("%f " , *(distance + i*m + j));
		}
	printf("\n" );
	}
	
	printf("\n" );
	double *kNeighbours = (double *)malloc(k * sizeof(double));
	for(int i=1; i<=k; i++){
	
		printf("oi %d neibhor ine: %f ",i , kthSmallest(distance, indexes,0, m*n-1, i));
		//loipon eho teleiosei nomizo me to distances alla des to k esu ligo genika dn eimai 100% sigouros oti einai sosto
		//i kthSmallest pou eho gyrnaei double ta k mikrotera distances edo. autos thelei na gurname k ta index
		//https://www.geeksforgeeks.org/how-to-return-multiple-values-from-a-function-in-c-or-cpp/
		//me tis texnikes edo logika pame sto na epistrefei domi i an eheis esu kamia kaluteri lusi ennoeite
		//kai nomizo etsi teleionei to seiriako alla ksanaleo dn eimai sigouros oti ta apotelesmata pou vgazo einai sosta
		//genika an einai check mia ton tropo se sxesi me ton typo D pou dinei gia na vroume tis apostaseis k an pisteueis oti einai sostos apla krata ta indexes k teleiose to seiriako
		*(kNeighbours + i) = kthSmallest(distance, indexes, 0, m-1, i); //edo allaksa tin kthSmallest(prosthesa os parametro to indexes-pou tha einai ta index tou pinaka distance 
																		//k ekana k mia swap gia integers pou tha xreiastei sigoura an pame me auti ti logiki)
		
		printf("\n");
		//na to valo na gyrnaei k ta indices tous-kamia domi gia to dist??
	}
	

	
}





