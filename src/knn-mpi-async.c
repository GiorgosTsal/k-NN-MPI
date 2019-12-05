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

		//	printf("Distance of Y[ %d] and X[ %d] is: %lf \n",i,j,distance[i*n+j]);

		}

	}

	

	knnres.ndist=(double *)malloc(m*k * sizeof(double));

	knnres.nidx=(int *)malloc(m*k * sizeof(int));

	//Calculates the minimum distance of each point of y from X dataset

	for(int i=0;i<m;i++){

	//	printf("\n Gia i= %d : \n",i);

		for(int j=0;j<k;j++){

	//		printf("MPAINEI");

			knnres.ndist[i*k+j]=kthSmallestWithIndex(&distance[i*n], &indexes[i*n],0, n-1, j+1);

			

			knnres.nidx[i*k+j]=indexes[i*n+j];

			

	//		printf("O geitonas %d exei apostasi: %lf \n",knnres.nidx[i*k+j],knnres.ndist[i*k+j]);

			

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

	int to = (process_rank+1)%pr;          // Sending source (next one)

	int from = (process_rank-1+pr)%pr; // Receiving source (previous)

	// Y matrix for receiving points
	double* Y =  (double *)malloc(n * d  *sizeof(double));
	double* Ytemp =  (double *)malloc(n * d  *sizeof(double));

	// MPI variables used for async communication
    MPI_Request reqSend, reqReceive;
    
	// sending data asynchronously
    MPI_Isend(X, n*d, MPI_DOUBLE, to, tag, MPI_COMM_WORLD, &reqSend);
    MPI_Irecv(Ytemp, n*d, MPI_DOUBLE, from, tag, MPI_COMM_WORLD, &reqReceive);




	int branks = process_rank;




	// run the kNN code to get the first result

	double * Xtemp = (double *)malloc(n * d  *sizeof(double));

	memcpy(Xtemp, X, n * d * sizeof(double));

	knnresult res = kNN(Xtemp, X, n, n, d, k);





	for(int i=0; i<n; i++){

		for(int j=0;j<k; j++){

			      *(res.nidx+i*k+j) = *(res.nidx+i*k+j)+n*((process_rank-1+pr)%pr);

		}

	}


	knnresult tempRes;



	//double * Y = (double *)malloc(n * d  *sizeof(double));

	memcpy(Y, X, n * d * sizeof(double));

	//double * Ytemp = (double *)malloc(n * d  *sizeof(double));



	int idxRes = 0;

	int idxNres = 0;

	double* tempResDis = (double *)malloc(n * k  *sizeof(double));

	int* tempResIDs = (int *)malloc(n * k  *sizeof(int));


	// Variables to hold minimum and maximum reduction
	double min, max, mintmp, maxtmp;


	for(int i = 0; i < pr-1; i++){


		MPI_Wait(&reqSend, MPI_STATUS_IGNORE);
		MPI_Wait(&reqReceive, MPI_STATUS_IGNORE);
		// if even, first send and then receive

		knnresult knnTemp;

		MPI_Isend(Y, n*d, MPI_DOUBLE, to, tag, MPI_COMM_WORLD, &reqSend);
		MPI_Irecv(Y, n*d, MPI_DOUBLE, from, tag, MPI_COMM_WORLD, &reqReceive);
		

		knnTemp = kNN(Y,X,n,n,d,k); 

		for(int x=0; x<n; x++){

			for(int j=0; j<k; j++){

				*(knnTemp.nidx+x*k+j) = *(knnTemp.nidx+x*k+j) + ((process_rank-2-i+pr)%pr)*n;

				

				

				if(*(res.ndist+x*k+k-1)>*(knnTemp.ndist+x*k+j)){

					*(res.ndist+x*k+k-1)=*(knnTemp.ndist+x*k+j);

					*(res.nidx+x*k+k-1)=*(knnTemp.nidx+x*k+j);

					for(int q=0;q<k;q++){

						kthSmallestWithIndex(&res.ndist[x*k], &res.nidx[x*k],0, k-1, q+1);

					}

				}

			}

		}

		free(Y);

		free(Ytemp);


	

	}
	
	// Find min max
	min = res.ndist[1];
	max = res.ndist[k-1];

	//reduce min and max
	for(int i = 0; i < n; i++){
		if(res.ndist[k*i] < min){
			min = res.ndist[k*i + 1];
		}	
		if(res.ndist[k*i + k - 1] > max){
			
			max = res.ndist[k*i + k - 1];
		}
	}
	

	MPI_Reduce(&min, &mintmp, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
	MPI_Reduce(&max, &max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);


	
	printf("Minimum distance is %f, maximum is %f\n", min, max);


	return res;

}

