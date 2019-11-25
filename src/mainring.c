// Example using MPI_Send and MPI_Recv to pass a message around in a ring.

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include "../inc/knn.h"

//to anevazo gia na to do apla spiti, mi tou dineis simasia

const int n = 23, d = 3, k=2, p = 5;

int main(int argc, char** argv) {
	// Initialize the MPI environment
	MPI_Init(NULL, NULL);
	
	
	
	
	
	//int rem = n%p;
	//int quan = d*p;
//	printf("To kvanto einai: %d kai h nea nRow diastash einai %d kai to modulo einai %d \n", quan, nRow, rem);
	double *X = (double *)malloc(n * d * sizeof(double));
	populateArray(X, n, d);
	printf("Array X is: ");
	printf("\n");
	printArray(X, n, d);
	printf("\n");
	int chunk_size = n*d/p;
	int bonus = n*d - chunk_size * p;
	printf("To chunk_size einai: %d kai to bonus einai %d \n", chunk_size, bonus);
	if (bonus != 0)
	{
		for (int i=0; i < p + 1; i++)
		{
			printf("Array X_%d is: ", i);
			printf("\n");
			if (p == i){
				printf("edo 1 kai ta rows einai: %d kai to bonus/d einai : \n", bonus);
				printArray(X + chunk_size*i -bonus*i, bonus*bonus, d); //pattern to use for splitting array
			}
			else{
				printf("edo 2\n");
				printArray(X + chunk_size*i -2*i, n/p, d); //pattern to use for splitting array
			}
			
		}
	}else
	{
		for (int i=0; i < p; i++)
		{
			printf("Array X_%d is: ", i);
			printf("\n");
			printArray(X + chunk_size*i , n/p, d); //pattern to use for splitting array
		}
	}
	

	
	
	
	//if (rem == 0)
	//{
		//for (int i = 0; i < p; i++)
		//{
			//printf("Array X_%d is: ", i);
			//printf("\n");
			//printArray(X + i*quan + d*i, nRow, d); //pattern to use for splitting array
		//}
	//}else
	//{
		//for (int i = 0; i < p+1; i++)
		//{
			//if (i != p)
			//{
				//printf("Array X_%d is: ", i);
				//printf("\n");
				//printArray(X + i*quan + d*i, nRow, d); 
			//}else
			//{
				//printf("Array X_%d is: ", i);
				//printf("\n");
				//printArray(X + i*quan + d*i, rem, d); //pattern to use for splitting array
			//}
		//}
	//}
	
	
	
	
	// Find out rank, size
	int process_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &process_rank); //the current ID of the running process
	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size); //total number of running MPI processes
	printf("to size einai: %d \n", world_size);

	int token;
	// Receive from the lower process and send to the higher process. Take care
	// of the special case when you are the first process to prevent deadlock.
	if (process_rank != 0) {
	MPI_Recv(&token, 1, MPI_INT, process_rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	printf("Process %d received token %d from process %d\n", process_rank, token, process_rank - 1);
	} else {
	// Set the token's value if you are process 0
	token = -1;
	}
	MPI_Send(&token, 1, MPI_INT, (process_rank + 1) % world_size, 0, MPI_COMM_WORLD);
	// Now process 0 can receive from the last process. This makes sure that at
	// least one MPI_Send is initialized before all MPI_Recvs (again, to prevent
	// deadlock)
	if (process_rank == 0) {
	MPI_Recv(&token, 1, MPI_INT, world_size - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	printf("Process %d received token %d from process %d\n", process_rank, token, world_size - 1);
	}
	MPI_Finalize();
}
