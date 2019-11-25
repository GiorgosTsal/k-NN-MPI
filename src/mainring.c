// Example using MPI_Send and MPI_Recv to pass a message around in a ring.

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include "../inc/knn.h"

//to anevazo gia na to do apla spiti, mi tou dineis simasia

const int n = 10, d = 3, k=2, p = 5, m = 2;

int main(int argc, char** argv) {
	printf("\n \n===============INIT===================\n \n");
	// Initialize the MPI environment
	MPI_Init(&argc, &argv); 
	
	int chunk_size = n/p;
	printf("to chunk ine: %d \n", chunk_size);
	double *data;
	
	// Find out rank, size
	int process_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &process_rank); //the current ID of the running process
	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size); //total number of running MPI processes
	printf("to size einai: %d \n", world_size);
	printf("to process_rank  einai: %d \n", process_rank);

	if (process_rank == 0)
	{
		printf("When process rank is 0...\n");
		 for (int pr= 0; pr < p; pr++){
			data = (double *)malloc(chunk_size * d * sizeof(double));
			populateArray(data, chunk_size, d);
			printf("Array data is: \n");
			printArray(data, chunk_size, d);
			printf("\n");
	
			int dest = pr+1;
			int tag = 0;
			
			if (pr == p-1)            // last chunk is mine
			break;
			// send to correct process
			printf("Send to correct process...\n");
			MPI_Send(data, chunk_size*d, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
			printf("Send to correct process done...\n");
			free(data);
		 }
	} else {                      // ..... other processes
	printf("When process rank is Not 0...\n");
	// from which process to I receive (master)? what tag?
	int recv = 0;
	int tag = 0;
	MPI_Status Stat;
	printf("3\n");
	data = (double * ) malloc( chunk_size*d * sizeof(double) );
	MPI_Recv(data, chunk_size*d, MPI_DOUBLE, recv, tag, MPI_COMM_WORLD, &Stat);
	printf("4\n");

	}

	//// Receive from the lower process and send to the higher process. Take care
	//// of the special case when you are the first process to prevent deadlock.
	//if (process_rank != 0) {
		//printf("When process rank not 0...\n");
		//MPI_Recv(&X, chunk_size*d, MPI_INT, process_rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		//printf("2\n");
		//printf("Process %d received token %lf from process %d\n", process_rank, *X, process_rank - 1);
	//} else {
		//printf("nothing, already setted...\n");
	//// Set the token's value if you are process 0
		
	//}
	//printf("Before Send to next process...\n");
	//MPI_Send(&X, 1, MPI_DOUBLE, (process_rank + 1) % world_size, 0, MPI_COMM_WORLD);
	//printf("Send  to next process done...\n");
	//// Now process 0 can receive from the last process. This makes sure that at
	//// least one MPI_Send is initialized before all MPI_Recvs (again, to prevent
	//// deadlock)
	//if (process_rank == 0) {
		//printf("When process rank is 0...\n");
		//MPI_Recv(&X, chunk_size*d, MPI_DOUBLE, world_size - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		//printf("Recieved...\n");
		//printf("Process %d received token %lf from process %d\n", process_rank, *X, world_size - 1);
		//printf("Ok...\n");
	//}

	printf("Finalizing...\n");
	MPI_Finalize();
	
	return 0;
}
