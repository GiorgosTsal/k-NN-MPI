#include <stdio.h>
#include <string.h>
#include "mpi.h"

int main(int argc, char** argv){
  int pid;
  int nproc;
  int n = 4, d = 2, m = 4, k=2;
  char message[100];
  MPI_Status status;
  int tag = 50;
 double *X = (double *)malloc(n * d * sizeof(double));
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &pid);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  if (pid = 0) {
    sprintf(message, "Hello from process %d!", pid);

    printf("%s\n",message);

    int dest = 0;
    MPI_Send(X, n*d+1, MPI_DOUBLE, 1,
	     tag, MPI_COMM_WORLD);
  } else {
    
   
      
      MPI_Recv(X, n*d, MPI_DOUBLE, MPI_ANY_SOURCE, tag,
	       MPI_COMM_WORLD, &status);
      
      printf("Received >>%s<<\n", message);
 
     
    
  }
  MPI_Finalize();

  return(0);
}
