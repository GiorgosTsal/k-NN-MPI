/*!
  \file   tester.c
  \brief  Validate kNN ring implementation.
  \author Dimitris Floros
  \date   2019-11-13
*/


#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <assert.h>
#include "../inc/knn.h"

#include "tester_helper.h"

int main()
{

  int n=12235;                    // corpus
  int m=456;                    // query 
  int d=11;                      // dimensions
  int k=3;                     // # neighbors

  double  * corpus = (double * ) malloc( n*d * sizeof(double) );
  double  * query  = (double * ) malloc( m*d * sizeof(double) );

  for (int i=0;i<n*d;i++)
    corpus[i] = ( (double) (rand()) ) / (double) RAND_MAX;

  for (int i=0;i<m*d;i++)
    query[i]  = ( (double) (rand()) ) / (double) RAND_MAX;

	double  total_time;
	struct timeval start, end;
	gettimeofday (&start, NULL);

  knnresult knnres = kNN( corpus, query, n, m, d, k );
	
	gettimeofday (&end, NULL);
	total_time = (double)((end.tv_usec - start.tv_usec)/1.0e6 + end.tv_sec - start.tv_sec);
	printf("Time taken by kNN function: : %f sec\n",total_time);

  int isValidC = validateResultColMajor( knnres, corpus, query, n, m, d, k );

  int isValidR = validateResultRowMajor( knnres, corpus, query, n, m, d, k );
  
  printf("Tester validation: %s NEIGHBORS\n",
         STR_CORRECT_WRONG[isValidC||isValidR]);

  free( corpus );
  free( query );

  return 0;
  
}
