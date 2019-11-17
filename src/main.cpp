
// CPP code for bubble sort 
// using template function 
#include <iostream> 
#include "../inc/knn.h"

using namespace std; 


const int n = 10, d = 2;


// Driver code 
int main() 
{ 
	
	double *X = (double *)malloc(n * d * sizeof(double));
	populateArray(X, n, d);
	printArray(X, n, d);

	int *idx = (int *)malloc(n *sizeof(int));
	idx = indexArray(n);
	printIndicesArray(idx, n);
	
   
    ///*Testing Point*/
    //Point p; 
    //p.x = 2.5; 
    //p.y = 9; 
  
    //// Parameter to decide group of the testing point 
    //int k = 1; 
    //printf ("The value classified to unknown point"
            //" is %d.\n", classifyAPoint(arr, n, k, p)); 
    return 0; 
} 


