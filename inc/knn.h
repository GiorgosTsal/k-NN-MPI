// Definition of the kNN result struct
typedef struct knnresult {
    int* nidx;		//!< Indices (0-based) of nearest neighbors [m - by - k]
    double* ndist;	//!< Distance of nearest neighbors			[m - by - k]
    int m;			//!< Number of query points					[scalar]
    int k;			//!< Number of nearest neighbors			[scalar]
} knnresult;

    
    
    
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
knnresult kNN(double* X, double* Y, int n, int m, int d, int k);

void populateArray(double *array, int n, int d);
void printArray(double *array, int n, int d);
void printIndicesArray(int * array, int size);
int * indexArray(int size);
int partitionWithIndex(double *arr,int *idx, int l, int r);
int partition(double *arr, int l, int r);
double kthSmallestWithIndex(double *arr,int *idx,int l, int r, int k);
double kthSmallest(double *arr, int l, int r, int k);
void swap_d(double *p, double *q);
void swap_i(int *p, int *q);
