
// C++ program to find groups of unknown 
// Points using K nearest neighbour algorithm. 

#include <bits/stdc++.h> 
#include "../inc/knn.h"

using namespace std; 

//Populate array/dataset with random double values
void populateArray(double *array, int n, int d){
	  cout << __func__ << endl;
	  for (int i = 0; i < n; i++) {
			for (int j = 0; j < d; j++) {
				*(array + i*d + j) = rand()/double(RAND_MAX)*24.f+1.f;	//rand() % 100;
			}
	  }
}
//Helper to visualize results
void printArray(double *array, int n, int d)
{
    cout << __func__ << endl;
    for (int i = 0; i < n; ++i)
    {
        cout << i << ": ";
        for (int j = 0; j < d; ++j)
            cout <<  *(array + i*d + j)<< endl;
       cout << endl;
    }
}
//Print Indices of an array
void printIndicesArray(int *array, int size)
{
    cout << __func__ << endl;
    for (int i = 0; i < size; ++i)
    {
            cout <<  *(array + i)<< endl;
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
 
double dist(){
	
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
	//// Fill distances of all points from p 
    for (int i = 0; i < n; i++) 
        dista = sqrt((arr[i].x - p.x) * (arr[i].x - p.x) + (arr[i].y - p.y) * (arr[i].y - p.y)); 
}


//This function calculates distances for all points 
double calcDistanceMatrix(double *distance, double *x, double *y, int n, int m){
	double *distances = (double *)malloc(size *sizeof(double));
	for(int i=0; i<size; i++){
		*(distances + i) = 
		*(distance+i) = this->calculateDistance(this->data+(*(index+i))*d);
	}
}

// This function calculates euclidean distances between d points
double calculateDistance(double *xRow, double *yRow, int d)
{
	double sum = 0;
	for (int i = 0; i < d; i++) {
		sum = sum + pow(*(xRow+i) - *(yRow+i), 2);
	}
	return sqrt(sum);
}



//// Used to sort an array of points by increasing 
//// order of distance 
//bool comparison(Point a, Point b) 
//{ 
    //return (a.distance < b.distance); 
//} 
  
//// This function finds classification of point p using 
//// k nearest neighbour algorithm. It assumes only two 
//// groups and returns 0 if p belongs to group 0, else 
//// 1 (belongs to group 1). 
//int classifyAPoint(Point arr[], int n, int k, Point p) 
//{ 
    //// Fill distances of all points from p 
    //for (int i = 0; i < n; i++) 
        //arr[i].distance = 
            //sqrt((arr[i].x - p.x) * (arr[i].x - p.x) + 
                 //(arr[i].y - p.y) * (arr[i].y - p.y)); 
  
    //// Sort the Points by distance from p 
    //sort(arr, arr+n, comparison); 
  
    //// Now consider the first k elements and only 
    //// two groups 
    //int freq1 = 0;     // Frequency of group 0 
    //int freq2 = 0;     // Frequency of group 1 
    //for (int i = 0; i < k; i++) 
    //{ 
        //if (arr[i].val == 0) 
            //freq1++; 
        //else if (arr[i].val == 1) 
            //freq2++; 
    //} 
  
    //return (freq1 > freq2 ? 0 : 1); 
//} 
  

