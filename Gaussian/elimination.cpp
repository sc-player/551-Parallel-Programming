//elimination.cpp
//Philip Strachan
//A parallel program that will solve a matrix of size n through Gaussian 
//Elimination with partial pivoting. It stores the matrix A in an array of 
//pointers to arrays. This allows us to swap the matrix rows easily.
//Usage: elimination <n> <threads>

#include<iostream>
#include<cstdlib>
#include<cmath>
#include<omp.h>
using namespace std;

//main
//parameters:
//  argc: must be 3
//  argv: elimination, n, threads
int main(int argc, char*argv[]){

  //Declare variables
  double** original;
  double* boriginal;
  double** A; 
  double* b;
  double* r;
  double swap;
  double* swapPtr;
  int n;
  int high;
  double highVal;
  double lSquaredNorm;
  double conversionFactor;
  int threads;
  double startTime;

  //Check usage and get the command line arguments.
  if(argc!=3){
    cout << "Usage: elimination <Matrix Size> <Amount of Threads>\n";
    exit(1);
  }
  n=atoi(argv[1]);
  threads=atoi(argv[2]);

  //Display information and seed generator.
  cout << "Number of cores available: " << omp_get_num_procs() << endl << "Number of threads: " << threads << endl;
  srand48(time(NULL));

  //Allocate all matrix rows and vectors.
  original=new double*[n];
  A=new double*[n];
  boriginal=new double[n];
  b=new double[n];
  r=new double[n];

  //Allocate matrix columns and generate numbers.
  for(int i=0; i<n; ++i){ 
    A[i]=new double[n];
    original[i]=new double[n];
    for(int j=0; j<n; ++j){
      A[i][j]=original[i][j]=drand48()*2*pow(10,6)-pow(10,6);
    }
    b[i]=boriginal[i]=drand48()*2*pow(10,6)-pow(10,6);
    r[i]=0;
  }

  //Take beginning time and begin elimination.
  startTime=omp_get_wtime();
  for(int i=0; i<n-1; ++i){
/*  Enable to print matrix. Don't on large matrices.

    cout << "Randomly Generated Numbers\n";
    for(int ind=0; ind<n; ++ind){
      for(int jnd=0; jnd<n; ++jnd){
        cout << A[ind][jnd] << " ";
      }
      cout << b[ind] << endl;     
    }*/
    
    //Partial pivoting.
    high=i;
    highVal=A[i][i];;
    for(int k=i; k<n; ++k){
      if(abs(A[k][i])>highVal){
        high=k;
        highVal=abs(A[k][i]);
      }
    }
    if(high!=i){
      swapPtr=A[i];
      A[i]=A[high];
      A[high]=swapPtr;
      swap=b[i];
      b[i]=b[high];
      b[high]=swap;
    }

    //Begin parallel elimination of the column. Each thread will take care
    //of one row in the column denoted by i.  My conversionfactor must be 
    //private to avoid overwriting.
#   pragma omp parallel for num_threads(threads) private(conversionFactor)
    for(int j=i+1; j<n; ++j){
      conversionFactor=A[j][i]/A[i][i];
      for(int k=i; k<n; ++k){
        A[j][k]-=A[i][k]*conversionFactor;
        if(k==i) A[j][k]=0;
      }
      b[j]-=b[i]*conversionFactor;
    }
  }

  //Back substitution.
  for(int i=n-1; i>=0; --i){    
    for(int j=n-1; j>i; --j){
      b[i]-=A[i][j]*b[j];
    }
    b[i]=b[i]/A[i][i];
  } 

  //Print the final time.
  cout << "Time elapsed: " << omp_get_wtime()-startTime << endl;

  //Calculate residuals.
  for(int i=0; i<n; ++i){
    r[i]=0;
    for(int j=0; j<n; ++j){
      r[i]+=original[i][j]*b[j];
    }
    r[i]-=boriginal[i];
  }

  //Calculate the lsquared norm, the square root of the sum of the squared 
  //residuals, to check correctness.
  lSquaredNorm=0;
  for(int i=0; i<n; ++i){
    lSquaredNorm+=pow(r[i], 2);
  }
  lSquaredNorm=sqrt(lSquaredNorm);
  cout << "lSquaredNorm: " << lSquaredNorm << endl;
}
