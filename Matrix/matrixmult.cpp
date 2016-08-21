/* Philip Strachan
 * matrixmult.cpp
 * A parallel MPI program to calculate a matrix.
 * Allows use of ijk, ikj, or kij forms.
 * First array accessing rows is scattered
 * across processes with each process receiving
 * n/p rows.
 * The second array is broadcasted, as it is accessed
 * by columns and all columns are required by each row.
 * The result is computed into a local_c array of the
 * same size as local_a (n/p). 
 * Finally, the result is gathered into the master process.
 */


//Includes================================================================
#include<string>
#include<iostream>
#include<stdio.h>
#include<cstdlib>
#include<ctime>
#include<cstring>
#include<unistd.h>
#include<mpi.h>
using namespace std;
//========================================================================



/* Philip Strachan
 * void get_data
 * parameters:
 *   char *f: pointer to form array (ijk, kij, ikj), decides which form 
 *   to use.
 *   char *flg: pointer to flag character (I, R), decides to input or 
 *   generate matrix.
 *   int* n: pointer to the size of the matrix.
 *   int** aIn: pointer to the first matrix.
 *   int** bIn: pointer to the second matrix.
 * This function reads the form, flag, and size from standard input.
 * It then allocates the arrays and  generates or reads in the matrices 
 * according to the flag.
 * All parameters are passed in by referance in order to modify them.
 */
void get_data(char* f, char* flg, int* n, int** aIn, int** bIn){
  
  //Form
  cin >> f;
  if(!((strcmp(f, "ijk")==0)||(strcmp(f, "kij")==0)||(strcmp(f,"ikj")==0))){
    cerr << "Invalid form.\n";
    exit(1);
  }

  //Flag
  cin >> flg;
  if((*flg!='R')&&(*flg!='I')){
    cerr << "Invalid flag.\n"; 
    exit(1);
  }

  //Size
  cin >> *n;
  if(*n<1){
    cerr <<"Invalid size.\n";
    exit(1);
  }

  //Allocate arrays
  *aIn=(int*)malloc((*n)*(*n)*sizeof(int));
  *bIn=(int*)malloc((*n)*(*n)*sizeof(int));

  //If flag is 'I', read in arrays
  if(*flg=='I'){
    for(int i=0; i<(*n)*(*n); ++i){
      cin >> (*aIn)[i];
    }
    for(int i=0; i<(*n)*(*n); ++i){
      cin >> (*bIn)[i];
    }
  }

  //If flag is 'R', generate arrays
  else if(*flg=='R'){
    srand(time(NULL));
    for(int col=0; col<(*n); ++col){
      for(int row=0; row<(*n); ++row){
        (*aIn)[row+col*(*n)]=rand()%100;
        (*bIn)[row+col*(*n)]=rand()%100;
      }
    }
  }
}

/* Philip Strachan
 * int main
 * no parameters
 * Main function ran in parallel. Rank 0 process will read input, and then
 * starts timing. The variables and matrix b are broadcasted, while matrix
 * a is scattered into evenly divided blocks of rows. Depending on the 
 * form string, the local result is calculated using inner, middle, or
 * outer product. Finally, the result is gathered to rank 0, where the 
 * timing information is printed out, and the result if the matrices were
 * inputted.
 */
int main(){
  char*form; //algorithm to use
  char flag; //input or generate matrices?
  int n; //size of matrices (nxn)
  //matrices: a*b=c
  int*a;
  int*b;
  int*c;
  //local portions of a and c.
  int*local_a;
  int*local_c;
  int local_n; //local size (rows of a and c)
  //MPI process information variables.
  int comm_sz; 
  int rank;
  double startTime; //Time

  //Initialize MPI and grab variables
  MPI_Init(NULL, NULL);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
  form=(char*)malloc(3);
  //Get data from input if rank is 0, otherwise wait
  if(rank==0){
    get_data(form, &flag, &n, &a, &b);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  
  //Start time for rank 0
  if(rank==0) startTime=MPI_Wtime();
  
  //Broadcast size and algorithm
  MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(form, 4, MPI_CHAR, 0, MPI_COMM_WORLD);

  //Allocate b if we are not rank 0, otherwise we are rank 0 so allocate c.
  //b is required for all proccesses. c is the where we gather the solution, so it is not required.
  if(rank!=0){
    b=(int*)malloc(n*n*sizeof(int)); 
  }
  else{
    c=(int*)malloc(n*n*sizeof(int));
  }

  //Allocate local arrays
  local_n=n/comm_sz;
  local_a=(int*)malloc(local_n*n*sizeof(int));
  local_c=(int*)malloc(local_n*n*sizeof(int));
  for(int i=0; i<local_n*n; ++i){
    local_a[i]=0;
    local_c[i]=0;
  }

  //Broadcast b, as it is required by all proccesses, and scatter a.
  MPI_Bcast(b, n*n, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Scatter(a, n*local_n, MPI_INT, local_a, n*local_n, MPI_INT, 0, MPI_COMM_WORLD);

  //Calculate local solution using the specified form.
  if(strcmp(form, "ijk")==0){
    for(int i=0; i<local_n; ++i){
      for(int j=0; j<n; ++j){
        for(int k=0; k<n; ++k){
          local_c[i*n+j]+=local_a[i*n+k]*b[k*n+j]; 
        }
      }
    }
  }
  else if(strcmp(form,"ikj")==0){
    for(int i=0; i<local_n; ++i){
      for(int k=0; k<n; ++k){
        for(int j=0; j<n; ++j){
          local_c[i*n+j]+=local_a[i*n+k]*b[k*n+j];
        }
      }
    }
  }
  else if(strcmp(form, "kij")==0){
    for(int k=0; k<n; ++k){
      for(int i=0; i<local_n; ++i){
        for(int j=0; j<n; ++j){
          local_c[i*n+j]+=local_a[i*n+k]*b[k*n+j];
        }
      }
    }
  }

  //Gather local results to proccess 0, and print timing/solution.
  MPI_Gather(local_c, n*local_n, MPI_INT, c, n*local_n, MPI_INT, 0, MPI_COMM_WORLD); 
  if(rank==0){
    cout << "running on " << comm_sz << " processors\n" << "elapsed time = "<< MPI_Wtime()-startTime << " seconds\n";
    if(flag=='I'){
      for(int col=0; col<n; ++col){
        for(int row=0; row<n; ++row){
          cout << c[col*n+row] << " ";
        }
        cout << endl;
      }
    }
  }

  //Clean up 
  free(local_a);
  free(local_c);
  free(form);
  if(rank==0) free(a);
  free(b);
  if(rank==0) free(c);
  MPI_Finalize();
}
