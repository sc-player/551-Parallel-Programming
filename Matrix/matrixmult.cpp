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
  cout << "Enter form: 'ikj', 'kij', or 'ijk'.\n";
  cin >> f;
  if(!((strcmp(f, "ijk")==0)||(strcmp(f, "kij")==0)||(strcmp(f,"ikj")==0))){
    cerr << "Invalid form.\n";
    exit(1);
  }

  //Flag
  cout << "Type I to input matrices, R to randomly generate.\n";
  cin >> flg;
  if(*flg!='R' && *flg!='I'){
    cerr << "Invalid flag.\n";
    exit(1);
  }

  //Size
  cout << "Enter size of one side of the matrices.\n";
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
    cout << "A\n";
    for(int i=0; i<(*n)*(*n); ++i){
      cin >> (*aIn)[i];
    }
    cout << "B\n";
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
  int* sendCounts; //Local sizes of each process.
  int* displ;
  //MPI process information variables.
  int comm_sz; 
  int rank;
  int remainder;
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
  //b is required for all proccesses. c is the where we gather the solution,
  //so it is not required.
  if(rank!=0){
    b=(int*)malloc(n*n*sizeof(int)); 
  }
  else{
    c=(int*)malloc(n*n*sizeof(int));
  }

  //Calculate local sizes.
  sendCounts = (int*)malloc(comm_sz*sizeof(int));
  displ=(int*)malloc(comm_sz*sizeof(int));
  if(rank==0){
    remainder=n%comm_sz;
    displ[0]=0;
    for(int i=0; i<comm_sz; ++i){
      sendCounts[i]=n/comm_sz;
      if(i<remainder){
        sendCounts[i]++;
      }
      sendCounts[i]=sendCounts[i]*n;
      if(i>0){
        displ[i]=displ[i-1]+sendCounts[i-1];
      }
    }
    for(int i=1; i<comm_sz; ++i){
      displ[i]=sendCounts[i-1]+displ[i-1];
    }
  }
  MPI_Bcast(sendCounts, comm_sz, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(displ, comm_sz, MPI_INT, 0, MPI_COMM_WORLD);

  //Allocate local arrays.
  local_a=(int*)malloc(sendCounts[rank]*sizeof(int));
  local_c=(int*)malloc(sendCounts[rank]*sizeof(int));
  for(int i=0; i<sendCounts[rank]; ++i){
    local_a[i]=0;
    local_c[i]=0;
  }
  
  //Broadcast b, as it is required by all proccesses, and scatter a.
  MPI_Bcast(b, n*n, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Scatterv(a, sendCounts, displ, MPI_INT, local_a, sendCounts[rank], MPI_INT, 0, MPI_COMM_WORLD);

  //Calculate local solution using the specified form.
  if(strcmp(form, "ijk")==0){
    for(int i=0; i<sendCounts[rank]/n; ++i){
      for(int j=0; j<n; ++j){
        for(int k=0; k<n; ++k){
          local_c[i*n+j]+=local_a[i*n+k]*b[k*n+j]; 
        }
      }
    }
  }
  else if(strcmp(form,"ikj")==0){
    for(int i=0; i<sendCounts[rank]/n; ++i){
      for(int k=0; k<n; ++k){
        for(int j=0; j<n; ++j){
          local_c[i*n+j]+=local_a[i*n+k]*b[k*n+j];
        }
      }
    }
  }
  else if(strcmp(form, "kij")==0){
    for(int k=0; k<n; ++k){
      for(int i=0; i<sendCounts[rank]/n; ++i){
        for(int j=0; j<n; ++j){
          local_c[i*n+j]+=local_a[i*n+k]*b[k*n+j];
        }
      }
    }
  }

  //Gather local results to proccess 0, and print timing/solution.
  MPI_Gatherv(local_c, sendCounts[rank], MPI_INT, c, sendCounts, displ, MPI_INT, 0, MPI_COMM_WORLD);  
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
  free(sendCounts);
  if(rank==0) free(a);
  free(b);
  if(rank==0) free(c);
  MPI_Finalize();
}
