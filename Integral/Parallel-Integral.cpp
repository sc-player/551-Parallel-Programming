#include<iostream>
#include<string>
#include<climits>
#include<cmath>
#include<cstdlib>
#include<mpi.h>
using namespace std;

#define trueValue 4003.7209001513268265
#define STOPCRIT 0.000000000000005

void GetData(int comm_sz, int rank, long double* a_ptr, long double* b_ptr, long int* n_ptr); 
long double Trap(long double local_a, long double local_b, long double local_n, long double h);
long double f(long double x);

int main(){
	long int n;
  long double b, a;
	long double error;
  int comm_sz;
  int rank;
  long double aproxValue;
  long double h;	
  long int local_n;
  long double local_a;
  long double local_b;
  long double total_approx;
  double startTime;
  double endTime;
  long int remainder;
  MPI_Init(NULL, NULL);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  GetData(comm_sz, rank, &a, &b, &n);
  h=(b-a)/n;
	local_n=n/((long int)comm_sz);
  remainder=n%((long int)comm_sz);
  if(rank<remainder){
    local_a=a+((long int)rank)*(local_n+1)*h;
    local_b=local_a+local_n*h;
    local_n++;
  } else {
    local_a=a+((long int)rank)*local_n*h+remainder*h;
    local_b=local_a+(local_n-1)*h;
  } 
  MPI_Barrier(MPI_COMM_WORLD); 
  startTime=MPI_Wtime(); 
	aproxValue=Trap(local_a, local_b, local_n, h);  
  MPI_Reduce(&aproxValue, &total_approx, 1, MPI_LONG_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  endTime=MPI_Wtime(); 
  if(rank==0){
	  error=fabs((trueValue-total_approx)/trueValue);
    cout.precision(6);
    cout << std::scientific << endl << "Running on " << comm_sz << " cores." << endl;
    cout << "Elapsed time = " << endTime-startTime<< " seconds" << endl;
    cout << "With n = " << n << " trapezoids, our estimate\nof the integral from "<< std::fixed << a << " to " << b << " = ";
    cout.precision(13);
    cout << std::scientific << total_approx << endl;
    cout.precision(19);
	  cout << "true value = " << trueValue << endl;
	  cout << "absolute relative true error: " << error << endl;
	  cout << " is " << ((error<STOPCRIT)?"":"NOT ") << "less than criteria = " << STOPCRIT << endl; 
  }
  MPI_Finalize();
}

void GetData(int p, int rank, long double* a_ptr, long double* b_ptr, long int* n_ptr){
  if(rank==0){
    cout << "Enter a, b, and n\n";
    cin >> *a_ptr;
    cin >> *b_ptr;
    cin >> *n_ptr;
  } 
  MPI_Bcast(a_ptr, 1, MPI_LONG_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(b_ptr, 1, MPI_LONG_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(n_ptr, 1, MPI_LONG_DOUBLE, 0, MPI_COMM_WORLD);
}    

long double Trap(long double local_a, long double local_b, long double local_n, long double h){
  long double res=(f(local_a)+f(local_b))/2.0l;
  for(int i=1; i<=local_n-1; i++){
    res+=f(local_a+i*h);
  }
  res=res*h;
  return res;
}

long double f(long double x)
{
	return cos(x/3.0l)-2.0l*cos(x/5.0l)+5.0l*sin(x/4.0l)+8.0l;
}
