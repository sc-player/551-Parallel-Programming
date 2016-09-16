/* Philip Strachan
 * Parallel-Integral.cpp
 * A parallel implementation of the trapezoidal rule for estimating
 * integrals. The function and all variables are hardcoded. The true value
 * of the function and the stopping criteria error are defined below, so as
 * to be easily changeable if the start and end times are changed or if you 
 * would like a different criteria.
 */

//Includes===============================================================
#include<iostream>
#include<string>
#include<climits>
#include<cmath>
#include<cstdlib>
#include<mpi.h>
using namespace std;
//=======================================================================

//Defined constants
#define trueValue 4003.7209001513268265
#define STOPCRIT 0.000000000000005

//Forward Declarations
void GetData(long double* a_ptr, long double* b_ptr, long int* n_ptr); 
long double Trap(long double local_a, long double local_b, long double local_n, long double h);
long double f(long double x);

/* Philip Strachan
 * int main
 * no parameters
 * After initializing MPI and grabbing communicator variables, we query the
 * user for data if we are process 0, otherwise we wait. We start the timing
 * after the data is entered and broadcast the start and end points and the
 * number of trapezoids. Each process calculate the length of the trapezoid,
 * and then each process calculates its start and end points and local size.
 * Each process sums the areas of its local trapezoids, and then the rank 0 
 * process reduces the approximate values into a global sum. The final time
 * is taken and the rank 0 process prints out the timing information.
 */
int main(){

  //Number of Trapezoids.
	long int n;

  //End, Start.
  long double b, a;

  //Calculated relative true error.
	long double error;

  //Number of processes.
  int comm_sz;

  //MPI Proccess ID.
  int rank;

  //Current calculated approximate value.
  long double aproxValue;

  //Size of one trapezoid.
  long double h;	

  //Number of trapezoids to each process.
  long int local_n;

  //Local start.
  long double local_a;

  //Local end.
  long double local_b;

  //Final reduced approximate value.
  long double total_approx;

  //Timings.
  double startTime;
  double endTime;

  //Leftover trapezoids when divided between processes.
  long int remainder;

  //Initialize MPI and grab variables
  MPI_Init(NULL, NULL);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if(rank==0) GetData(&a, &b, &n);

  //Start Timing
  MPI_Barrier(MPI_COMM_WORLD);
  startTime=MPI_Wtime();

  //Broadcast a, b, and n.
  MPI_Bcast(&a, 1, MPI_LONG_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&b, 1, MPI_LONG_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&n, 1, MPI_LONG_DOUBLE, 0, MPI_COMM_WORLD);

  //Calculate trapezoid length (x axis).
  h=(b-a)/n;

  //Calculate local size and leftover trapezoids.
	local_n=n/((long int)comm_sz);
  remainder=n%((long int)comm_sz);
  
  //Allocate remaining trapezoids, and calculate local sizes.
  if(rank<remainder){
    local_n++;
    local_a=a+local_n*rank*h;
  } else {
    local_a=a+(remainder*(local_n+1)+(rank-remainder)*local_n)*h;
  }
  local_b=local_a+local_n*h;
 
  //Calculate local approximate.
	aproxValue=Trap(local_a, local_b, local_n, h);  

  //Reduce local approximates into a total approximate on proccess 0.
  MPI_Reduce(&aproxValue, &total_approx, 1, MPI_LONG_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  //Take end time.
  endTime=MPI_Wtime(); 

  //Print results.
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

/* Philip Strachan
 * void GetData
 * parameters:
 *   a_ptr: pointer to where start point is stored.
 *   b_ptr: pointer to where end point is stored.
 *   n_ptr: pointer to where the number of trapezoids is stored.
 *   Parameters passed by reference are set using cin.
 */
void GetData(long double* a_ptr, long double* b_ptr, long int* n_ptr){
  cout << "Enter a, b, and n\n";
  cin >> *a_ptr;
  cin >> *b_ptr;
  cin >> *n_ptr;
}    

/* Philip Strachan
 * long double Trap
 * parameters:
 *   long double local_a: starting point
 *   long double local_b: ending point
 *   long double local_n: trapezoids to calculate
 *   long double h: width of a trapezoid (x axis)
 * A function that estimates the value of an integral using the trapezoidal
 * method from local_a to local_b, using local_n trapezoids of h width.
 * Returns the estimation of the integral from local_a to local_b.
 */
long double Trap(long double local_a, long double local_b, long double local_n, long double h){

  //Start with the height of the start and endpoints divided by two.
 // long double res=0;
  long double res=(f(local_a)+f(local_b))/2.0l;
  //Sum all points in between.
  for(int i=1; i<local_n; i++){
    res+=f(local_a + i*h);
  }

  //Multiply times the width.
  res=res*h;
  return res;
}

//Hardcoded function
long double f(long double x)
{
	return cosl(x/3.0l)-2.0l*cosl(x/5.0l)+5.0l*sinl(x/4.0l)+8.0l;
}
