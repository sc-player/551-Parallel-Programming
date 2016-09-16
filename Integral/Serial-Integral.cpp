/* Philip Strachan
 * Serial-Integral.cpp
 * A program that calculates the minimum number of trapezoids required to 
 * estimate the true value of a function with a certain degree of error. The
 * function and all variables are hardcoded. The true value of the function
 * and the stopping criteria error are defined below, so as to be easily 
 * changeable if the start and end times are changed or if you would like a 
 * different criteria.
 *
 * T found for error=5*10^-15: 10500410
 */

//Includes================================================================
#include<iostream>
#include<string>
#include<climits>
#include<cmath>
#include<cstdlib>
#include<mpi.h>
using namespace std;
//========================================================================

//Defined constants
#define trueValue 4003.7209001513268265
#define STOPCRIT 0.000000000000005

//Function to estimate.
long double f(long double x);

/* Philip Strachan
 * int main
 * no parameters
 * We take the hardcoded values and approximate the integral with
 * increasing numbers of trapezoids. While the error is larger than the
 * stopping criteria, the exponent will count up, starting at 1. Once the 
 * error is smaller than the stopping criteria, the exponent is lowered by 
 * 1, and then we begin increasing the top digit. When we have again passed
 * the stopping criteria, we will lower the exponent again by one and 
 * increase this digit. Repeating this process, we will find the minimum
 * number of trapezoids required to reach this stopping criteria.
 */
int main(){

  //Current trapezoid estimation, endpoint, and startpoint.
	long double t=1, b=600.0l, a=100.0l;
	long double error = 1.0l; //Current error.
	long double h; //Size of one trapezoid.
	long double aproxValue; //Current approximate value.
  bool digit=false; //Have we reached the maximum digit?
  int exp=0; //What digit are we looking at currently?

  //While the error is greater than criteria or we are not at the smallest
  //digit...
  cout.precision(20);
  while((error>STOPCRIT)||(exp>=0)){

    //Calculate length of one trapezoid (x direction).
    h=(b-a)/t;

    //Calculate approximate value for this number of trapezoids.
    aproxValue=f(a)+f(b);
	  for(int i=1; i<=t-1; ++i){
		  aproxValue+= 2.0l*f(a+i*h);
  	}
	  aproxValue=h/2.0l*aproxValue;

    //Calculate error.
	  error=abs((trueValue-aproxValue)/trueValue);

    //If we have not reached the max digit yet...
    if(!digit){

      //If we are still above the criteria, add another exponent.
      if(error>STOPCRIT){
        t=t*10;
        ++exp;
      }

      //Otherwise, go back to the last exponent and let the program know 
      //that we have passed criteria.
      else
      {
        t=t/10;
        digit=true;
        exp-=2;
      }
    } 

    //If we have reached the max digit...
    else{

      //If we haven't reached the criteria, add another to that digit.
      if(error>STOPCRIT){
        t+=pow(10,exp);
      }
      
      //Otherwise, subtract one from the digit and go to the next lowest
      //digit.
      else{
        t-=pow(10,exp);
        exp-=1;
      }
    }

    //Print current estimation.
    cout << "Current t estimation: " << t<< endl;
  }

  //Print Results.
	cout << "True: " << trueValue << endl;
	cout << "Approximate: " << aproxValue << endl;
	cout << "Absolute Relative True Error: " << fabs(error) << endl;
	cout << "Stopping Criteria: " << STOPCRIT << endl;
  cout << "Min trapezoids required (t):" << t << endl;
}

//Hardcoded function.
long double f(long double x)
{
	return cos(x/3.0l)-2.0l*cos(x/5.0l)+5.0l*sin(x/4.0l)+8.0l;
}
