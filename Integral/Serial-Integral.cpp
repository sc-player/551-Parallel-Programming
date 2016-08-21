//T found: 10500410

#include<iostream>
#include<string>
#include<climits>
#include<cmath>
#include<cstdlib>
#include<mpi.h>
using namespace std;

#define trueValue 4003.7209001513268265
#define STOPCRIT 0.000000000000005

long double f(long double x);

int main(int argc, char**argv){
  cout.precision(20);

	long double t=1, b=600.0l, a=100.0l;
	long double error = 1.0l;
	long double h;
	long double aproxValue;
  int digit=0;
  int exp=0;
  while((error>STOPCRIT)||(exp>=0)){
    h=(b-a)/t;
    aproxValue=f(a)+f(b);
	  for(int i=1; i<=t-1; ++i){
		  aproxValue+= 2.0l*f(a+i*h);
  	}
	  aproxValue=h/2.0l*aproxValue;
	  error=abs((trueValue-aproxValue)/trueValue);
    if(digit==0){
      if(error>STOPCRIT){
        t=t*10;
        ++exp;
      }
      else
      {
        t=t/10;
        digit+=1;
        exp-=2;
      }
    } 
    else{
      if(error>STOPCRIT){
        t+=pow(10,exp);
      }
      else{
        t-=pow(10,exp);
        exp-=1;
      }
    }
    cout << "Current t estimation: " << t<< endl;
  }

	cout << "True: " << trueValue << endl;
	cout << "Approximate: " << aproxValue << endl;
	cout << "Absolute Relative True Error: " << fabs(error) << endl;
	cout << "Stopping Criteria: " << STOPCRIT << endl;
  cout << "Min trapezoids required (t):" << t << endl;
}

long double f(long double x)
{
	return cos(x/3.0l)-2.0l*cos(x/5.0l)+5.0l*sin(x/4.0l)+8.0l;
}
