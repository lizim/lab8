// Lorenz.cxx
// Runge Kutta 4th order
// GL, 4.12.2015
//--------------------
#include <cmath>
#include <iostream>
#include <fstream>
//--------------------
void f(double* const y0);
void RKstep(double* const yn, const double* const y0, const double x, const double dx, double* k1, double* k2, double* k3, double* k4);
double bisection( double y0, double xn, double dx, double ynp1 ,double* k1, double* k2, double* k3, double* k4);
//--------------------
using namespace std;
//--------------------

int main(void)
{ 
	//ofstream out("solution");
	const int dim = 2;
	double dx = 0.1;
	double k1[dim], k2[dim], k3[dim], k4[dim];
	const double L = 100;

	double yn[dim];
  ofstream out1("Periode");	
for (int i = 1; i <50; i++){
	double y0[dim] = {0.1*i, 0};
	double ynp1 = y0[1]; // temporary y2 to compare
	//out << x << "\t" << y0[0] << "\t" << y0[1] << endl;
	double x=0;
	while(x<=L)  // Beim Abbruch dieser schleife ist die Ableitung (y2) ynp1 <0 und y0 (= yn ) >0
	{
		x += dx;
		
		RKstep(yn, y0, x, dx, k1, k2, k3, k4);
		ynp1= yn[1];
		if (ynp1 < 0 && y0[1] >0 ){break;}
		  for(int i=0; i<dim; i++) {y0[i] = yn[i];}
		//out << x << "\t" << y0[0] << "\t" << y0[1] << endl;
	}
	double period = bisection( y0[1], x, dx, ynp1 ,k1, k2, k3,k4);
	
	
	out1 <<y0[0]<< "\t" <<period << endl; // das funzt hier, weil yn nach einer periode wieder y0 ist
	
  
}
  //out.close();
  out1.close();
	return(0);
}
//-------------------
void RKstep(double* const yn, const double* const y0,
            const double x, const double dx, double* k1, double* k2, double* k3, double* k4)
{
	const int dim = 2;
	

  for(int i=0;i<dim; i++) k1[i] = y0[i];
	f(k1);

	for(int i=0;i<dim; i++) k2[i] = y0[i] + 0.5 * dx * k1[i];
	f(k2);

	for(int i=0;i<dim; i++) k3[i] = y0[i] + 0.5 * dx * k2[i];
	f(k3);

  for(int i=0;i<dim; i++) k4[i] = y0[i] + dx * k3[i];
	f(k4);

	for(int i=0;i<dim; i++)
	 yn[i] = y0[i] + 1./6.*dx*(k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);
}
//-------------------
// motion in E field
void f(double* const y0)
{
  double ytemp = y0[0];
  y0[0] = y0[1];
  y0[1] = - ytemp/(sqrt(1+ ytemp*ytemp));
  
}
double bisection( double y0, double xn, double dx, double ynp1 ,double* k1, double* k2, double* k3, double* k4){
  double yn = y0;
 
  double thetal = 0;
  double thetar = 1;
  double theta = (thetal + thetar) /2;
  
  double ym=1;
  while (abs(ym) > 1e-5){
      theta = (thetal + thetar) /2;
      double b1 = theta - 3.0/2.0 *theta*theta + 2.0/3.0 * theta* theta* theta;
      double b2 = theta*theta - 2.0/3.0* theta* theta* theta;
      double b3 = b2;
      double b4 = -0.5* theta*theta + 2.0/3.0 * theta* theta* theta;
      
      ym = y0 + dx*(b1*k1[1] + b2*k2[1] +b3*k3[1] +b4*k4[1]);
       
    if (ym < 0) 
      thetar = theta ; 
    
    else 
      thetal = theta ; 
      
      
    
    
   }
    double xm = (xn+ theta*dx);
  return xm;
}