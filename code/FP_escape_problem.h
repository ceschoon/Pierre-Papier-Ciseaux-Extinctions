////////////////////////////////////////////////////////////////////////////
//                                                                        //
//    This file is part of a project on extinctions in a cyclic           //
//    predation model like the "rock-paper-scissors" game.                //
//                                                                        //
//    It contains source code to compute the mean escape time from an     //
//    interval, given a homogeneous 1D stochastic process.                //
//                                                                        //
//    Author: Cedric Schoonen                                             //
//                                                                        //
////////////////////////////////////////////////////////////////////////////


#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

using namespace std;


// The stochastic process must obey a 1D homogeneous Fokker-Planck equation
//    d/dt p(x,t) = - d/dx( A(x) p(x,t) ) + 1/2 d^2/dx^2( B(x) p(x,t) )

// x=a is an absorbing boundary
// x=b is a reflective boundary 
// We assume a<b

// For more information, see Gardiner's Handbook of stochastic methods,
// section 5.2

// It assumes that the drift and diffusion coefficients A,B are defined 
// somewhere else

double A(double x);
double B(double x);

double a_boundary;
double b_boundary;
double integration_dx;


double integrate(double (*f)(double), double x_min, double x_max, double dx)
{
	// Numerical integration using trapezoid method
	// Makes an error quadratic in dx
	
	double I = 0;
	double f0 = f(x_min);
	double f1;
	
	for (double x=x_min; x<x_max; x+=dx) 
	{
		f1 = f(x+dx);
		I += dx * (f0+f1)/2.0;
		f0 = f1;
	}
	
	return I;
}


double f_phi(double x) {return 2*A(x)/B(x);}

double compute_phi(double x)
{
	#ifdef DEBUG_ESCAPE_LVL3
		cout << "Computing phi ..." << endl;
	#endif
	
	double phi_log = integrate(f_phi,a_boundary,x,integration_dx);
	
	#ifdef DEBUG_ESCAPE_LVL3
		cout << " computed phi(" << x << ") = " << exp(phi_log) << endl;
	#endif
	
	return exp(phi_log);
}


double f_exit1(double z) 
{
	#ifdef DEBUG_ESCAPE_LVL3
		cout << "Computing f_exit1 ..." << endl;
	#endif
	
	double phi_z = compute_phi(z);
	
	#ifdef DEBUG_ESCAPE_LVL3
		cout << " computed f_exit1(" << z << ") = " << phi_z/B(z) << endl;
	#endif
	
	return phi_z/B(z);
}

double f_exit2(double y)
{
	#ifdef DEBUG_ESCAPE_LVL2
		cout << "Computing f_exit2 ..." << endl;
	#endif
	
	double phi_y = compute_phi(y);
	
	double I = integrate(f_exit1,y,b_boundary,integration_dx);
	
	#ifdef DEBUG_ESCAPE_LVL2
		cout << " computed f_exit2(" << y << ") = " << I/phi_y << endl;
	#endif
	
	return I/phi_y;
}

double compute_exit_time(double x)
{
	#ifdef DEBUG_ESCAPE_LVL1
		cout << "Computing mean exit time ..." << endl;
	#endif
	
	double I = integrate(f_exit2,a_boundary,x,integration_dx);
	
	#ifdef DEBUG_ESCAPE_LVL1
		cout << " computed t_exit(" << x << ") = " << 2*I << endl;
	#endif
	
	return 2*I;
}


double debug_func(double x) {return x*x;}

void debug_integration()
{
	double x_min = 0;
	double x_max = 1;
	double dx = 0.01;
	
	double integral_ana = 1.0/3;
	double integral_num = integrate(debug_func,x_min,x_max,dx);
	double integral_err = abs(integral_ana-integral_num);
	
	cout << endl;
	cout << "--- debugging numerical integration ---" << endl;
	cout << "integrating y=x^2 from x=0 to x=1 by steps of " << dx << endl;
	cout << "integral_ana = " << integral_ana << endl;
	cout << "integral_num = " << integral_num << endl;
	cout << "integral_err = " << integral_err << endl;
	cout << endl;
}



