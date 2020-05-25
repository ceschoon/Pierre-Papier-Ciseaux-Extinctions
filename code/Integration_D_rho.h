////////////////////////////////////////////////////////////////////////////
//                                                                        //
//    This file is part of a project on extinctions in a cyclic           //
//    predation model like the "rock-paper-scissors" game.                //
//                                                                        //
//    It contains source code to compute D(rho) by numerical integration. //
//                                                                        //
//    Author: Cedric Schoonen                                             //
//                                                                        //
////////////////////////////////////////////////////////////////////////////


#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

using namespace std;


// dynamics
double at(double a, double b, double c) {return a*(b-c);}
double bt(double a, double b, double c) {return b*(c-a);}
double ct(double a, double b, double c) {return c*(a-b);}


void step_euler(double &a, double &b, double &c, double dt)
{
	double a_new = a + at(a,b,c) * dt;
	double b_new = b + bt(a,b,c) * dt;
	double c_new = c + ct(a,b,c) * dt;
	
	a = a_new;
	b = b_new;
	c = c_new;
}


double find_max_a(double rho)
{
	// The max value of a occurs when b=c is minimum. We find the value
	// of a iteratively until we have a value obeying
	//    rho = a*b*c = a*b^2 = 1/4 * a*(1-a)^2
	
	double a = 1.0;
	double diff = rho - 1.0/4 * a*(1-a)*(1-a);
	double step = -1e-2;
	double precision = 1e-8; 
	
	while (abs(step)>precision)
	{
		double a_new = a + step;
		double diff_new = rho - 1.0/4 * a_new*(1-a_new)*(1-a_new);
		
		// If we go in the wrong direction,
		// we turn back with a lower step
		
		if (abs(diff_new) > abs(diff))
		{
			step *= -1;
			step /= 2;
		}
		
		// In any case, we update the values of a and diff
		
		a = a_new;
		diff = diff_new;
	}
	
	return a;
}


double compute_D_rho(double rho)
{
	// The strategy is to integrate the dynamics numerically and 
	// record D(a,b,c) at each step. D(rho) is computed on-the-fly by
	// averaging D(a,b,c) along the trajectory. 
	
	// We start with values for a,b,c such that a is maximum. The value of
	// a will thus decrease initially. When a starts increasing, it means 
	// we have completed half a period of the deterministic orbit. We only
	// do half a period because the function D(a,b,c) is symmetric under
	// the exchange of b and c.
	
	double t = 0;
	double dt = 1e-3;
	double a = find_max_a(rho);
	double b = 1.0/2 * (1-a);
	double c = 1.0/2 * (1-a);
	double D_rho = 0;
	double T_rho = 0;
	bool completed_half_period = false;
	
	while (!completed_half_period)
	{
		double a_old = a;
		step_euler(a,b,c,dt);
		t += dt;
		
		double D_abc = rho*rho*(-9+1/a+1/b+1/c);
		D_rho += D_abc * dt;
		
		if (a>a_old) 
		{
			completed_half_period = true;
			T_rho = 2.0*t;
		}
	}
	
	D_rho *= 2.0/T_rho;
	
	# ifdef DEBUG_DRHO
	cout << endl;
	cout << "--- Completed D(rho) computation ---" << endl;
	cout << "rho =   " << rho << endl;
	cout << "a+b+c = " << a+b+c << endl;
	cout << "a*b*c = " << a*b*c << endl;
	cout << "T_rho = " << T_rho << endl;
	cout << "D_rho = " << D_rho << endl;
	cout << endl;
	# endif
	
	return D_rho;
}


void debug()
{
	#ifndef DEBUG_DRHO
	#define DEBUG_DRHO
	#endif
	
	double rho = 0.02;
	double a_max = find_max_a(rho);
	double D_rho = compute_D_rho(rho);
	
	cout << "rho =   " << rho << endl;
	cout << "a_max = " << a_max << endl;
	cout << "D_rho = " << D_rho << endl;
}



