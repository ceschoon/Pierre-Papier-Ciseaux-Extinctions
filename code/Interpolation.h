////////////////////////////////////////////////////////////////////////////
//                                                                        //
//    This file is part of a project on extinctions in a cyclic           //
//    predation model like the "rock-paper-scissors" game.                //
//                                                                        //
//    It contains source code to construct the interpolation of a 1D      //
//    function. The goal is to speed up computations when the function    //
//    is computed many times.                                             //
//                                                                        //
//    Author: Cedric Schoonen                                             //
//                                                                        //
////////////////////////////////////////////////////////////////////////////


#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

using namespace std;

class Interpolation
{
	public:
		
		Interpolation(double x_min, double x_max, double dx)
		{
			x_min_ = x_min;
			x_max_ = x_max;
			dx_ = dx;
			
			N_ = int ( (x_max-x_min)/dx );
			x_ = vector<double>(N_,0);
			y_ = vector<double>(N_,0);
		}
		
		void init(double (*f)(double))
		{
			for (int i=0; i<N_; i++)
			{
				x_[i] = x_min_ + dx_*i;
				y_[i] = f(x_[i]);
			}
		}
		
		double eval(double x)
		{
			// find closer index in x_
			
			int i = int ( (x-x_min_)/dx_ );
			
			// check that it is not too close to the boundaries
			
			if (i<1 || i>N_-2)
			{
				cout << "ERROR: in interpolation.h: evaluation at point too"
				     << " close to interval boundaries" 
				     << endl;
				
				cout << "       x =      " << x << endl;
				cout << "       x_min_ = " << x_min_ << endl;
				cout << "       x_max_ = " << x_max_ << endl;
				cout << "       dx_ =    " << dx_ << endl;
			}
			
			// nearby data points
			
			double x1 = x_[i-1];
			double x2 = x_[i];
			double x3 = x_[i+1];
			
			double y1 = y_[i-1];
			double y2 = y_[i];
			double y3 = y_[i+1];
			
			// quadratic interpolation
			
			double det = (x1-x3)*(x3-x2)*(x2-x1);
			
			double a = - 1/det * (        (x3-x2)*y1 +         (x1-x3)*y2 +         (x2-x1)*y3);
			double b =   1/det * ((x3+x2)*(x3-x2)*y1 + (x1+x3)*(x1-x3)*y2 + (x2+x1)*(x2-x1)*y3);
			double c = - 1/det * (  x3*x2*(x3-x2)*y1 +   x1*x3*(x1-x3)*y2 +   x2*x1*(x2-x1)*y3);
			
			return a*x*x + b*x + c;
		}
		
	protected:
	
		double x_min_;
		double x_max_;
		double dx_;
		
		int N_;
		vector<double> x_;
		vector<double> y_;
};

