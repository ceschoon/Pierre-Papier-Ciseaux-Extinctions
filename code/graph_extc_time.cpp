////////////////////////////////////////////////////////////////////////////
//                                                                        //
//    This file is part of a project on extinctions in a cyclic           //
//    predation model like the "rock-paper-scissors" game.                //
//                                                                        //
//    Its purpuse is to make a graph of the extinction time vs rho.       //
//                                                                        //
//    Compile with $ g++ -O3 -o graph_extc_time graph_extc_time.cpp       //
//                                                                        //
//    Author: Cedric Schoonen                                             //
//                                                                        //
////////////////////////////////////////////////////////////////////////////


#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>
#include <stdlib.h>
#include "Integration_D_rho.h"
#include "FP_escape_problem.h"
#include "Interpolation.h"

using namespace std;



/////////////////////////// global parameters //////////////////////////////

int N = 1000;
double D_approx = 1e-3;

// the precision of the interpolation should be big enough compared to
// than the precision of the integrations in the escape problem
// more precisely:   dx_interp <= 1e4 * dx_integr^2    (see below)
Interpolation D_rho_interpol(0, 1.0/27, 1e-6); 

double A_exact(double rho) {return -3.0/N * rho;}

double B_comput(double rho) {return  1.0/N * compute_D_rho(rho);}
double B_interp(double rho) {return  1.0/N * D_rho_interpol.eval(rho);} 
double B_approx(double rho) {return  1.0/N * D_approx;}

////////////////////////////////////////////////////////////////////////////




int main(int argc, char **argv)
{
	//////////////////// option for estimation of D(rho) ///////////////////
	
	// default
	A = &A_exact;
	B = &B_interp;
	
	for (int i=1; i<argc; i++)
	{
		if ((string(argv[i])=="--comput")) B = &B_comput;
		if ((string(argv[i])=="--interp")) B = &B_interp;
		if ((string(argv[i])=="--approx")) B = &B_approx;
	}
	
	//////////////////////////// other options /////////////////////////////
	
	for (int i=1; i<argc; i++)
	{
		if ((string(argv[i])=="--N")) 
		{
			if (argc > i+1) N = stoi(string(argv[i+1]));
			else 
				throw runtime_error("You must specify the value of N, e.g. --N 100");
		}
		
		if ((string(argv[i])=="--D_approx")) 
		{
			if (argc > i+1) D_approx = stod(string(argv[i+1]));
			else 
				throw runtime_error("You must specify the value of D_approx, e.g. --D_approx 1e-3");
		}
	}
	
	////////////////////////// global parameters ///////////////////////////
	
	a_boundary = 0;          // absorbing barrier
	b_boundary = 1.0/27;     // reflective barrier
	integration_dx = 4e-5;   // discretisation length
	
	// we crop the interval to avoid problems due to D(a)=D(b)=0
	// as there are places where we divide by D.
	// also, the interpolation we use for D cannot be used too close to 
	// the boundaries
	// we crop the interval by an amount that is quadratic in dx in 
	// order to save the quadratic precision of the integration algorithm.
	
	a_boundary +=   integration_dx*integration_dx*1e4;
	b_boundary -= 2*integration_dx*integration_dx*1e4;
	
	// contruct interpolation of D(rho)
	
	D_rho_interpol.init(compute_D_rho);
	
	//////////// We compute T(rho) for several values of rho ///////////////
	
	//double tau = compute_exit_time(0.02);
	//cout << "exit time / N = " << tau/N << endl;
	
	double rho_min = 0;
	double rho_max = 1.0/27;
	double rho_step = 1e-3;
	
	// shift to avoid problems near boundaries
	rho_min +=   rho_step;
	rho_max -= 2*rho_step;
	
	int sysresult = system("mkdir -p extc_time");
	
	ofstream dataFile("extc_time/data.dat");
	dataFile << fixed << setprecision(8);
	dataFile << "# rho      extc_time" << endl;
	
	for (double rho=rho_min; rho<rho_max; rho+=rho_step)
	{
		double tau = compute_exit_time(rho);
		dataFile << rho << " " << tau/N << endl;
	}
	
	/////////////////////////// Gnuplot script /////////////////////////////
	
	ofstream plotFile("extc_time/plot.gp");
	
	plotFile << "#----------------------------------------" << endl;
	plotFile << "set terminal epslatex standalone color size 9cm,7cm" << endl;
	plotFile << "set output 'plot.tex'" << endl;
	plotFile << endl;
	plotFile << "#set title '$T(\\rho)$'" << endl;
	plotFile << "set xlabel '$\\rho$'" << endl;
	plotFile << "set ylabel '$T(\\rho)$'" << endl;
	plotFile << "set style data lines" << endl;
	plotFile << "set key off" << endl;
	plotFile << endl;
	plotFile << "set xtics 0.01" << endl;
	plotFile << "set ytics 0.1" << endl;
	plotFile << "set grid" << endl;
	plotFile << endl;
	plotFile << "plot 'data.dat'" << endl;
	plotFile << endl;
	
	plotFile.close();
	
	//////////////////// Bash script to plot everything ////////////////////
	
	ofstream shellFile("extc_time/plot.sh");
	
	shellFile << "#! /bin/bash" << endl;
	shellFile << "" << endl;
	shellFile << "gnuplot   plot.gp" << endl;
	shellFile << "" << endl;
	shellFile << "pdflatex  plot.tex" << endl;
	shellFile << "" << endl;
	shellFile << "rm   plot.aux   plot.log   plot.tex   plot-inc.eps   plot-inc-eps-converted-to.pdf" << endl;
	shellFile << "" << endl;
	shellFile << "mv plot.pdf extc_time.pdf" << endl;
	shellFile << "" << endl;
	
	shellFile.close();
	
	// execute it
	sysresult = system("chmod +x extc_time/plot.sh; cd extc_time; ./plot.sh");
	
	
	return 0;
}

