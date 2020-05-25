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

int N = 100;
double D_rho_approx = 0.001;

Interpolation D_rho_interpol(0, 1.0/27, 1e-6); 

double A(double rho) {return -3.0/N * rho;}
double B(double rho) {return  1.0/N * D_rho_interpol.eval(rho);} 
//double B(double rho) {return  1.0/N * compute_D_rho(rho);}
//double B(double rho) {return  1.0/N * D_rho_approx;}

////////////////////////////////////////////////////////////////////////////

int main()
{
	////////////////////////// global parameters ///////////////////////////
	
	a_boundary = 0;          // absorbing barrier
	b_boundary = 1.0/27;     // reflective barrier
	
	integration_dx = 8e-5;
	a_boundary += integration_dx;    // reduce interval to avoid problems
	b_boundary -= 2*integration_dx;  // because D(a)=D(b)=0 and there are 1/D
	
	// contruct interpolation of D(rho)
	
	D_rho_interpol.init(compute_D_rho);
	
	//double tau = compute_exit_time(0.02);
	//cout << "exit time / N = " << tau/N << endl;
	
	//////////// We compute D(rho) for several values of rho ///////////////
	
	double rho_min = 0;
	double rho_max = 1.0/27;
	double rho_step = 1e-3;
	
	// shift to avoid problems near boundaries
	rho_min += rho_step;
	rho_max -= 2*rho_step;
	
	int sysresult = system("mkdir -p exit_time");
	
	ofstream dataFile("exit_time/data.dat");
	dataFile << fixed << setprecision(8);
	dataFile << "# rho      exit_time" << endl;
	
	for (double rho=rho_min; rho<rho_max; rho+=rho_step)
	{
		double tau = compute_exit_time(rho);
		dataFile << rho << " " << tau/N << endl;
	}
	
	/////////////////////////// Gnuplot script /////////////////////////////
	
	ofstream plotFile("exit_time/plot.gp");
	
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
	
	ofstream shellFile("exit_time/plot.sh");
	
	shellFile << "#! /bin/bash" << endl;
	shellFile << "" << endl;
	shellFile << "gnuplot   plot.gp" << endl;
	shellFile << "" << endl;
	shellFile << "pdflatex  plot.tex" << endl;
	shellFile << "" << endl;
	shellFile << "rm   plot.aux   plot.log   plot.tex   plot-inc.eps   plot-inc-eps-converted-to.pdf" << endl;
	shellFile << "" << endl;
	shellFile << "mv plot.pdf exit_time.pdf" << endl;
	shellFile << "" << endl;
	
	shellFile.close();
	
	// execute it
	sysresult = system("chmod +x exit_time/plot.sh; cd exit_time; ./plot.sh");
	
	
	return 0;
}
