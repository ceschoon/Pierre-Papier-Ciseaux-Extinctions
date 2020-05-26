////////////////////////////////////////////////////////////////////////////
//                                                                        //
//    This file is part of a project on extinctions in a cyclic           //
//    predation model like the "rock-paper-scissors" game.                //
//                                                                        //
//    Its purpuse is to make the phase portrait of the deterministic      //
//    dynamics.                                                           //
//                                                                        //
//    Compile with $ g++ -O3 -o phase_portrait phase_portrait.cpp         //
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

using namespace std;


int main()
{
	///////// We generate trajectories for different values of rho /////////
	
	// parameters
	
	double T = 20;
	double dt = 1e-3;
	
	double rho_min = 0;
	double rho_max = 1.0/27;
	double rho_step = 5e-3;
	
	rho_min += rho_step;
	
	// declarations
	
	int sysresult = system("mkdir -p portrait");
	
	ofstream dataFile("portrait/data.dat");
	dataFile << fixed << setprecision(8);
	dataFile << "# x_traj1  y_traj1 x_traj2  y_traj2 ...  " << endl;
	
	int Ntraj = int((rho_max-rho_min)/rho_step);
	
	vector<double> rho(Ntraj,0);
	for (int i=0; i<Ntraj; i++) rho[i] = rho_min + i*rho_step;
	
	vector<double> a(Ntraj,0);
	vector<double> b(Ntraj,0);
	vector<double> c(Ntraj,0);
	
	// identify starting point
	
	for (int i=0; i<Ntraj; i++)
	{
		a[i] = find_max_a(rho[i],3);
		b[i] = 1.0/2 * (1-a[i]);
		c[i] = 1.0/2 * (1-a[i]);
	}
	
	// generate trajectory
	
	double t=0;
	
	while (t<T)
	{
		// time step 
		
		for (int i=0; i<Ntraj; i++)
		{
			double a0 = a[i];
			double b0 = b[i];
			double c0 = c[i];
			
			step_euler(a0,b0,c0,dt);
			
			a[i] = a0;
			b[i] = b0;
			c[i] = c0;
		}
		
		t += dt;
		
		// transform (a,b,c) to (x,y) coordinates in the (1,1,1) plane
		// and save the position in data file
		
		for (int i=0; i<Ntraj; i++)
		{
			double x = (c[i]-b[i]) * sqrt(3)/2;
			double y = a[i]-b[i]/2-c[i]/2;
			
			dataFile << x << " " << y << " ";
		}
		
		dataFile << endl;
	}
	
	/////////////////////////// Gnuplot script /////////////////////////////
	
	ofstream plotFile("portrait/plot.gp");
	
	plotFile << "#----------------------------------------" << endl;
	plotFile << "set terminal epslatex standalone color size 5cm,5cm" << endl;
	plotFile << "set output 'plot.tex'" << endl;
	plotFile << endl;
	plotFile << "#set title ''" << endl;
	plotFile << "#set xlabel '$\\frac{\\sqrt{3}}{2} (c-b)$'" << endl;
	plotFile << "#set ylabel '$a-c/2-b/2$'" << endl;
	plotFile << "set key off" << endl;
	plotFile << endl;
	plotFile << "unset xtics" << endl;
	plotFile << "unset ytics" << endl;
	plotFile << "unset border" << endl;
	plotFile << endl;
	plotFile << "set xrange [-0.85:0.80]" << endl;
	plotFile << "set yrange [-0.70:1.15]" << endl;
	plotFile << endl;
	plotFile << "set label 'a' at -0.05,1.1" << endl;
	plotFile << "set label 'b' at -1.03,-0.5" << endl;
	plotFile << "set label 'c' at 0.94,-0.5" << endl;
	plotFile << endl;
	plotFile << "set arrow from 0.866,-0.5 to -0.866,-0.5 nohead" << endl;
	plotFile << "set arrow from 0.866,-0.5 to 0,1 nohead" << endl;
	plotFile << "set arrow from -0.866,-0.5 to 0,1 nohead" << endl;
	plotFile << endl;
	
	plotFile << "plot 'data.dat' using 1:2 with lines linecolor 1 ,\\" << endl;
	for (int i=1; i<Ntraj; i++)
		plotFile << "     'data.dat' using " << 2*i+1 << ":" << 2*i+2 << " with lines linecolor 1 ,\\" << endl;
	
	plotFile << endl;
	
	plotFile.close();
	
	//////////////////// Bash script to plot everything ////////////////////
	
	ofstream shellFile("portrait/plot.sh");
	
	shellFile << "#! /bin/bash" << endl;
	shellFile << "" << endl;
	shellFile << "gnuplot   plot.gp" << endl;
	shellFile << "" << endl;
	shellFile << "pdflatex  plot.tex" << endl;
	shellFile << "" << endl;
	shellFile << "rm   plot.aux   plot.log   plot.tex   plot-inc.eps   plot-inc-eps-converted-to.pdf" << endl;
	shellFile << "" << endl;
	shellFile << "mv plot.pdf portrait.pdf" << endl;
	shellFile << "" << endl;
	
	shellFile.close();
	
	// execute it
	sysresult = system("chmod +x portrait/plot.sh; cd portrait; ./plot.sh");
	
	
	return 0;
}

