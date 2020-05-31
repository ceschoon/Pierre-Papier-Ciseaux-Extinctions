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
	
	ofstream dataFile("portrait.dat");
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
	
	return 0;
}

