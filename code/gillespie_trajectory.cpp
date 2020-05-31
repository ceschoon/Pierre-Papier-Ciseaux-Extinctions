////////////////////////////////////////////////////////////////////////////
//                                                                        //
//    This program performs the simulation of cyclic predation networks   //
//    like the "rock paper scissors" game.                                //
//                                                                        //
//    It has been designed for an exam of a course on stochastic          //
//    processes -- PHYS-F446 (academic year 2019-2020) at ULB             //
//                                                                        //
//    Compile with $ g++ -O3 -fopenmp -o gillespie_trajectory \           //
//                   gillespie_trajectory.cpp                             //
//                                                                        //
//    Author: Cedric Schoonen                                             //
//                                                                        //
////////////////////////////////////////////////////////////////////////////


#include <random>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>
#include <omp.h>
#include <stdexcept>
#include "Integration_D_rho.h"  // find_max_a

using namespace std;

double const PI = 3.14159265358979;



bool extinction(vector<int> species)
{
	for (int i=0; i<species.size(); i++)
		if (species[i] == 0)
			return true;
	
	return false;
}


double predation_rate(int predator, int pray, int Ns)
{
	//// no intra-species competition
	if (predator==pray) return 0;
	
	//// cyclic predation network
	if (predator+1==pray) return 1;
	else if (predator+1==Ns && pray==0) return 1;
	
	//// custom rules
	//if (predator==0 && pray==1) return 1;
	//if (predator==0 && pray==2) return 1;
	//if (predator==1 && pray==2) return 1;
	
	//// random predation network
	//uniform_real_distribution<double> dist01(0,1);
	//default_random_engine gen(predator*Ns+pray);
	//double ran = dist01(gen);
	//return dist01(gen);
	
	//// default if 0 (no reaction)
	return 0;
}


int main(int argc, char **argv)
{
	/////////////////////// Random number generator ////////////////////////
	
	random_device true_gen;
	int seed = true_gen();
	
	// We create the rng engine inside of the parallelised loop
	// otherwise the same random number sequence is used by all threads
	
	uniform_real_distribution<double> dist01(0,1);
	
	////////////////////////////// Parameters //////////////////////////////
	
	bool script_mode = false;
	bool cout_enabled = true;
	
	double T = -1;        // simulation time (negative means unlimited)
	double rho = 1.0/27;  // a*b*c invariant of the deterministic dynamics
	                      // rho controls to number of individuals per species
	int N = 1000;         // number of individuals (from any species)
	const int Ns = 3;     // number of species
	int SamplingRate = 1; // rate at which data is recorded
	
	ofstream dataFile("trajectory.dat");
	
	////////////////////////////// Options /////////////////////////////////
	
	for (int i=1; i<argc; i++)
	{
		if ((string(argv[i])=="--rho")) 
		{
			if (argc > i+1) rho = stod(string(argv[i+1]));
			else 
				throw runtime_error("You must specify the value of rho, e.g. --rho 0.02");
		}
		
		if ((string(argv[i])=="--N")) 
		{
			if (argc > i+1) N = stoi(string(argv[i+1]));
			else 
				throw runtime_error("You must specify the value of N, e.g. --N 100");
		}
		
		if ((string(argv[i])=="--SamplingRate")) 
		{
			if (argc > i+1) SamplingRate = stoi(string(argv[i+1]));
			else 
				throw runtime_error("You must specify the value of SamplingRate, e.g. --SamplingRate 100");
		}
	}
	
	////////////////////////////// Simulation //////////////////////////////
	
	// Report parameters
	
	if (cout_enabled)
	{
		cout << endl;
		cout << "seed = " << seed << endl;
		cout << endl;
		cout << "Parameters:" << endl;
		cout << endl;
		cout << "T =   " << T   << endl;
		cout << "rho = " << rho << endl;
		cout << "N = " << N << endl;
		cout << "Ns =  " << Ns  << endl;
		cout << "SamplingRate =  " << SamplingRate  << endl;
		cout << endl;
	}
	
	// init system 
	
	double t = 0;
	int reactionCounter = 0;
	vector<int> species(Ns,0);
	default_random_engine gen(seed);
	
	// compute number of individuals per species from rho
	// set b=c in the initial state
	
	double a = find_max_a(rho,Ns);
	double b = 1.0/(Ns-1) * (1-a);
	double c = 1.0/(Ns-1) * (1-a);
	
	double x = (c-b) * sqrt(3)/2;
	double y = a-b/2-c/2;
	
	species[0] = int(a*N);
	for (int i=1; i<species.size(); i++) species[i] = int(b*N);
	
	// init file
	
	dataFile << "#      t       A       B       C       x       y" 
	         << endl << endl;
	dataFile << fixed << setprecision(8);
	dataFile << t << " " << species[0] << " " << species[1] << " "
	         << species[2] << " " << x << " " << y << endl;
	
	// time loop
	
	while (t<T || T<0)
	{
		// compute total reaction rate
		
		double rate_tot = 0;
		
		for (int i=0; i<Ns; i++) 
		{
			for (int j=0; j<Ns; j++) 
			{
				double rate_ij = predation_rate(i,j,Ns)
				                 *species[i]*species[j];
				rate_tot += rate_ij;
				
				//cout << "Rate for i=" << i << " j=" << j << " is ";
				//cout << rate_ij << endl;
			}
		}
		
		// detect extinction
		
		// if (rate_tot==0)
		if (extinction(species)) break;
		
		// advance time
		
		double ran1 = dist01(gen);
		t += -1.0/rate_tot*log(ran1);
		
		// select reaction
		
		double ran2 = dist01(gen);
		double proba_cumul = 0;
		int predator = -1;
		int pray = -1;
		
		for (int i=0; i<Ns; i++) 
		{
			for (int j=0; j<Ns; j++) 
			{
				double rate_ij = predation_rate(i,j,Ns)
				                 *species[i]*species[j];
				proba_cumul += rate_ij/rate_tot;
				
				if (ran2 < proba_cumul)
				{
					predator = i;
					pray = j;
					break;
				}
			}
			
			if (ran2 < proba_cumul) break;
		}
		
		// perform reaction
		
		species[predator]++;
		species[pray]--;
		reactionCounter ++;
		
		// print in file
		
		if (reactionCounter%SamplingRate == 0)
		{
			// convert to mesoscopic variables
			
			double a = double(species[0])/N;
			double b = double(species[1])/N;
			double c = double(species[2])/N;
			
			// transform (a,b,c) to (x,y) coordinates in the (1,1,1) plane
			// and save the position in data file
			
			double x = (c-b) * sqrt(3)/2;
			double y = a-b/2-c/2;
			
			// record
			
			dataFile << t << " " << species[0] << " " << species[1] << " "
			         << species[2] << " " << x << " " << y << endl;
		}
	}
	
	////////////////////////////////////////////////////////////////////////
	
	return 0;
}



