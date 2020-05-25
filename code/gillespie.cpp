////////////////////////////////////////////////////////////////////////////
//                                                                        //
//    This program performs the simulation of cyclic predation networks   //
//    like the "rock paper scissors" game.                                //
//                                                                        //
//    It has been designed for an exam of a course on stochastic          //
//    processes -- PHYS-F446 (academic year 2019-2020) at ULB             //
//                                                                        //
//    Compile with $ g++ -O3 -fopenmp -o gillespie gillespie.cpp          //
//                                                                        //
//    Author: Cedric Schoonen                                             //
//                                                                        //
////////////////////////////////////////////////////////////////////////////


#include <random>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <omp.h>
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


int main()
{
	/////////////////////// Random number generator ////////////////////////
	
	random_device true_gen;
	int seed = true_gen();
	cout << "seed = " << seed << endl;
	cout << endl;
	
	// We create the rng engine inside of the parallelised loop
	// otherwise the same random number sequence is used by all threads
	
	uniform_real_distribution<double> dist01(0,1);
	
	////////////////////////////// Parameters //////////////////////////////
	
	double T = -1;      // simulation time (negative means unlimited)
	double rho = 0.02;  // a*b*c invariant of the deterministic dynamics
	                    // rho controls to number of individuals per species
	int N = 1000;       // number of individuals (from any species)
	int Ns = 3;         // number of species
	int M = 100;        // number of trajectories
	
	cout << endl;
	cout << "Parameters:" << endl;
	cout << endl;
	cout << "T =   " << T   << endl;
	cout << "rho = " << rho << endl;
	cout << "N = " << N << endl;
	cout << "Ns =  " << Ns  << endl;
	cout << "M =   " << M   << endl;
	cout << endl;
	
	////////////////////////////// Simulation //////////////////////////////
	
	// extinction time for all trajectories
	vector<double> t_extinction(M,0);
	
	#pragma omp parallel for 
	for (int m=0; m<M; m++)
	{
		// init system 
		
		double t = 0;
		vector<int> species(Ns,0);
		default_random_engine gen(seed+m);
		
		// compute number of individuals per species from rho
		// set b=c in the initial state
		
		double a = find_max_a(rho,Ns);
		double b = 1.0/(Ns-1) * (1-a);
		
		species[0] = int(a*N);
		for (int i=1; i<species.size(); i++) species[i] = int(b*N);
		
		
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
			if (extinction(species))
			{
				t_extinction[m] = t;
				break;
			}
			
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
			
			// Print final state of the system
			
			//cout << "State: ";
			//for (int i=0; i<Ns; i++) cout << species[i] << " ";
			//cout << " at t=" << t << endl;
		}
	}
	
	//////////////////////////// Report results ////////////////////////////
	
	// convert microscopic time (the one in the master equation)
	// to macroscopic time (the one in the Fokker-Planck equation)
	
	//for (int m=0; m<M; m++) 
	//	t_extinction[m] *= N;
	
	// mean extinction time
	
	double t_extinction_avg = 0;
	for (int m=0; m<M; m++) 
		t_extinction_avg += t_extinction[m]/M;
	
	// var and std of extinction time
	
	double t_extinction_var = 0;
	for (int m=0; m<M; m++) 
		t_extinction_var += pow(t_extinction[m]-t_extinction_avg,2)/(M-1);
	
	double t_extinction_std = sqrt(t_extinction_var);
	
	// Report
	
	cout << endl;
	cout << "Mean Extinction time (MS time) = " << t_extinction_avg;
	cout << " +/- " << t_extinction_std/sqrt(M) << endl;
	
	////////////////////////////////////////////////////////////////////////
	
	return 0;
}



