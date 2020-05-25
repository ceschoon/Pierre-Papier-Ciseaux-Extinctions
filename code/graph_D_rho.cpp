////////////////////////////////////////////////////////////////////////////
//                                                                        //
//    This file is part of a project on extinctions in a cyclic           //
//    predation model like the "rock-paper-scissors" game.                //
//                                                                        //
//    Its purpuse is to make a graph of D(rho).                           //
//                                                                        //
//    Compile with $ g++ -O3 -o graph_D_rho graph_D_rho.cpp               //
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
	//////////// We compute D(rho) for several values of rho ///////////////
	
	double rho_min = 0;
	double rho_max = 1.0/27;
	double rho_step = 1e-4;
	
	int sysresult = system("mkdir -p D_rho");
	
	ofstream dataFile("D_rho/data.dat");
	dataFile << fixed << setprecision(8);
	dataFile << "# rho      D_rho  " << endl;
	
	for (double rho=rho_min; rho<rho_max; rho+=rho_step)
	{
		double D_rho = compute_D_rho(rho);
		dataFile << rho << " " << D_rho << endl;
	}
	
	/////////////////////////// Gnuplot script /////////////////////////////
	
	ofstream plotFile("D_rho/plot.gp");
	
	plotFile << "#----------------------------------------" << endl;
	plotFile << "set terminal epslatex standalone color size 9cm,7cm" << endl;
	plotFile << "set output 'plot.tex'" << endl;
	plotFile << endl;
	plotFile << "#set title '$\\bar D(\\rho)$'" << endl;
	plotFile << "set xlabel '$\\rho$'" << endl;
	plotFile << "set ylabel '$\\bar D(\\rho)$'" << endl;
	plotFile << "set style data lines" << endl;
	plotFile << "set key off" << endl;
	plotFile << endl;
	plotFile << "set xtics 0.01" << endl;
	plotFile << "set ytics 0.0005" << endl;
	plotFile << "set grid" << endl;
	plotFile << endl;
	plotFile << "plot 'data.dat'" << endl;
	plotFile << endl;
	
	plotFile.close();
	
	//////////////////// Bash script to plot everything ////////////////////
	
	ofstream shellFile("D_rho/plot.sh");
	
	shellFile << "#! /bin/bash" << endl;
	shellFile << "" << endl;
	shellFile << "gnuplot   plot.gp" << endl;
	shellFile << "" << endl;
	shellFile << "pdflatex  plot.tex" << endl;
	shellFile << "" << endl;
	shellFile << "rm   plot.aux   plot.log   plot.tex   plot-inc.eps   plot-inc-eps-converted-to.pdf" << endl;
	shellFile << "" << endl;
	shellFile << "mv plot.pdf D_rho.pdf" << endl;
	shellFile << "" << endl;
	
	shellFile.close();
	
	// execute it
	sysresult = system("chmod +x D_rho/plot.sh; cd D_rho; ./plot.sh");
	
	
	return 0;
}

