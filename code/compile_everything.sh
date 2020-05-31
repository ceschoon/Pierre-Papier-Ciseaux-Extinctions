#! /bin/bash

g++ -O3 -fopenmp -o gillespie_extc_time gillespie_extc_time.cpp
g++ -O3 -fopenmp -o gillespie_trajectory gillespie_trajectory.cpp
g++ -O3 -o graph_D_rho graph_D_rho.cpp
g++ -O3 -o graph_extc_time graph_extc_time.cpp
g++ -O3 -o phase_portrait phase_portrait.cpp

