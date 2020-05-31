#! /bin/bash

N=1000
Ns=3
M=10000

echo "#      rho  extc_time        err   N  Ns  M" > data_simul.dat

../gillespie_extc_time --script --rho 0.005 --N $N --Ns $Ns --M $M >> data_simul.dat
../gillespie_extc_time --script --rho 0.010 --N $N --Ns $Ns --M $M >> data_simul.dat
../gillespie_extc_time --script --rho 0.015 --N $N --Ns $Ns --M $M >> data_simul.dat
../gillespie_extc_time --script --rho 0.020 --N $N --Ns $Ns --M $M >> data_simul.dat
../gillespie_extc_time --script --rho 0.025 --N $N --Ns $Ns --M $M >> data_simul.dat
../gillespie_extc_time --script --rho 0.030 --N $N --Ns $Ns --M $M >> data_simul.dat
../gillespie_extc_time --script --rho 0.035 --N $N --Ns $Ns --M $M >> data_simul.dat

