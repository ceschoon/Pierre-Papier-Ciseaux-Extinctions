#! /bin/bash

N=$1
Ns=3
M=10000

dataFile=data_simul.dat

echo "#      rho  extc_time        err   N  Ns  M" > $dataFile

../gillespie_extc_time --script --rho 0.005 --N $N --Ns $Ns --M $M >> $dataFile
../gillespie_extc_time --script --rho 0.010 --N $N --Ns $Ns --M $M >> $dataFile
../gillespie_extc_time --script --rho 0.015 --N $N --Ns $Ns --M $M >> $dataFile
../gillespie_extc_time --script --rho 0.020 --N $N --Ns $Ns --M $M >> $dataFile
../gillespie_extc_time --script --rho 0.025 --N $N --Ns $Ns --M $M >> $dataFile
../gillespie_extc_time --script --rho 0.030 --N $N --Ns $Ns --M $M >> $dataFile
../gillespie_extc_time --script --rho 0.035 --N $N --Ns $Ns --M $M >> $dataFile

