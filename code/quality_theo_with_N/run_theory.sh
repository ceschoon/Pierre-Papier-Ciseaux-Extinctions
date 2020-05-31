#! /bin/bash

N=$1

cd ..

./graph_extc_time --interp --N $N
cp extc_time/data.dat quality_theo_with_N/data_interp.dat

./graph_extc_time --approx --N $N
cp extc_time/data.dat quality_theo_with_N/data_approx.dat


