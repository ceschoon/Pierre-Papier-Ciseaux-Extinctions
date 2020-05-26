#! /bin/bash

cd ..

./graph_extc_time --interp
cp extc_time/data.dat comparison_extc/data_interp.dat

./graph_extc_time --approx
cp extc_time/data.dat comparison_extc/data_approx.dat


