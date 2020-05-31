#! /bin/bash

N=10
./run_simulations.sh $N
./run_theory.sh $N
./plot.sh $N

N=20
./run_simulations.sh $N
./run_theory.sh $N
./plot.sh $N

N=50
./run_simulations.sh $N
./run_theory.sh $N
./plot.sh $N

N=100
./run_simulations.sh $N
./run_theory.sh $N
./plot.sh $N

N=200
./run_simulations.sh $N
./run_theory.sh $N
./plot.sh $N

N=500
./run_simulations.sh $N
./run_theory.sh $N
./plot.sh $N

N=1000
./run_simulations.sh $N
./run_theory.sh $N
./plot.sh $N