#----------------------------------------
set terminal epslatex standalone color size 9cm,7cm
set output 'plot.tex'

#set title '$T(\rho)$'
set xlabel '$\rho$'
set ylabel '$T(\rho)$'
set key bottom right

set xtics 0.01
set ytics 0.1
set grid

plot 'data_interp.dat' with lines  linecolor 1 title 'interp' ,\
     'data_approx.dat' with lines  linecolor 2 title 'approx' ,\
     'data_simul.dat'  using 1:2:($2-$3):($2+$3) with yerrorbars pointtype 1 linecolor 8 title 'simul' ,\

