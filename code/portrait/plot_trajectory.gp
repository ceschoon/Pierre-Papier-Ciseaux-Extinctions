#----------------------------------------
set terminal epslatex standalone color size 5cm,5cm
set output 'plot_trajectory.tex'

#set title ''
#set xlabel '$\frac{\sqrt{3}}{2} (c-b)$'
#set ylabel '$a-c/2-b/2$'
set key off

unset xtics
unset ytics
unset border

set xrange [-0.85:0.80]
set yrange [-0.70:1.15]

set label 'a' at -0.05,1.1
set label 'b' at -1.03,-0.5
set label 'c' at 0.94,-0.5

set arrow from 0.866,-0.5 to -0.866,-0.5 nohead
set arrow from 0.866,-0.5 to 0,1 nohead
set arrow from -0.866,-0.5 to 0,1 nohead

plot 'portrait.dat' using 1:2 with lines linecolor 1 ,\
     'portrait.dat' using 3:4 with lines linecolor 1 ,\
     'portrait.dat' using 5:6 with lines linecolor 1 ,\
     'portrait.dat' using 7:8 with lines linecolor 1 ,\
     'portrait.dat' using 9:10 with lines linecolor 1 ,\
     'portrait.dat' using 11:12 with lines linecolor 1 ,\
     'trajectory.dat' using 5:6 with linespoints pointsize 0.3 linecolor 2

