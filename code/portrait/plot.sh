#! /bin/bash

gnuplot   plot.gp
gnuplot   plot_trajectory.gp

pdflatex  plot.tex
pdflatex  plot_trajectory.tex

rm   *.aux   *.log   *.tex   *-inc.eps   *-inc-eps-converted-to.pdf

mv plot.pdf portrait.pdf
mv plot_trajectory.pdf portrait_trajectory.pdf

