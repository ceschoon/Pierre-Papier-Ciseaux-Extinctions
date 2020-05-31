#! /bin/bash

N=$1

gnuplot   plot.gp

pdflatex  plot.tex

rm   plot.aux   plot.log   plot.tex   plot-inc.eps   plot-inc-eps-converted-to.pdf

mv plot.pdf extc_time_comparison_N$N.pdf

