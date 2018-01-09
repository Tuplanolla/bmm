set xlabel 't'
set ylabel ''
plot 'est.data' using 1 : ($2 + $3 - $4) with linespoints title 'E', \
  '' using 1 : 2 with linespoints title 'E+', \
  '' using 1 : 3 with linespoints title 'E*', \
  '' using 1 : 4 with linespoints title 'E-', \
  '' using 1 : 5 with linespoints title '\mu', \
  '' using 1 : 6 with linespoints title '\mu_b', \
  '' using 1 : 7 with linespoints title 'vx', \
  '' using 1 : 8 with linespoints title 'vy', \
  '' using 1 : 9 with linespoints title 'Fx', \
  '' using 1 : 10 with linespoints title 'Fy'
