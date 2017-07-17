# set terminal png size 1600, 800
# set output 'kdes.png'
set xrange [0.0 : 0.5]
set yrange [0 : *]
set xlabel 'r'
set ylabel 'g'
set dummy r
set style fill transparent solid 0.5
low(x) = (1.0 - 0.5) ** x
high(x) = (1.0 + 0.5) ** x
plot 'kde-amorphous.data' every 2 using 1 : (rand(0) * rand(0)) with points linetype 1 pointtype 0 notitle, \
  'kde-crystalline.data' every 2 using 1 : (rand(0) * rand(0)) with points linetype 3 pointtype 0 notitle, \
  'rubbish-amorphous.data' using 1 : ($2 * low($1)) : ($2 * high($1)) with filledcurves linetype 1 title 'Amorphous', \
  '' using 1 : 2 with lines linetype 1 notitle, \
  'rubbish-crystalline.data' using 1 : ($2 * low($1)) : ($2 * high($1)) with filledcurves linetype 3 title 'Cystalline', \
  '' using 1 : 2 with lines linetype 3 notitle, \
  1.0 with lines linetype 2 notitle
