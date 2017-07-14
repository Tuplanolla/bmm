set xrange [0.0 : 0.5]
stats 'kde.data' using 2 name 'stats'
set table 'kde-table.data'
set samples 256
plot 'kde.data' using 1 : ($2 / stats_sum) smooth kdensity bandwidth 3.0e-3
unset table
set dummy r
plot 'kde.data' every 2 using 1 : (rand(0) * rand(0)) with points linetype 1 pointtype 0 notitle, \
  'kde-table.data' using 1 : ($2 / (2*pi * $1)) with lines linetype 1 notitle, \
  1.0 with lines linetype 3 notitle
