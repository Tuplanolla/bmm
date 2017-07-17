set xrange [0.0 : 0.5]
stats 'kde-amorphous.data' using (1.0 / $2) name 'stats'
set table 'kde-amorphous-table.data'
set samples 256
plot 'kde-amorphous.data' using 1 : (1.0 / ($2 * stats_sum)) smooth kdensity bandwidth 3.0e-3
unset table
stats 'kde-crystalline.data' using 2 name 'stats'
set table 'kde-crystalline-table.data'
set samples 256
plot 'kde-crystalline.data' using 1 : (1.0 / ($2 * stats_sum)) smooth kdensity bandwidth 3.0e-3
unset table
set dummy r
plot 'kde-amorphous.data' every 2 using 1 : (rand(0) * rand(0)) with points linetype 1 pointtype 0 notitle, \
  'kde-crystalline.data' every 2 using 1 : (rand(0) * rand(0)) with points linetype 2 pointtype 0 notitle, \
  'kde-amorphous-table.data' using 1 : ($2 * (3.0 / 2.0) / (2*pi * $1)) with lines linetype 1 title 'Amorphous', \
  'kde-crystalline-table.data' using 1 : ($2 * (2.0 / 3.0) / (2*pi * $1)) with lines linetype 2 title 'Cystalline', \
  1.0 with lines linetype 3 notitle
