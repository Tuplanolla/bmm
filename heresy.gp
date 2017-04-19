# set terminal epslatex
# set output 'proj.tex'
# stats 'heresy.data' using 4 name 'stats_periodic'
# periodic = stats_periodic_mean
periodic = 1
# stats 'heresy.data' using 5 name 'stats_L'
# L = stats_L_mean
L = 1.0
set xlabel '$x$'
set ylabel '$y$'
if (periodic) {
  set xrange [-L / 2.0 : L + L / 2.0]
  set yrange [-L / 2.0 : L + L / 2.0]
  set object 1 rectangle from 0.0, 0.0 to L, L fillstyle empty
  plot for [dx = -L : L : L] for [dy = -L : L : L] \
    'heresy.data' using ($1 + dx) : ($2 + dy) : ($3 * 2.0) \
    with ellipses linetype 1 notitle
} else {
  set xrange [0.0 : L]
  set yrange [0.0 : L]
  plot 'heresy.data' using 1 : 2 : ($3 * 2.0) \
    with ellipses linetype 1 notitle
}
