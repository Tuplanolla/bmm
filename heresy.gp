set terminal epslatex
set output 'figures/heresy.tex'
# stats 'heresy.data' using 4 name 'stats_periodic'
# periodic = stats_periodic_mean
periodic = 1
# stats 'heresy.data' using 5 name 'stats_L'
# L = stats_L_mean
L = 1.0
set xlabel '$x \: [\si \meter]$'
set ylabel '$y \: [\si \meter]$'
set size ratio -1
unset key
nearbyint(x) = floor(x + 0.5)
swrap(x, b) = x - b * nearbyint(x / b)
fx(i, x, y) = i == 0 ? \
  (x2 = NaN, y2 = NaN, x1 = x, y1 = y, x2) : \
  (x2 = x1, y2 = y1, x1 = x, y1 = y, x2)
fy(i, x, y) = y2
if (periodic) {
  dfx(i, x, y) = swrap(x1 - x2, L)
  dfy(i, x, y) = swrap(y1 - y2, L)
  set xrange [-L / 2.0 : L + L / 2.0]
  set yrange [-L / 2.0 : L + L / 2.0]
  set object 1 rectangle from 0.0, 0.0 to L, L fillstyle empty
  plot for [dx = -L : L : L] for [dy = -L : L : L] \
    'heresy.data' using ($2 + dx) : ($3 + dy) : ($4 * 2.0) \
    with ellipses linetype 1, \
    for [dx = -L : L : L] for [dy = -L : L : L] \
    'more-heresy.data' \
    using (fx($1, $2 + dx, $3 + dy)) : (fy($1, $2 + dx, $3 + dy)) : \
    (dfx($1, $2 + dx, $3 + dy)) : (dfy($1, $2 + dx, $3 + dy)) \
    with vectors linetype 2 filled
} else {
  dfx(i, x, y) = x1 - x2
  dfy(i, x, y) = y1 - y2
  set xrange [-L / 2.0 : L / 2.0]
  set yrange [-L / 2.0 : L / 2.0]
  plot 'heresy.data' using 2 : 3 : ($4 * 2.0) \
    with ellipses linetype 1, \
    'more-heresy.data' \
    using (fx($1, $2, $3)) : (fy($1, $2, $3)) : \
    (dfx($1, $2, $3)) : (dfy($1, $2, $3)) \
    with vectors linetype 2 filled
}
