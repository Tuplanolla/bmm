# First fix the remainder line.
# There is an infinite family of possible lines,
# but if we demand that they are parametrized
# by step functions placed at the origin and
# have to be zero or one apart,
# only five candidates remain.

# set terminal pngcairo size 1600, 640 fontscale 0.8
# set output 'quot.png'
set xlabel ''
set ylabel ''
set xrange [-4.0 : 4.0]
set yrange [-4.0 : 5.0]
set samples 256
set key center top
set dummy x
set multiplot layout 2, 5
even(x) = x % 2 == 0
rem(x) = x - floor(x) - 1.0 / 2.0
quot(x) = x - rem(x)
stepn(x, y) = x < 0.0 ? y : 0.0
stepp(x, y) = x > 0.0 ? y : 0.0
ratio(n, d) = n == 0 || n == d ? sprintf('%d', n) : sprintf('%d / %d', n, d)
minus(n) = n < 0 ? '-' : ''
do for [i = 0 : 1] {
  s = i == 0 ? -1.0 : 1.0
  do for [j = 0 : 2] {
    do for [k = 0 : 2] {
      if (even(k - j)) {
        nt2 = j - 1.0
        pt2 = k - 1.0
        n = (1.0 / 2.0) * nt2
        p = (1.0 / 2.0) * pt2
        set title \
          sprintf('s = %d, n = %s, p = %s', s, ratio(nt2, 2), ratio(pt2, 2))
        plot s * (quot(x) - stepn(x, n) - stepp(x, p)) with lines title \
          sprintf('q = x / %sy', minus(s)), \
          rem(x) + stepn(x, n) + stepp(x, p) with lines title \
          sprintf('r = x %% %sy', minus(s)), \
          quot(x) + rem(x) with lines title \
          sprintf('x = q * %sy + r', minus(s))
      }
    }
  }
}
unset multiplot
