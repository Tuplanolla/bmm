# First fix the remainder line.
# There is an infinite family of possible lines,
# but if we demand that they are parametrized
# by integer step functions placed at the origin,
# only four candidates remain.
# Eight more appear if the periodicity is doubled.

# set terminal pngcairo size 1200, 1600 fontscale 0.8
# set output 'quot.png'
set xlabel ''
set ylabel ''
set xrange [-4.0 : 4.0]
set yrange [-4.0 : 5.0]
set samples 256
set grid
set key center top
set dummy x
set multiplot layout 6, 4
rem(x) = x - floor(x) - 1.0 / 2.0
quot(x) = x - rem(x)
stepn(x, y) = x < 0.0 ? y : 0.0
stepp(x, y) = x > 0.0 ? y : 0.0
ratio(n, d) = n == 0 || n == d ? sprintf('%d', n) : sprintf('%d / %d', n, d)
plusorminus(x) = x < 0.0 ? sprintf('%d', x) : sprintf('+%d', x)
minus(x) = x < 0.0 ? '-' : ''
do for [i = 0 : 1] {
  s = 1.0 - i * 2.0
  do for [j = 0 : 1] {
    do for [k = 0 : 1] {
      nt2 = j * 2.0 - 1.0
      pt2 = k * 2.0 - 1.0
      n = (1.0 / 2.0) * nt2
      p = (1.0 / 2.0) * pt2
      set title sprintf('s = %s, n = %s, p = %s', \
          plusorminus(s), ratio(nt2, 2), ratio(pt2, 2))
      plot s * (quot(x) - stepn(x, n) - stepp(x, p)) with lines title \
        sprintf('q = x / %sy', minus(s)), \
        rem(x) + stepn(x, n) + stepp(x, p) with lines title \
        sprintf('r = x %% %sy', minus(s)), \
        quot(x) + rem(x) with lines title \
        sprintf('x = q * %sy + r', minus(s))
    }
  }
}
rem2(x) = x - floor(x / 2.0) * 2.0 - 1.0
quot2(x) = x - rem2(x)
plusandminus(x) = x == 0.0 ? sprintf('%d', x) : sprintf('+-%d', x)
do for [i = 0 : 1] {
  z = i
  do for [j = 0 : 1] {
    s = 1.0 - j * 2.0
    do for [k = 0 : 1] {
      do for [l = 0 : 1] {
        b = k
        a = l
        set title sprintf('s = %s, b = %s, a = %s, z = %s', \
          plusorminus(s), plusandminus(b), plusandminus(a), plusandminus(z))
        plot s * (quot2(x + stepn(x, b) + stepp(x, a) - z) - \
          stepn(x, b) - stepp(x, a) + z) with lines title \
        sprintf('q = x / %sy', minus(s)), \
        rem2(x + stepn(x, b) + stepp(x, a) - z) with lines title \
        sprintf('r = x %% %sy', minus(s)), \
        quot2(x) + rem2(x) with lines title \
        sprintf('x = q * %sy + r', minus(s))
      }
    }
  }
}
unset multiplot
