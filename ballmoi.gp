set dummy d
set xlabel 'd'
set ylabel 'j'
set xrange [0 < * :]
set yrange [0 < * :]
plot 'ballmoi.data' using 2 : 4 : 5 \
  with yerrorbars linestyle 1 notitle, \
  for [i = -1 : 1 : 2] '' every ::1 using 2 : ($4 + i * $5) : (1 / (1 + $5)) \
  with lines linestyle 2 smooth acsplines notitle, \
  for [p = 0 : 1] (d - (1 - p)) / (d + 2) \
  with lines linestyle 3 notitle
