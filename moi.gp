set xlabel 'd'
set ylabel 'j'
set dummy d
plot 'moi.data' using 1 : 2 : 3 \
  with yerrorbars linestyle 1 notitle, \
  for [i = -1 : 1 : 2] '' every ::1 using 1 : ($2 + i * $3) : (1 / (1 + $3)) \
  with lines linestyle 2 smooth acsplines notitle, \
  (d - 1) / (d + 2) with lines linestyle 3 notitle
