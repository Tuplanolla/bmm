#! /bin/sh

f=ballmoi.trace && \
t="$(mktemp)" && \
trap 'rm -f "$t"' EXIT && \
n=$(echo 2 ^ 32 | bc) && \
for d in $(seq 0 32)
do for p in 0 1
do nice -n 10 ./ballmoi $p $d $n >> "$t" &
done
done && \
wait && \
sort -V "$t" > "$f"
