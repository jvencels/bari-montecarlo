load 'scripts/plot-common.gp'


set output srcdir.'convergence'.ext
set termoption dashed

set xlabel "Time [ns]"
set ylabel "Convergence"
set logscale y

set format x "%.0s"
set format y "%.0te%T"
set grid

set xtics mirror
set ytics mirror


plot srcdir.'convergence.dat' \
        u 1:2 w lines notitle

