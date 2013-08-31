load 'scripts/plot-common.gp'


set output srcdir.'totals'.ext
set termoption dashed

set xlabel "Time [ns]"
set ylabel "Average energy [eV]" textcolor rgb "red"
set y2label "Number of electrons" rotate by -90 textcolor rgb "green"
set yrange [0:]
set y2range [0:]

set format x "%.0s"
set format y "%.0te%T"
set format y2 "%.0te%T"
set grid

set xtics nomirror
set ytics nomirror tc rgb "red" # 2e-9
set y2tics nomirror tc rgb "green"

set key outside above 

#ts = system("awk '{ print $1 }' ".srcdir."norm-e-distr.dat | uniq | awk '{for (i=1; i<20; i++) {if (NR==i*i+1) {print $1}}}'")

#set for [t in ts] arrow from t,graph(0,0) to t,graph(1,1) nohead

plot srcdir.'totals.dat' \
        u 1:2 w lines lt 1 lc rgb "red" notitle, \
    ''  u 1:3 axes x1y2 w lines lt 1 lc rgb "green" notitle, \
    ''  u 1:($4*20) axes x1y2 w lines lt 1 lc rgb "blue" notitle, \
    ''  u 1:5 axes x1y2 w lines lt 1 lc rgb "black" notitle
