load 'scripts/plot-common.gp'

#tiny=1
minimum=1E-5

set output srcdir.'eedf-evolution'.ext
set logscale y
set yrange [minimum:1]

#set key center top

#ts = system("awk '{ print $1 }' ".srcdir."eedf.dat | uniq | awk '{for (i=1; i<20; i++) {if (NR==i*i+1) {print $1}}}'")

ts = system("awk '{ print $1 }' ".srcdir."eedf.dat | uniq | awk 'NR==2; END{print}'")

thistime(t) = '<awk ''/^ '.t.'/'' '.srcdir.'eedf.dat'

set xrange [-10:1010]
set yrange [1e-9:1]

set ylabel 'EEDF'

set xlabel  'Energy [eV]'


plot for [t in ts] thistime(t) u 2:3 w linespoints notitle

