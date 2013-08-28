load 'scripts/plot-common.gp'

tiny=1
minimum=1E-5

set output srcdir.'norm-e-distr'.ext
set logscale y
set yrange [minimum:1]

#set key center top
#set xlabel 'Energy [eV]'


#ts = system("awk '{ print $1 }' ".srcdir."norm-e-distr.dat | uniq | awk '{for (i=1; i<20; i++) {if (NR==i*i+1) {print $1}}}'")

ts = system("awk '{ print $1 }' ".srcdir."norm-e-distr.dat | uniq | awk 'NR==2; END{print}'")

thistime(t) = '<awk ''/^ '.t.'/'' '.srcdir.'norm-e-distr.dat'

set xrange [-10:1010]
set yrange [1e-7:1]

set ylabel 'Probability [%]'

set xlabel  'Energy [eV]'


plot for [t in ts] thistime(t) u 2:3 w linespoints notitle



#set multiplot

#set border 1+2+4
#set lmargin at screen 0.20
#set rmargin at screen 0.55
#set bmargin at screen  0.15
#set tmargin at screen 0.95
#set xrange [-10:100]
#set xtics ("0" 0, "50" 50)
#set ytics nomirror
#set arrow 1 from 100-tiny, minimum/2 to 100+tiny, minimum*2 nohead
#set arrow 2 from 100-tiny, 1.0/2 to 100+tiny, 1*2 nohead
#set ylabel 'Probability [%]'

#plot for [t in ts] thistime(t) u 2:3 w linespoints title sprintf("%11.2e  [s]",t+0)

#unset ylabel
#unset ytics
#set border 1+4+8
#set key right
#set lmargin at screen 0.57
#set rmargin at screen 0.95
#set bmargin at screen 0.15
#set tmargin at screen 0.95
#set label 1 'Energy [eV]' at screen 0.5, 0.05 centre
#set xrange [900:1010]
#set xtics ("950" 950, "1000" 1000)
#set arrow 1 from 900-tiny, minimum/2 to 900+tiny, minimum*2 nohead
#set arrow 2 from 900-tiny, 1.0/2 to 900+tiny, 1*2 nohead

#plot for [t in ts] thistime(t) u 2:3 w linespoints notitle
#unset multiplot
