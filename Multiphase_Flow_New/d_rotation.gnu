#set term postscript eps enhanced "Helvetica" 10
#set title "rarefraction wave a two different angles"

set key center top

set pointsize 0.5

set xlabel "Position"

set yrange [700:820]
set ylabel "Density"

#set size 0.5,0.5
#set key 1.5,32
#set logscale x
#set xrange [690:840]

#set xtics ( 3, 4, 5 )
#set format y "%e"
#set logscale y
#set ytics ( 1, 2, 4, 8, 16, 32 )
#set grid

plot "./d_100x200_Lax-Friedrich_2d_split1_IC0_div5till5_8Steps_d1" title "0 degrees", \
     "./d_100x200_Lax-Friedrich_2d_split1_IC1_div5till5_8Steps_d1" title "60 degrees",
pause -1
set term post enh 16
set output 'density_rotation.eps
replot
#quit
