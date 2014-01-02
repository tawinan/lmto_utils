 set term postscript enh color font 'Times-Roman,24'
 set output 'bnds.ps'
 set multiplot
 set noxzeroaxis
 set tics out
 set noxtics
 set nokey
 set xtics ( ' L '   0.0000 ,' {/Symbol G} '   1.0000 ,' X '   2.1547 ,' W '   2.7321 ,' {/Symbol G} '   4.0230 )
 set lmargin screen 0.12
 set rmargin screen 0.85
 set bmargin screen 0.10
 set tmargin screen 0.90
 set pm3d map
 set hidden3d
 set palette defined ( 0 "white", .2 "#70FEBD"  ,40 "blue",70 "red", 100 "white")
 set title 'SmN'
 set yrange [ -5.0000 :   5.0000] 
 set ylabel 'Energy (eV)'
 splot 'WCOLOR.DAT'
 unset xtics
 unset ylabel
 unset ytics
 unset title
 set style line 1 lt 1 lw 2 lc rgb 'red'
 set style line 2 lt 2 lw 2 lc rgb 'green'
 set style line 3 lt 1 lw 1 lc rgb 'black'
 plot for[i=2:13] 'BNDS.DAT' using 1:i with line ls 1, \
		for[i=2:13] 'BNDS2.DAT' using 1:i with line ls 2,\
		'FERMI.DAT' using 1:2 with line ls 3
 unset multiplot
