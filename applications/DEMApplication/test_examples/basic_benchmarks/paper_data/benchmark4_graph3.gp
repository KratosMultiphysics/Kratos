set grid
set key right center
set style line 1 pt 8 lt -1 ps 3
set style line 2 pt 9 lt  3 ps 3
#set style line 2 pt 9 lt -1 ps 3
plot [0:90][-30:90] 'benchmark4_graph3.dat' index 0 w lp ls 1 t 'Al. oxide', \
		  '' index 1 w lp ls 2 t 'Al. alloy', \
		  '' index 2 w p pt 7 ps 2 lt -1 t 'Experiment'
