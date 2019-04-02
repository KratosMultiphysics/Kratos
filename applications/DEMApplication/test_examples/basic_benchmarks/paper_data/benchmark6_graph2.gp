set grid
set key right center
set style line 1 pt 8 lt -1 ps 3
set style line 2 pt 9 lt  3 ps 3
plot [0:8][-2:8] 'benchmark6_graph2.dat' index 0 w lp ls 1 t 'Al. alloy', \
		  '' index 1 w lp ls 2 t 'Nylon'
