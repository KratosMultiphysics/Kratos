set grid
set key bottom center
set style line 1 pt 8 lt -1 ps 3
set style line 2 pt 9 lt  3 ps 3
plot [0:25][-1:.6] 'benchmark6_graph1.dat' index 0 w lp ls 1 t 'Al. alloy', \
		  '' index 1 w lp ls 2 t 'Nylon'
