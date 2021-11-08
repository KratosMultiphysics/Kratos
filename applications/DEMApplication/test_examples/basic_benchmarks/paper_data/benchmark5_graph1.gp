set grid
set key right center
set style line 1 pt 8 lt -1 ps 3
set style line 2 pt 9 lt  3 ps 3
set style line 3 pt 7 lt -1 ps 2
plot [0:14][-4:6] 'benchmark5_graph1.dat' index 0 w lp ls 1 t 'Steel', \
		  '' index 1 w lp ls 2 t 'Polyethylene', \
		  '' index 2 w lp ls 3 t 'FEM'
