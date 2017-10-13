set grid
set key right center
set style line 1 pt 8 lt -1 ps 3
set style line 2 pt 9 lt  3 ps 3
#set style line 2 pt 9 lt -1 ps 3
plot [0:20][-6:0] 'benchmark5_graph2.dat' index 0 w lp ls 1 t 'Steel', \
		  '' index 1 w lp ls 2 t 'Polyethylene', \
		  '' index 2 w p pt 7 ps 2 lt -1 t 'FEM'
