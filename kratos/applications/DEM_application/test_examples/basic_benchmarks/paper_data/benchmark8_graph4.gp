set grid
set key right center
set style line 1 pt 8 lt -1 ps 2
set style line 2 pt 6 lt -1 ps 2
set style line 3 pt 4 lt -1 ps 2
plot [0:8][-2:10] 'benchmark8_graph4.dat' index 0 w lp ls 1 t 'Two spheres impact (case A)', \
		  '' index 1 w lp ls 2 t 'Two spheres impact (case B)', \
		  '' index 2 w lp ls 3 t 'Sphere vs rigid surface (Test 6)'
