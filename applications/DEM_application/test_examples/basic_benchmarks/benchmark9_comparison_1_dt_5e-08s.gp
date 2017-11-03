set grid
set key left bottom
set xlabel 'Coefficient of restitution'
set ylabel 'Damping ratio'
set style line 1 pt 8 lt -1 ps 3
set style line 2 pt 9 lt  3 ps 3
plot [0:1][0:1] 'benchmark9_dt_5e-08_restitution_numbers_vector_list_data.dat' w lp lt 1 lw 1.5 ps 2 pt 5,\
'paper_data/benchmark9_graph1.dat' w lp ls 1 t 'Al. oxide',\
'paper_data/benchmark9_graph1.dat' w lp ls 2 t 'Cast iron'
