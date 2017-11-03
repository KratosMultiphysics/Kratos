set grid
set key left bottom
set style line 1 pt 8 lt -1 ps 3
set style line 2 pt 9 lt 3 ps 3
plot [0:25][-10:10] 'benchmark7_dt_4.4614e-07_final_tangential_center_vel_list_data.dat' w lp lt 1 lw 1.5 ps 2 pt 5,\
'paper_data/benchmark7_graph1.dat' w lp ls 1 t 'Al. alloy',\
'paper_data/benchmark7_graph1.dat' w lp ls 2 t 'Copper'
