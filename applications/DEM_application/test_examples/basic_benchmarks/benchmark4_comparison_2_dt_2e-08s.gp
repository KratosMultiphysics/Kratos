set grid
set key left bottom
set xlabel 'Incident angle (deg)'
set ylabel 'Final angular velocity (rad/s)'
set style line 1 pt 8 lt -1 ps 3
set style line 2 pt 9 lt  3 ps 3
plot [0:90][-750:0] 'benchmark4_dt_2e-08_final_angular_vel_list_data.dat' w lp lt 1 lw 1.5 ps 2 pt 5,\
'paper_data/benchmark4_graph2.dat' index 0 w lp ls 1 t 'Al. oxide',\
'paper_data/benchmark4_graph2.dat' index 1 w lp ls 2 t 'Al. alloy',\
'paper_data/benchmark4_graph2.dat' index 2 w p pt 7 ps 2 lt -1 t 'Experimental'
