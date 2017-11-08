set grid
set key left bottom
set xlabel 'Normalized incident angle'
set ylabel 'Normalized final angular velocity'
set style line 1 pt 8 lt -1 ps 3
set style line 2 pt 9 lt 3 ps 3
plot [0:14][-6:0] 'benchmark11_dt_1e-07_normalized_rebound_angular_velocity_list_data.dat' w lp lt 1 lw 1.5 ps 2 pt 5,\
'paper_data/bench_10_norm_reb_ang_vel_e_090.dat' index 0 w lp ls 1 t 'Paper data'
