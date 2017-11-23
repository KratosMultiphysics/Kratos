set grid
set key left bottom
set xlabel 'Normalized incident angle'
set ylabel 'Normalized rebound tangential surface velocity'
set style line 1 pt 8 lt -1 ps 3
set style line 2 pt 9 lt 3 ps 3
plot [0:10][-2:3] 'benchmark10_dt_2e-08_normalized_rebound_tangential_surface_vel_list_data.dat' w lp lt 1 lw 1.5 ps 2 pt 5,\
'paper_data/bench_10_norm_reb_tang_e_090.dat' index 1 w lp ls 1 t 'Paper data'
