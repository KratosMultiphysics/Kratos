set grid
set key left bottom
set xlabel 'Normalized incident angle'
set ylabel 'Tangential coefficient of restitution'
set style line 1 pt 8 lt -1 ps 3
set style line 2 pt 9 lt 3 ps 3
plot [0:10][0.5:1.0] 'benchmark11_dt_1e-07_tangential_coefficient_of_restitution_list_data.dat' w lp lt 1 lw 1.5 ps 2 pt 5,\
'paper_data/bench_10_tang_coeff_rest_e_090.dat' index 0 w lp ls 1 t 'Paper data'
