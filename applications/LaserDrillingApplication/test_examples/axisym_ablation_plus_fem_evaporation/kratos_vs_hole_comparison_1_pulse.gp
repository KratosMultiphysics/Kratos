set grid
set key right center font "Sans,12"
set xtics font "Sans,14"
set ytics font "Sans,14"
set xlabel "wide (um)" font "Sans,17"
set ylabel "depth (um)" font "Sans,17"
set term pngcairo size 1500,600
set output "laser_hole.png"
plot 'list_of_decomposed_nodes_coords_no_evap.txt' u 2:(-$1) w lp pt 2 lw 2 lc rgb 'black' t 'Kratos no evap', \
     'list_of_decomposed_nodes_coords_no_evap.txt' u (-$2):(-$1) w lp pt 2 lw 2 lc rgb 'black' t '', \
     'list_of_decomposed_nodes_coords_with_evap.txt' u 2:(-$1) w lp pt 2 lw 2 lc rgb 'red' t 'Kratos with evap', \
     'list_of_decomposed_nodes_coords_with_evap.txt' u (-$2):(-$1) w lp pt 2 lw 2 lc rgb 'red' t '', \
     'experimental_hole_5_pulses.txt' u (($1 - 24.0)*1e-3):((-$2 + 1.0)*0.2e-3) w lp pt 6 lc rgb 'blue' t 'experimental 1 pulse'
set output
