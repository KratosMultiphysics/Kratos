set grid
set key right center font "Sans,12"
set xtics font "Sans,14"
set ytics font "Sans,14"
set xlabel "wide (um)" font "Sans,17"
set ylabel "depth (um)" font "Sans,17"
set term pngcairo size 1500,600
set output "laser_hole_thermal_depths_".ARG1."_".ARG2.".png"
first_filename  = 'hole_coords_l_s=0.0006_F_th=0.003_H_ev=400000.0_l_th=0.00012_alpha_ion=0.95_'.ARG1.'_'.ARG2.'.txt'
second_filename = 'hole_coords_l_s=0.0006_F_th=0.003_H_ev=400000.0_l_th=0.00025_alpha_ion=0.95_'.ARG1.'_'.ARG2.'.txt'
third_filename  = 'hole_coords_l_s=0.0006_F_th=0.003_H_ev=400000.0_l_th=0.0005_alpha_ion=0.95_'.ARG1.'_'.ARG2.'.txt'
plot [][-0.0035:] 'experimental_hole_5_pulses.txt' u (($1 - 24.0)*1e-3):((-$2 + 1.0)*0.2e-3) w lp pt 6 lc rgb 'black' t 'experimental one pulse', \
     first_filename u 2:(-$1) w l lw 2 lc rgb 'green' t 'Kratos with l_{th} = 0.00012 mm', \
     first_filename u (-$2):(-$1) w l lw 2 lc rgb 'green' t '', \
     second_filename u 2:(-$1) w l lw 2 lc rgb 'red' t 'Kratos with l_{th} = 0.00025 mm', \
     second_filename u (-$2):(-$1) w l lw 2 lc rgb 'red' t '', \
     third_filename u 2:(-$1) w l lw 2 lc rgb 'blue' t 'Kratos with l_{th} = 0.0005 mm', \
     third_filename u (-$2):(-$1) w l lw 2 lc rgb 'blue' t ''
set output
