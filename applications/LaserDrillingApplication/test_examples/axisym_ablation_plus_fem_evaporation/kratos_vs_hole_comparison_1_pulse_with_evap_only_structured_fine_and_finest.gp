set grid
#set size ratio -1
set key right center font "Times-Roman,17"
set xtics font "Times-Roman,17"
set ytics font "Times-Roman,17"
set xlabel "wide (um)" font "Times-Roman,17"
set ylabel "depth (um)" font "Times-Roman,17"
set term pngcairo size 1500,600
set output "laser_hole_structured_fine_and_finest.png"
plot 'hole_coords_l_s=0.0006_F_th=0.003_H_ev=400000.0_l_th=0.00025_alpha_ion=0.95_structured_finer.txt' u 2:(-$1) w l lw 2 lc rgb 'red' t 'Kratos structured fine', \
     'hole_coords_l_s=0.0006_F_th=0.003_H_ev=400000.0_l_th=0.00025_alpha_ion=0.95_structured_finer.txt' u (-$2):(-$1) w l lw 2 lc rgb 'red' t '', \
     'hole_coords_l_s=0.0006_F_th=0.003_H_ev=400000.0_l_th=0.00025_alpha_ion=0.95_structured_finest.txt' u 2:(-$1) w l lw 2 lc rgb 'blue' t 'Kratos structured finest', \
     'hole_coords_l_s=0.0006_F_th=0.003_H_ev=400000.0_l_th=0.00025_alpha_ion=0.95_structured_finest.txt' u (-$2):(-$1) w l lw 2 lc rgb 'blue' t '', \
     'experimental_hole_5_pulses.txt' u (($1 - 24.0)*1e-3):((-$2 + 1.0)*0.2e-3) w lp pt 6 lc rgb 'black' t 'experimental 1 pulse'
set output
