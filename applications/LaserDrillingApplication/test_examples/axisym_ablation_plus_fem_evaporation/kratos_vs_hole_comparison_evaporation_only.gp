set grid
set key right center font "Sans,12"
set xtics font "Sans,14"
set ytics font "Sans,14"
set xlabel "wide (um)" font "Sans,17"
set ylabel "depth (um)" font "Sans,17"
set term pngcairo size 1500,600
set output "laser_hole_ablation_only.png"
filename = 'hole_coords_l_s=3.089478307077969e-05_F_th=0.003_H_ev=400000.0_l_th=0.00025_alpha_ion=0.915_structured_fine.txt'
plot [][-0.0035:] 'experimental_hole_5_pulses.txt' u (($1 - 24.0)*1e-3):((-$2 + 1.0)*0.2e-3) w lp pt 6 lc rgb 'black' t 'experimental one pulse', \
     filename u 2:(-$1) w l lw 2 lc rgb 'blue' t 'Kratos', \
     filename u (-$2):(-$1) w l lw 2 lc rgb 'blue' t ''
set output
