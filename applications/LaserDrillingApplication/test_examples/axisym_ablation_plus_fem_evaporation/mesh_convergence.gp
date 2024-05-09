set grid lw 2
set key right center font "Times-Roman,17"
set xtics font "Times-Roman,15"
set ytics font "Times-Roman,15"
set xlabel "Number of elements" font "Times-Roman,17"
set ylabel "Temperature (K)" font "Times-Roman,17"
set term pngcairo size 900,600
set output "convergence_data.png"
set format x "%1.1E"
set rmargin 6
#plot [][600:] 'convergence_data.txt' w lp pt 6 lw 2 lc rgb 'blue' t 'Max temperature at t=1e-6s'
#plot [][600:] 'convergence_data.txt' using 1:2:(sprintf("(%d, %d)", $1, $2)) with labels point pt 7 offset char 1,1 notitle with lines
plot [][600:] "convergence_data.txt" w lp pt 6 lw 2 lc rgb 'blue' t 'Max temperature at t=1e-6s' \
   , '' u 1:2:(sprintf("(%d, %d)", $1, $2)) with labels point pt 7 offset char 1,1 notitle
set output
