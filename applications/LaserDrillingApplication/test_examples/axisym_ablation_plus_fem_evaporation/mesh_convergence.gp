set grid lw 2
set key right center font "Times-Roman,16"
set xtics font "Times-Roman,15"
set ytics font "Times-Roman,15"
set xlabel "Number of elements" font "Times-Roman,17"
set ylabel "Temperature (K)" font "Times-Roman,17"
set term pngcairo size 900,600
set output "convergence_data_both.png"
set format x "%1.1E"
set rmargin 6
set style line 1 \
    linecolor rgb '#0060ad' \
    linetype 1 linewidth 2 \
    pointtype 7 pointsize 1.5
set style line 2 \
    linecolor rgb '#dd181f' \
    linetype 1 linewidth 2 \
    pointtype 5 pointsize 1.5
plot [][600:] "convergence_data_structured.txt" with linespoints linestyle 1 t 'Structured, Tmax at t=1e-6s' \
   , '' u 1:2:(sprintf("(%d, %d)", $1, $2)) with labels offset char -1,1 notitle \
   , "convergence_data_unstructured.txt" with linespoints linestyle 2 t 'Unstructured, Tmax at t=1e-6s' \
   , '' u 1:2:(sprintf("(%d, %d)", $1, $2)) with labels offset char -1,1 notitle
set output
