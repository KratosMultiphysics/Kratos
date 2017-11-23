set grid
set key left bottom
set xlabel 'Tangent of incident angle'
set ylabel 'Tangent of recoil angle'
set style line 1 pt 8 lt -1 ps 3
set style line 2 pt 9 lt 3 ps 3
plot [0:7][-2:8] 'benchmark6_dt_1e-06_Vst_prima_div_Vcn_prima_data.dat' w lp lt 1 lw 1.5 ps 2 pt 5,\
'paper_data/benchmark6_graph2.dat' index 0 w lp ls 1 t 'Al. alloy',\
'paper_data/benchmark6_graph2.dat' index 1 w lp ls 2 t 'Nylon'
