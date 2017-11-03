set grid
set key left bottom
set style line 1 pt 8 lt -1 ps 3
set style line 2 pt 9 lt 3 ps 3
plot [0:14][-4:6] 'benchmark5_dt_3.6e-11_Vst_prima_div_mu_per_Vcn_prima_list_data.dat' w lp lt 1 lw 1.5 ps 2 pt 5,\
'paper_data/benchmark5_graph1.dat' index 0 w lp ls 1 t 'Steel',\
'paper_data/benchmark5_graph1.dat' index 1 w lp ls 2 t 'Polyethylene',\
'paper_data/benchmark5_graph1.dat' index 2 w p pt 7 ps 2 lt -1 t 'FEM'
