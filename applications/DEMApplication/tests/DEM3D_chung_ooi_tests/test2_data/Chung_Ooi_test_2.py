from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.DEMApplication.DEM_analysis_stage import DEMAnalysisStage
import KratosMultiphysics.DEMApplication.plot_variables as plot_variables
import KratosMultiphysics.DEMApplication.Chung_Ooi_class as COC

remove_all_results = True

def print_results(dt):
    gnuplot_script_name1 = 'benchmark2_dt_' + str(dt) + 's_time_vs_normal_force.gp'
    gnuplot_outfile1 = open(gnuplot_script_name1, 'w')
    gnuplot_outfile1.write("set grid; set key center right; set xlabel 'Time (us)'; set ylabel 'Normal contact force (KN)'; plot [0:800][0:12] 'variables_for_node_1.dat' every 20 u (1e6*$1):(1e-3*$9) w lp ls 1 ps 1.5 pt 5 t 'DEM',\\\n")
    gnuplot_outfile1.write("'paper_data/benchmark2_graph2.dat' w lp ls 2 ps 1.5 pt 9 t 'reference',\\\n")
    gnuplot_outfile1.close()
    COC.print_gnuplot_files_on_screen(gnuplot_script_name1)
    gnuplot_script_name2 = 'benchmark2_dt_' + str(dt) + 's_indentation_vs_normal_force.gp'
    gnuplot_outfile2 = open(gnuplot_script_name2, 'w')
    gnuplot_outfile2.write("set grid; set key center right; set xlabel 'Normal contact displacement (um)'; set ylabel 'Normal contact force (KN)'; plot [0:60][0:12] 'variables_for_node_1.dat' every 20 u (-1e6*$3):(1e-3*$9) w lp ls 1 ps 1.5 pt 5 t 'DEM',\\\n")
    gnuplot_outfile2.write("'paper_data/benchmark2_graph1.dat' w lp ls 2 ps 1.5 pt 9 t 'reference',\\\n")
    gnuplot_outfile2.close()
    COC.print_gnuplot_files_on_screen(gnuplot_script_name2)

def GetInputParameters():
    file_name = "ProjectParameters2.json"
    mat_name = "Materials2.json"
    with open(file_name, 'r') as parameters_file:
        parameters = Parameters(parameters_file.read())
    with open(mat_name,'r') as mat_parameter_file:
        mat_parameters = Parameters(mat_parameter_file.read())
    return parameters, mat_parameters

class DEMAnalysisStageForChungOoiTest2(DEMAnalysisStage):

    def Initialize(self):
        super().Initialize()
        self.plotter = plot_variables.variable_plotter(self.spheres_model_part, [1])

    def Finalize(self):
        super().Finalize()
        self.plotter.close_files()

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        self.plotter.plot_variables(self.time)

parameters, _ = GetInputParameters()
DEMAnalysisStageForChungOoiTest2(Model(), parameters).Run()
print_results(parameters["MaxTimeStep"].GetDouble())
if remove_all_results:
    plot_variables.delete_archives()
