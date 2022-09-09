from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.DEMApplication.DEM_analysis_stage import DEMAnalysisStage
import KratosMultiphysics.DEMApplication.plot_variables as plot_variables
import KratosMultiphysics.DEMApplication.Chung_Ooi_class as COC

remove_all_results = True

def GetInputParameters(iteration):
    file_name = "ProjectParameters3" + str(iteration) + ".json"
    mat_name = "Materials3" + str(iteration) + ".json"
    with open(file_name, 'r') as parameters_file:
        parameters = Parameters(parameters_file.read())
    with open(mat_name,'r') as mat_parameter_file:
        mat_parameters = Parameters(mat_parameter_file.read())
    return parameters, mat_parameters

class Test3:

    def __init__(self):
        self.restitution_numbers_list = []
        self.initial_normal_vel = -3.9
        self.number_of_points_in_graph = 6

    def get_final_data(self, modelpart):

        for node in modelpart.Nodes:
            final_vel = node.GetSolutionStepValue(VELOCITY_Y)
        restitution_coefficient = -final_vel / self.initial_normal_vel
        self.restitution_numbers_list.append(restitution_coefficient)

    def print_results(self, dt=0):

        self.restitution_numbers_vector_list_outfile_name = "benchmark3_dt_" + str(dt) + '_restitution_numbers_vector_list_data.dat'
        self.restitution_numbers_vector_list_outfile = open(self.restitution_numbers_vector_list_outfile_name, 'w')

        for i in range(0, self.number_of_points_in_graph):
            first_col = 1 / (self.number_of_points_in_graph - 1) * i
            self.restitution_numbers_vector_list_outfile.write("%6.4f %11.8f" % (first_col, self.restitution_numbers_list[i]) + '\n')
        self.restitution_numbers_vector_list_outfile.close()

        gnuplot_script_name = 'benchmark3_dt_' + str(dt) + 's_CORs.gp'
        self.gnuplot_outfile = open(gnuplot_script_name, 'w')
        self.gnuplot_outfile.write("set grid; set key center right;  set xlabel 'Restitution coefficient'; set ylabel 'Vn_{final}/Vn_{initial}'; plot '" + self.restitution_numbers_vector_list_outfile_name + "' u 1:2 w lp lt 3 lw 1.5 ps 2 pt 4 t 'DEM',\\\n")
        self.gnuplot_outfile.write("'paper_data/benchmark3_graph1.dat' w lp ls 2 ps 1.5 pt 9 t 'reference',\\\n")
        self.gnuplot_outfile.close()
        COC.print_gnuplot_files_on_screen(gnuplot_script_name)

test = Test3()

class DEMAnalysisStageForChungOoiTest3(DEMAnalysisStage):

    def __init__(self, model, DEM_parameters, iteration, number_of_points_in_the_graphic):
        super().__init__(model, DEM_parameters)

        self.iteration = iteration
        self.number_of_points_in_the_graphic = number_of_points_in_the_graphic

    def Initialize(self):
        super().Initialize()
        self.plotter = plot_variables.variable_plotter(self.spheres_model_part, [1])

    def Finalize(self):
        test.get_final_data(self.spheres_model_part)
        self.plotter.close_files()
        super().Finalize()

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        self.plotter.plot_variables(self.time)

list_of_coeffs_of_restitution = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
number_of_points_in_the_graphic = len(list_of_coeffs_of_restitution)
iteration = 1
for coeff in list_of_coeffs_of_restitution:
    parameters, mat_parameters = GetInputParameters(iteration)
    DEMAnalysisStageForChungOoiTest3(Model(), parameters, iteration, number_of_points_in_the_graphic).Run()
    iteration += 1
test.print_results(parameters["MaxTimeStep"].GetDouble())
if remove_all_results:
    plot_variables.delete_archives()
