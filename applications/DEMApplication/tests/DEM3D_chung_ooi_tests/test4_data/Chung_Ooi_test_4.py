from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.DEMApplication.DEM_analysis_stage import DEMAnalysisStage
import KratosMultiphysics.DEMApplication.plot_variables as plot_variables
import KratosMultiphysics.DEMApplication.Chung_Ooi_class as COC
from math import pi, sin, cos, atan

remove_all_results = True

def GetInputParameters():
    file_name = "ProjectParameters4.json"
    mat_name = "Materials4.json"
    with open(file_name, 'r') as parameters_file:
        parameters = Parameters(parameters_file.read())
    with open(mat_name,'r') as mat_parameter_file:
        mat_parameters = Parameters(mat_parameter_file.read())
    return parameters, mat_parameters

class Test4:

    def __init__(self):
        self.initial_module_vel = 3.9
        self.initial_tangential_vel = 0
        self.radius = 0.0025
        self.degrees = 0
        self.angles_list = []
        self.tangential_restitution_coefficient_list = []
        self.final_angular_vel_list = []
        self.rebound_angle_list = []
        self.final_angular_vel_list_outfile = None
        self.rebound_angle_list_outfile = None
        self.tangential_restitution_coefficient_list_outfile = None

    def set_initial_data(self, modelpart, iteration, number_of_points_in_the_graphic):

        self.degrees = 90 / (number_of_points_in_the_graphic + 1) * iteration
        self.initial_tangential_vel =  -self.initial_module_vel * sin(self.degrees * pi / 180.0)
        initial_normal_vel = -self.initial_module_vel * cos(self.degrees * pi / 180.0)
        for node in modelpart.Nodes:
            node.SetSolutionStepValue(VELOCITY_Y, initial_normal_vel)
            node.SetSolutionStepValue(VELOCITY_Z, self.initial_tangential_vel)

    def get_final_data(self, modelpart):

        for node in modelpart.Nodes:
            final_angular_vel = node.GetSolutionStepValue(ANGULAR_VELOCITY_X)
            final_normal_center_velocity = node.GetSolutionStepValue(VELOCITY_Y)
            final_tangential_center_velocity = node.GetSolutionStepValue(VELOCITY_Z)
            final_tangential_contact_velocity = final_tangential_center_velocity - final_angular_vel * self.radius
            rebound_angle = 180 / pi * atan(final_tangential_contact_velocity / -final_normal_center_velocity)
            tangential_restitution_coefficient = final_tangential_center_velocity / self.initial_tangential_vel
        self.final_angular_vel_list.append(final_angular_vel)
        self.rebound_angle_list.append(rebound_angle)
        self.tangential_restitution_coefficient_list.append(tangential_restitution_coefficient)
        self.angles_list.append(self.degrees)

    def print_results(self, number_of_points_in_the_graphic, dt=0):

        self.tangential_restitution_coefficient_list_outfile_name = "benchmark4_dt_" + str(dt) + '_tangential_restitution_coefficient_list_data.dat'
        self.final_angular_vel_list_outfile_name = "benchmark4_dt_" + str(dt) + '_final_angular_vel_list_data.dat'
        self.rebound_angle_list_outfile_name = "benchmark4_dt_" + str(dt) + '_rebound_angle_list_data.dat'
        self.tangential_restitution_coefficient_list_outfile = open(self.tangential_restitution_coefficient_list_outfile_name, 'w')
        self.final_angular_vel_list_outfile = open(self.final_angular_vel_list_outfile_name, 'w')
        self.rebound_angle_list_outfile = open(self.rebound_angle_list_outfile_name, 'w')
        for i in range(0, number_of_points_in_the_graphic):
            self.tangential_restitution_coefficient_list_outfile.write("%14.8f %14.8f" % (self.angles_list[i], self.tangential_restitution_coefficient_list[i]) + '\n')
            self.final_angular_vel_list_outfile.write("%14.8f %14.8f" % (self.angles_list[i], self.final_angular_vel_list[i]) + '\n')
            self.rebound_angle_list_outfile.write("%14.8f %14.8f" % (self.angles_list[i], self.rebound_angle_list[i]) + '\n')
        self.tangential_restitution_coefficient_list_outfile.close()
        self.final_angular_vel_list_outfile.close()
        self.rebound_angle_list_outfile.close()
        self.create_gnuplot_scripts(self.tangential_restitution_coefficient_list_outfile_name, self.final_angular_vel_list_outfile_name,\
                                    self.rebound_angle_list_outfile_name, dt)

    def create_gnuplot_scripts(self, tangential_restitution_coefficient_list_outfile_name, final_angular_vel_list_outfile_name,\
                               rebound_angle_list_outfile_name, dt):

        gnuplot_script_name_1 = 'benchmark4_comparison_1_dt_' + str(dt) + 's.gp'
        self.gnuplot_outfile = open(gnuplot_script_name_1, 'w')
        self.gnuplot_outfile.write("set grid\nset key right center\nset xlabel 'Incident angle (deg)'\nset ylabel 'Tangential restitution coefficient'\nset style line 1 pt 8 lt -1 ps 3\nset style line 2 pt 9 lt 3 ps 3\n")
        self.gnuplot_outfile.write("plot [0:90][.4:1] '" + tangential_restitution_coefficient_list_outfile_name + "' w lp lt 1 lw 1.5 ps 2 pt 5 t 'DEM',\\\n")
        self.gnuplot_outfile.write("'paper_data/benchmark4_graph1.dat' index 0 w lp ls 1 t 'reference',\\\n")
        self.gnuplot_outfile.close()
        COC.print_gnuplot_files_on_screen(gnuplot_script_name_1)
        gnuplot_script_name_2 = 'benchmark4_comparison_2_dt_' + str(dt) + 's.gp'
        self.gnuplot_outfile = open(gnuplot_script_name_2, 'w')
        self.gnuplot_outfile.write("set grid\nset key right center\nset xlabel 'Incident angle (deg)'\nset ylabel 'Final angular velocity (rad/s)'\nset style line 1 pt 8 lt -1 ps 3\nset style line 2 pt 9 lt 3 ps 3\n")
        self.gnuplot_outfile.write("plot [0:90][-800:0] '" + final_angular_vel_list_outfile_name + "' w lp lt 1 lw 1.5 ps 2 pt 5 t 'DEM',\\\n")
        self.gnuplot_outfile.write("'paper_data/benchmark4_graph2.dat' index 0 w lp ls 1 t 'reference',\\\n")
        self.gnuplot_outfile.close()
        COC.print_gnuplot_files_on_screen(gnuplot_script_name_2)
        gnuplot_script_name_3 = 'benchmark4_comparison_3_dt_' + str(dt) + 's.gp'
        self.gnuplot_outfile = open(gnuplot_script_name_3, 'w')
        self.gnuplot_outfile.write("set grid\nset key right center\nset xlabel 'Incident angle (deg)'\nset ylabel 'Rebound angle (deg)'\nset style line 1 pt 8 lt -1 ps 3\nset style line 2 pt 9 lt 3 ps 3\n")
        self.gnuplot_outfile.write("plot [0:90][-30:90] '" + rebound_angle_list_outfile_name + "' w lp lt 1 lw 1.5 ps 2 pt 5 t 'DEM',\\\n")
        self.gnuplot_outfile.write("'paper_data/benchmark4_graph3.dat' index 0 w lp ls 1 t 'reference',\\\n")
        self.gnuplot_outfile.close()
        COC.print_gnuplot_files_on_screen(gnuplot_script_name_3)

test = Test4()

class DEMAnalysisStageForChungOoiTest4(DEMAnalysisStage):

    def __init__(self, model, DEM_parameters, iteration = 0, number_of_points_in_the_graphic = 0):
        super().__init__(model, DEM_parameters)
        self.iteration = iteration
        self.number_of_points_in_the_graphic = number_of_points_in_the_graphic

    def Initialize(self):
        super().Initialize()
        self.plotter = plot_variables.variable_plotter(self.spheres_model_part, [1])

    def ReadModelParts(self):
        super().ReadModelParts()
        test.set_initial_data(self.spheres_model_part, self.iteration, self.number_of_points_in_the_graphic)

    def Finalize(self):
        test.get_final_data(self.spheres_model_part)
        self.plotter.close_files()
        super().Finalize()

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        self.plotter.plot_variables(self.time)

number_of_points_in_the_graphic = 17
for iteration in range(1, number_of_points_in_the_graphic + 1):
    parameters, _ = GetInputParameters()
    DEMAnalysisStageForChungOoiTest4(Model(), parameters, iteration, number_of_points_in_the_graphic).Run()
test.print_results(number_of_points_in_the_graphic, parameters["MaxTimeStep"].GetDouble())
if remove_all_results:
    plot_variables.delete_archives()
