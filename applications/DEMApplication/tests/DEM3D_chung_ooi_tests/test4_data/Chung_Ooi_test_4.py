import KratosMultiphysics as Kratos
from KratosMultiphysics.DEMApplication.DEM_analysis_stage import DEMAnalysisStage
import KratosMultiphysics.DEMApplication.plot_variables as plot_variables
import KratosMultiphysics.DEMApplication.Chung_Ooi_class as COC
import KratosMultiphysics.KratosUnittest as KratosUnittest
from math import pi, sin, cos, atan
import os

class ChungOoiTest4(KratosUnittest.TestCase):

    def __init__(self, *args, **kwargs):
        super(KratosUnittest.TestCase, self).__init__(*args, **kwargs)
        self.remove_all_results = True
        self.initial_tangential_vel = 0
        self.degrees = 0
        self.angles_list = []
        self.tangential_restitution_coefficient_list = []
        self.final_angular_vel_list = []
        self.rebound_angle_list = []
        self.final_angular_vel_list_outfile = None
        self.rebound_angle_list_outfile = None
        self.tangential_rest_coeff_list_outfile = None

    def GetInputParameters(self):
        os.chdir(os.path.dirname(os.path.realpath(__file__)))
        file_name = "ProjectParameters4.json"
        with open(file_name, 'r') as parameters_file:
            parameters = Kratos.Parameters(parameters_file.read())
        return parameters

    def test_Run(self):
        parameters = self.GetInputParameters()
        for iteration in range(1, 18):
            iteration_case = DEMAnalysisStageForChungOoiTest4(Kratos.Model(), parameters, iteration)
            iteration_case.Run()
            self.final_angular_vel_list.append(iteration_case.final_angular_vel)
            self.rebound_angle_list.append(iteration_case.rebound_angle)
            self.tangential_restitution_coefficient_list.append(iteration_case.tangential_restitution_coefficient)
            self.angles_list.append(iteration_case.degrees)

        self.PrintResultsAfterAllComputations(parameters["MaxTimeStep"].GetDouble())
        self.CheckResults()
        if self.remove_all_results:
            plot_variables.delete_archives()

    def CheckResults(self):
        self.ComputeErrors()
        self.assertAlmostEqual(self.error1, 0.0, delta = 1.6e-2)
        self.assertAlmostEqual(self.error2, 0.0, delta = 2.2e-2)
        self.assertAlmostEqual(self.error3, 0.0, delta = 2.5e-2)

    def PrintResultsAfterAllComputations(self, dt):
        self.tangential_restitution_coefficient_list_outfile_name = "benchmark4_dt_" + str(dt) + '_tangential_restitution_coefficient_list_data.dat'
        self.final_angular_vel_list_outfile_name = "benchmark4_dt_" + str(dt) + '_final_angular_vel_list_data.dat'
        self.rebound_angle_list_outfile_name = "benchmark4_dt_" + str(dt) + '_rebound_angle_list_data.dat'
        self.tangential_rest_coeff_list_outfile = open(self.tangential_restitution_coefficient_list_outfile_name, 'w')
        self.final_angular_vel_list_outfile = open(self.final_angular_vel_list_outfile_name, 'w')
        self.rebound_angle_list_outfile = open(self.rebound_angle_list_outfile_name, 'w')
        for i in range(0, 17):
            self.tangential_rest_coeff_list_outfile.write("%14.8f %14.8f" % (self.angles_list[i], self.tangential_restitution_coefficient_list[i]) + '\n')
            self.final_angular_vel_list_outfile.write("%14.8f %14.8f" % (self.angles_list[i], self.final_angular_vel_list[i]) + '\n')
            self.rebound_angle_list_outfile.write("%14.8f %14.8f" % (self.angles_list[i], self.rebound_angle_list[i]) + '\n')
        self.tangential_rest_coeff_list_outfile.close()
        self.final_angular_vel_list_outfile.close()
        self.rebound_angle_list_outfile.close()
        self.create_gnuplot_scripts(self.tangential_restitution_coefficient_list_outfile_name, self.final_angular_vel_list_outfile_name,\
                                    self.rebound_angle_list_outfile_name, dt)

    def ComputeErrors(self):
        import numpy as np
        reference_data = np.loadtxt('paper_data/benchmark4_graph1_ref.dat', usecols = (0, 1))
        dem_results = np.loadtxt(self.tangential_restitution_coefficient_list_outfile_name, usecols = (0, 1))
        i_meas, i_truth = np.where(np.isclose(reference_data[:, None, 0], dem_results[:, 0], atol = 0.5))
        a = reference_data[i_meas][:, 1]
        b = dem_results[i_truth][:, 1]
        self.error1 = COC.ComputeRelativeError(a, b)
        reference_data = np.loadtxt('paper_data/benchmark4_graph2_ref.dat', usecols = (0, 1))
        dem_results = np.loadtxt(self.final_angular_vel_list_outfile_name, usecols = (0, 1))
        i_meas, i_truth = np.where(np.isclose(reference_data[:, None, 0], dem_results[:, 0], atol = 0.5))
        a = reference_data[i_meas][:, 1]
        b = dem_results[i_truth][:, 1]
        self.error2 = COC.ComputeRelativeError(a, b)
        reference_data = np.loadtxt('paper_data/benchmark4_graph3_ref.dat', usecols = (0, 1))
        dem_results = np.loadtxt(self.rebound_angle_list_outfile_name, usecols = (0, 1))
        i_meas, i_truth = np.where(np.isclose(reference_data[:, None, 0], dem_results[:, 0], atol = 0.5))
        a = reference_data[i_meas][:, 1]
        b = dem_results[i_truth][:, 1]
        self.error3 = COC.ComputeRelativeError(a, b)

    def create_gnuplot_scripts(self, tangential_restitution_coefficient_list_outfile_name, final_angular_vel_list_outfile_name,\
                               rebound_angle_list_outfile_name, dt):

        gnuplot_script_name_1 = 'benchmark4_comparison_1_dt_' + str(dt) + 's.gp'
        self.gnuplot_outfile = open(gnuplot_script_name_1, 'w')
        self.gnuplot_outfile.write("set grid\nset key right center\nset xlabel 'Incident angle (deg)'\nset ylabel 'Tangential restitution coefficient';\\\n")
        self.gnuplot_outfile.write("set style line 1 pt 8 lt -1 ps 3\nset style line 2 pt 9 lt 3 ps 3\n")
        self.gnuplot_outfile.write("plot [0:90][.4:1] '" + tangential_restitution_coefficient_list_outfile_name + "' w lp lt 1 lw 1.5 ps 2 pt 5 t 'DEM',\\\n")
        self.gnuplot_outfile.write("'paper_data/benchmark4_graph1.dat' index 0 w lp ls 1 t 'reference',\\\n")
        self.gnuplot_outfile.close()
        #COC.print_gnuplot_files_on_screen(gnuplot_script_name_1)
        gnuplot_script_name_2 = 'benchmark4_comparison_2_dt_' + str(dt) + 's.gp'
        self.gnuplot_outfile = open(gnuplot_script_name_2, 'w')
        self.gnuplot_outfile.write("set grid\nset key right center\nset xlabel 'Incident angle (deg)'\nset ylabel 'Final angular velocity (rad/s)';\\\n")
        self.gnuplot_outfile.write("set style line 1 pt 8 lt -1 ps 3\nset style line 2 pt 9 lt 3 ps 3\n")
        self.gnuplot_outfile.write("plot [0:90][-800:0] '" + final_angular_vel_list_outfile_name + "' w lp lt 1 lw 1.5 ps 2 pt 5 t 'DEM',\\\n")
        self.gnuplot_outfile.write("'paper_data/benchmark4_graph2.dat' index 0 w lp ls 1 t 'reference',\\\n")
        self.gnuplot_outfile.close()
        #COC.print_gnuplot_files_on_screen(gnuplot_script_name_2)
        gnuplot_script_name_3 = 'benchmark4_comparison_3_dt_' + str(dt) + 's.gp'
        self.gnuplot_outfile = open(gnuplot_script_name_3, 'w')
        self.gnuplot_outfile.write("set grid\nset key right center\nset xlabel 'Incident angle (deg)'\nset ylabel 'Rebound angle (deg)';\\\n")
        self.gnuplot_outfile.write("set style line 1 pt 8 lt -1 ps 3\nset style line 2 pt 9 lt 3 ps 3\n")
        self.gnuplot_outfile.write("plot [0:90][-30:90] '" + rebound_angle_list_outfile_name + "' w lp lt 1 lw 1.5 ps 2 pt 5 t 'DEM',\\\n")
        self.gnuplot_outfile.write("'paper_data/benchmark4_graph3.dat' index 0 w lp ls 1 t 'reference',\\\n")
        self.gnuplot_outfile.close()
        #COC.print_gnuplot_files_on_screen(gnuplot_script_name_3)

class DEMAnalysisStageForChungOoiTest4(DEMAnalysisStage):

    def __init__(self, model, DEM_parameters, iteration):
        super().__init__(model, DEM_parameters)
        self.iteration = iteration

    def Initialize(self):
        super().Initialize()
        self.set_initial_data(self.spheres_model_part, self.iteration)
        self.plotter = plot_variables.variable_plotter(self.spheres_model_part, [1])

    def set_initial_data(self, modelpart, iteration):
        self.degrees = (90 / 18) * iteration
        self.initial_tangential_vel =  -3.9 * sin(self.degrees * pi / 180.0)
        initial_normal_vel = -3.9 * cos(self.degrees * pi / 180.0)
        for node in modelpart.Nodes:
            node.SetSolutionStepValue(Kratos.VELOCITY_Y, initial_normal_vel)
            node.SetSolutionStepValue(Kratos.VELOCITY_Z, self.initial_tangential_vel)

    def GetFinalData(self, modelpart):
        for node in modelpart.Nodes:
            self.final_angular_vel = node.GetSolutionStepValue(Kratos.ANGULAR_VELOCITY_X)
            self.final_normal_center_velocity = node.GetSolutionStepValue(Kratos.VELOCITY_Y)
            self.final_tangential_center_velocity = node.GetSolutionStepValue(Kratos.VELOCITY_Z)
            self.final_tangential_contact_velocity = self.final_tangential_center_velocity - self.final_angular_vel * 0.0025
            self.rebound_angle = 180 / pi * atan(self.final_tangential_contact_velocity / -self.final_normal_center_velocity)
            self.tangential_restitution_coefficient = self.final_tangential_center_velocity / self.initial_tangential_vel

    def Finalize(self):
        self.GetFinalData(self.spheres_model_part)
        self.plotter.close_files()
        super().Finalize()

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        self.plotter.plot_variables(self.time)


if __name__ == "__main__":
    Kratos.Logger.GetDefaultOutput().SetSeverity(Kratos.Logger.Severity.WARNING)
    KratosUnittest.main()
