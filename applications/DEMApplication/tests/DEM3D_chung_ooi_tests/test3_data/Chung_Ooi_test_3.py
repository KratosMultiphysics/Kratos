import KratosMultiphysics as Kratos
from KratosMultiphysics.DEMApplication.DEM_analysis_stage import DEMAnalysisStage
import KratosMultiphysics.DEMApplication.plot_variables as plot_variables
import KratosMultiphysics.KratosUnittest as KratosUnittest
import os
class ChungOoiTest3(KratosUnittest.TestCase):

    def __init__(self, *args, **kwargs):
        super(KratosUnittest.TestCase, self).__init__(*args, **kwargs)
        self.remove_all_results = True
        self.measured_restitution_numbers = []

    def GetInputParameters(self, iteration):
        os.chdir(os.path.dirname(os.path.realpath(__file__)))
        file_name = "ProjectParameters3" + str(iteration) + ".json"
        with open(file_name, 'r') as parameters_file:
            parameters = Kratos.Parameters(parameters_file.read())
        return parameters

    def test_Run(self):
        for iteration in range(1, 7):
            parameters = self.GetInputParameters(iteration)
            iteration_case = DEMAnalysisStageForChungOoiTest3(Kratos.Model(), parameters)
            iteration_case.Run()
            self.measured_restitution_numbers.append(iteration_case.restitution_coefficient)

        self.PrintResultsAfterAllComputations(parameters["MaxTimeStep"].GetDouble())
        self.CheckResults()
        if self.remove_all_results:
            plot_variables.delete_archives()

    def CheckResults(self):
        ideal_CR_list = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
        for i in range(0, 6):
            self.assertAlmostEqual(ideal_CR_list[i], self.measured_restitution_numbers[i], delta = 1e-3)

    def PrintResultsAfterAllComputations(self, dt):
        self.restitution_numbers_vector_list_outfile_name = "benchmark3_dt_" + str(dt) + '_restitution_numbers_vector_list_data.dat'
        self.restitution_numbers_vector_list_outfile = open(self.restitution_numbers_vector_list_outfile_name, 'w')
        number_of_points_in_graph = len(self.measured_restitution_numbers)
        for i in range(0, number_of_points_in_graph):
            first_col = 1 / (number_of_points_in_graph - 1) * i
            self.restitution_numbers_vector_list_outfile.write("%6.4f %11.8f" % (first_col, self.measured_restitution_numbers[i]) + '\n')
        self.restitution_numbers_vector_list_outfile.close()
        gnuplot_script_name = 'benchmark3_dt_' + str(dt) + 's_CORs.gp'
        self.gnuplot_outfile = open(gnuplot_script_name, 'w')
        self.gnuplot_outfile.write("set grid; set key center right;  set xlabel 'Restitution coefficient'; set ylabel 'Vn_{final}/Vn_{initial}';\\\n")
        self.gnuplot_outfile.write("plot '" + self.restitution_numbers_vector_list_outfile_name + "' u 1:2 w lp lt 3 lw 1.5 ps 2 pt 4 t 'DEM',\\\n")
        self.gnuplot_outfile.write("'paper_data/benchmark3_graph1.dat' w lp ls 2 ps 1.5 pt 9 t 'reference',\\\n")
        self.gnuplot_outfile.close()
        #COC.print_gnuplot_files_on_screen(gnuplot_script_name)

class DEMAnalysisStageForChungOoiTest3(DEMAnalysisStage):

    def Initialize(self):
        super().Initialize()
        self.plotter = plot_variables.variable_plotter(self.spheres_model_part, [1])

    def Finalize(self):
        for node in self.spheres_model_part.Nodes:
            final_vel = node.GetSolutionStepValue(Kratos.VELOCITY_Y)
        self.restitution_coefficient = final_vel / 3.9
        self.plotter.close_files()
        super().Finalize()

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        self.plotter.plot_variables(self.time)


if __name__ == "__main__":
    Kratos.Logger.GetDefaultOutput().SetSeverity(Kratos.Logger.Severity.WARNING)
    KratosUnittest.main()
