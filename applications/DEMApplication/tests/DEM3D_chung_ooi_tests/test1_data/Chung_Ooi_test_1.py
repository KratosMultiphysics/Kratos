import KratosMultiphysics as Kratos
from KratosMultiphysics.DEMApplication.DEM_analysis_stage import DEMAnalysisStage
import KratosMultiphysics.DEMApplication.plot_variables as plot_variables
import KratosMultiphysics.DEMApplication.Chung_Ooi_class as COC
import KratosMultiphysics.KratosUnittest as KratosUnittest
import os
class ChungOoiTest1(KratosUnittest.TestCase):

    def __init__(self, *args, **kwargs):
        super(KratosUnittest.TestCase, self).__init__(*args, **kwargs)
        self.remove_all_results = True

    def GetInputParameters(self):
        os.chdir(os.path.dirname(os.path.realpath(__file__)))
        file_name = "ProjectParameters1.json"
        with open(file_name, 'r') as parameters_file:
            parameters = Kratos.Parameters(parameters_file.read())
        return parameters

    def test_Run(self):
        parameters = self.GetInputParameters()
        DEMAnalysisStageForChungOoiTest1(Kratos.Model(), parameters).Run()
        self.PrintResultsAfterAllComputations(parameters["MaxTimeStep"].GetDouble())
        self.CheckResults()
        if self.remove_all_results:
            plot_variables.delete_archives()

    def CheckResults(self):
        self.ComputeErrors()
        self.assertAlmostEqual(self.error1, 0.0, delta = 3.1e-3)
        self.assertAlmostEqual(self.error2, 0.0, delta = 2e-3)

    def PrintResultsAfterAllComputations(self, dt):
        gnuplot_script_name1 = 'benchmark1_dt_' + str(dt) + 's_time_vs_normal_force.gp'
        gnuplot_outfile1 = open(gnuplot_script_name1, 'w')
        gnuplot_outfile1.write("set grid; set key center right; set xlabel 'Time (us)'; set ylabel 'Normal contact force (KN)';\\\n")
        gnuplot_outfile1.write("plot [0:60][0:12] 'variables_for_node_1.dat' every 20 u (1e6*$1):(1e-3*$9) w lp ls 1 ps 1.5 pt 5 t 'DEM',\\\n")
        gnuplot_outfile1.write("'paper_data/benchmark1_graph2.dat' w lp ls 2 ps 1.5 pt 9 t 'reference',\\\n")
        gnuplot_outfile1.close()
        #COC.print_gnuplot_files_on_screen(gnuplot_script_name1)
        gnuplot_script_name2 = 'benchmark1_dt_' + str(dt) + 's_indentation_vs_normal_force.gp'
        gnuplot_outfile2 = open(gnuplot_script_name2, 'w')
        gnuplot_outfile2.write("set grid; set key center right; set xlabel 'Normal contact displacement (um)'; set ylabel 'Normal contact force (KN)';\\\n")
        gnuplot_outfile2.write("plot [0:400][0:12] 'variables_for_node_1.dat' every 20 u (-2e6*$3):(1e-3*$9) w lp ls 1 ps 1.5 pt 5 t 'DEM',\\\n")
        gnuplot_outfile2.write("'paper_data/benchmark1_graph1.dat' w lp ls 2 ps 1.5 pt 9 t 'reference',\\\n")
        gnuplot_outfile2.close()
        #COC.print_gnuplot_files_on_screen(gnuplot_script_name2)

    def ComputeErrors(self):
        import numpy as np
        reference_data = np.loadtxt('paper_data/benchmark1_graph2.dat', usecols = (0, 1))
        dem_results = np.loadtxt('variables_for_node_1.dat', usecols = (0, 8), skiprows = 1)
        dem_results[:, 0] *= 1e6
        dem_results[:, 1] *= 1e-3
        i_meas, i_truth = np.where(np.isclose(reference_data[:, None, 0], dem_results[:, 0], atol = 0.015))
        a = reference_data[i_meas][:, 1]
        b = dem_results[i_truth][:, 1]
        self.error1 = COC.ComputeRelativeError(a, b)
        reference_data = np.loadtxt('paper_data/benchmark1_graph1.dat', usecols = (0, 1))
        dem_results = np.loadtxt('variables_for_node_1.dat', usecols = (2, 8), skiprows = 1)
        dem_results[:, 0] *= -2e6
        dem_results[:, 1] *= 1e-3
        i_meas, i_truth = np.where(np.isclose(reference_data[:, None, 0], dem_results[:, 0], atol = 0.05))
        a = reference_data[i_meas][:, 1]
        b = dem_results[i_truth][:, 1]
        self.error2 = COC.ComputeRelativeError(a, b)

class DEMAnalysisStageForChungOoiTest1(DEMAnalysisStage):

    def Initialize(self):
        super().Initialize()
        self.plotter = plot_variables.variable_plotter(self.spheres_model_part, [1])

    def Finalize(self):
        super().Finalize()
        self.plotter.close_files()

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        self.plotter.plot_variables(self.time)


if __name__ == "__main__":
    Kratos.Logger.GetDefaultOutput().SetSeverity(Kratos.Logger.Severity.WARNING)
    KratosUnittest.main()
