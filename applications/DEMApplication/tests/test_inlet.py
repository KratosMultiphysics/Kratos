import os
import numpy as np
import KratosMultiphysics as Kratos
from Kratos import Logger
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.DEMApplication.DEM_analysis_stage as dem_analysis

import random_variable_tests_files.random_variable as rv

debug_mode = False

class TestDEMInletAnalysis(dem_analysis.DEMAnalysisStage):

    @staticmethod
    def Say(*args):
        Logger.PrintInfo("DEM", *args)
        Logger.Flush()

    @staticmethod # This error is bounded by 1
    def Error(values, values_ref, interval_widths):
        return 0.5 * sum(abs((values_ref - values)) * interval_widths)

    @classmethod
    def GetMainPath(self):
        return os.path.join(os.path.dirname(os.path.realpath(__file__)), "piecewise_linear_inlet_tests_files")

    def __init__(self, model, DEM_parameters):
        self.tolerance = 0.05
        self.error = self.tolerance + 1
        self.n_bins = 25 # for the histogram

        super().__init__(model, DEM_parameters)

    def Initialize(self):
        self.samples = dict()
        return super().Initialize()

    def InitializeSolutionStep(self):
        return super().InitializeSolutionStep()

    def HasConverged(self):
        return self.error < self.tolerance

    def KeepAdvancingSolutionLoop(self):
        return super().KeepAdvancingSolutionLoop() and not self.HasConverged()

    def GetProblemNameWithPath(self):
        return os.path.join(self.main_path, self.DEM_parameters["problem_name"].GetString())

    def FinalizeSolutionStep(self):
        particle_nodes = [node for node in self.spheres_model_part.Nodes if node.IsNot(Kratos.BLOCKED)]

        # Mapping to Id so it works even if some particles are removed before the end of the simulation
        for node in particle_nodes:
            self.samples[node.Id] = node.GetSolutionStepValue(Kratos.RADIUS)

        if particle_nodes:
            self.empirical_pdf, interval_boundaries = np.histogram(list(self.samples.values()), self.n_bins, density=True)
            self.centers = 0.5 * (interval_boundaries[1:] + interval_boundaries[:-1])
            interval_widths = interval_boundaries[1:] - interval_boundaries[:-1]
            breakpoints = self.DEM_parameters['dem_inlets_settings']['Inlet_inlet']['random_variable_settings']["pdf_breakpoints"].GetVector()
            pdf_values = self.DEM_parameters['dem_inlets_settings']['Inlet_inlet']['random_variable_settings']["pdf_values"].GetVector()
            rv_parameters = Kratos.Parameters("""{}""")
            rv_parameters.AddVector("breakpoints", breakpoints)
            rv_parameters.AddVector("values", pdf_values)
            expected_rv = rv.PiecewiseLinearRV(rv_parameters)
            self.pdf_expected = np.array([expected_rv.InterpolateHeight(x) for x in self.centers])

            self.error = TestDEMInletAnalysis.Error(self.empirical_pdf, self.pdf_expected, interval_widths)
            TestDEMInletAnalysis.Say('Error = ' + '{:.4f}'.format(self.error), '( self.tolerance = ' + '{:.4f}'.format(self.tolerance), ')')

        super().FinalizeSolutionStep()

    def Finalize(self):
        self.samples = np.array(list(self.samples.values()))

        self.procedures.RemoveFoldersWithResults(str(self.main_path), str(self.problem_name), '')

        if debug_mode:
            import matplotlib.pyplot as plt
            plt.plot(self.centers, self.empirical_pdf, label='obtained (' + str(len(self.samples)) + ' samples)')
            plt.hist(self.samples, bins = self.n_bins, density=True, color = "skyblue", lw=0, alpha=0.25)
            plt.plot(self.centers, self.pdf_expected, label='desired')
            plt.legend()
            plt.savefig('histogram.pdf')

        super().Finalize()

class TestDEMInlet(KratosUnittest.TestCase):

    def setUp(self):
        pass

    @staticmethod
    def Say(*args):
        Logger.PrintInfo("DEM", *args)
        Logger.Flush()

    @classmethod
    def test_piecewise_linear_inlet(self):
        path = TestDEMInletAnalysis.GetMainPath()
        parameters_file_name = os.path.join(path, "ProjectParametersDEM.json")
        model = Kratos.Model()

        with open(parameters_file_name, 'r') as parameter_file:
            project_parameters = Kratos.Parameters(parameter_file.read())

        TestDEMInletAnalysis(model, project_parameters).Run()

if __name__ == "__main__":
    if debug_mode:
        Logger.GetDefaultOutput().SetSeverity(Logger.Severity.INFO)
    else:
        Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)

    KratosUnittest.main()
