import os
import numpy as np
import KratosMultiphysics as Kratos
from Kratos import Logger
Logger.GetDefaultOutput().SetSeverity(Logger.Severity.DETAIL)
import KratosMultiphysics.DEMApplication as DEM
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.DEMApplication.DEM_analysis_stage

import random_variable_tests_files.random_variable as rv

this_working_dir_backup = os.getcwd()

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class TestDEMInlet(Kratos.DEMApplication.DEM_analysis_stage.DEMAnalysisStage, KratosUnittest.TestCase):
    debug_mode = True

    @staticmethod
    def Say(*args):
        Logger.PrintInfo("DEM", *args)
        Logger.Flush()

    @staticmethod
    def Error(values, values_ref, interval_widths):
        return sum(abs((values_ref - values)) * interval_widths)

    def __init__(self, model, DEM_parameters):
        self.tolerance = 0.05
        self.error = self.tolerance + 1
        self.n_bins = 25 # for the histogram
        self.mu = 0.001 / 2
        self.sigma = 0.00025
        self.centers = []
        self.empirical_pdf = []
        self.pdf_expected = []
        super().__init__(model, DEM_parameters)

    def setUp(self):
        if TestDEM3DContact.debug_mode:
            Logger.GetDefaultOutput().SetSeverity(Logger.Severity.DETAIL)
        else:
            Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)

    def Initialize(self):
        self.samples = dict()
        return super().Initialize()

    def InitializeSolutionStep(self):
        return super().InitializeSolutionStep()

    def HasConverged(self):
        return self.error < self.tolerance

    def KeepAdvancingSolutionLoop(self):
        return super().KeepAdvancingSolutionLoop() and not self.HasConverged()

    @classmethod
    def GetMainPath(self):
        return os.path.join(os.path.dirname(os.path.realpath(__file__)), "piecewise_linear_inlet_tests_files")

    def GetProblemNameWithPath(self):
        return os.path.join(self.main_path, self.DEM_parameters["problem_name"].GetString())

    def FinalizeSolutionStep(self):
        particle_nodes = [node for node in self.spheres_model_part.Nodes if node.IsNot(Kratos.BLOCKED)]

        for node in particle_nodes:
            self.samples[node.Id] = node.GetSolutionStepValue(Kratos.RADIUS)

        if particle_nodes:
            self.empirical_pdf, interval_boundaries = np.histogram(list(self.samples.values()), self.n_bins, density=True)
            self.centers = 0.5 * (interval_boundaries[1:] + interval_boundaries[:-1])
            interval_widths = interval_boundaries[1:] - interval_boundaries[:-1]
            range_x = np.linspace(interval_boundaries[0], interval_boundaries[-1], num=self.n_bins)
            parameters = Kratos.Parameters("""{}""")
            breakpoints = Kratos.Vector(list(range_x))
            pdf_values = Kratos.Vector(list(rv.GaussianPDF(np.array(range_x), mu=self.mu, sigma=self.sigma)))
            parameters.AddVector("breakpoints", breakpoints)
            parameters.AddVector("values", pdf_values)
            expected_rv = rv.PiecewiseLinearRV(parameters)
            self.pdf_expected = np.array([expected_rv.InterpolateHeight(x) for x in self.centers])

            self.error = TestDEMInlet.Error(self.empirical_pdf, self.pdf_expected, interval_widths)
            TestDEMInlet.Say('Probability error:', self.error)
        super().FinalizeSolutionStep()

    def Finalize(self):
        self.samples = np.array(list(self.samples.values()))
        TestDEMInlet.Say('samples:', self.samples)

        self.procedures.RemoveFoldersWithResults(str(self.main_path), str(self.problem_name), '')

        TestDEMInlet.Say('Error = ' + '{:.4f}'.format(self.error), '( self.tolerance = ' + '{:.4f}'.format(self.tolerance), ')')

        if TestDEMInlet.debug_mode:
            import matplotlib.pyplot as plt
            plt.plot(self.centers, self.empirical_pdf)
            plt.hist(self.samples, bins = self.n_bins, density=True, label='obtained (' + str(len(self.samples)) + ' samples)')
            plt.plot(self.centers, self.pdf_expected, label='desired')
            plt.legend()
            plt.savefig('piecewise_pdf.pdf')
        super().Finalize()


class TestDEM3DContact(KratosUnittest.TestCase):
    debug_mode = False

    def setUp(self):
        pass

    @staticmethod
    def Error(values, values_ref, interval_widths):
        return np.sqrt(sum((values_ref - values)**2 * interval_widths) / sum(interval_widths))

    @staticmethod
    def Say(*args):
        Logger.PrintInfo("DEM", *args)
        Logger.Flush()


    @classmethod
    def test_piecewise_linear_inlet(self):
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "piecewise_linear_inlet_tests_files")
        parameters_file_name = os.path.join(path, "ProjectParametersDEM.json")
        model = Kratos.Model()

        # Test parallel computation.
        with open(parameters_file_name,'r') as parameter_file:
            project_parameters = Kratos.Parameters(parameter_file.read())
        TestDEMInlet(model, project_parameters).Run()


if __name__ == "__main__":
    Logger.GetDefaultOutput().SetSeverity(Logger.Severity.DETAIL)
    KratosUnittest.main()
