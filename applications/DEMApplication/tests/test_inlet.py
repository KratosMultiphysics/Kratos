import os
import numpy as np
import KratosMultiphysics
from KratosMultiphysics import Logger
Logger.GetDefaultOutput().SetSeverity(Logger.Severity.DETAIL)
import KratosMultiphysics.DEMApplication as DEM
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.DEMApplication.DEM_analysis_stage

import auxiliary_functions_for_tests

this_working_dir_backup = os.getcwd()

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

def Gaussian(x, mu=0, sigma=1):
	return 1/(sigma * np.sqrt(2 * np.pi)) * np.exp( - (x - mu)**2 / (2 * sigma**2))


class TestDEMInlet(KratosMultiphysics.DEMApplication.DEM_analysis_stage.DEMAnalysisStage, KratosUnittest.TestCase):
    debug_mode = True

    @staticmethod
    def Say(*args):
        Logger.PrintInfo("DEM", *args)
        Logger.Flush()

    @staticmethod
    def Error(values, values_ref, interval_widths):
        return np.sqrt(sum((values_ref - values)**2 * interval_widths) / sum(interval_widths))


    def setUp(self):
        if TestDEM3DContact.debug_mode:
            Logger.GetDefaultOutput().SetSeverity(Logger.Severity.DETAIL)
        else:
            Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)

    def InitializeSolutionStep(self):
        self.n_bins = 25 # for the histogram
        self.tolerance = 0.005
        self.n_max_experiments = int(1e6) # should be enough with the self.tolerance=0.005
        self.mu = 0.0001 / 2
        self.sigma = 0.000025

        return super().InitializeSolutionStep()

    @classmethod
    def GetMainPath(self):
        return os.path.join(os.path.dirname(os.path.realpath(__file__)), "piecewise_linear_inlet_tests_files")

    def GetProblemNameWithPath(self):
        return os.path.join(self.main_path, self.DEM_parameters["problem_name"].GetString())

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        pass

    def Finalize(self):
        particle_nodes = [node for node in self.spheres_model_part.Nodes if node.IsNot(KratosMultiphysics.BLOCKED)]
        self.sample = np.zeros(len(particle_nodes))

        for i, node in enumerate(particle_nodes):
            self.sample[i] = node.GetSolutionStepValue(KratosMultiphysics.RADIUS)
        TestDEMInlet.Say(self.sample)

        self.procedures.RemoveFoldersWithResults(str(self.main_path), str(self.problem_name), '')
        empirical_pdf, interval_boundaries = np.histogram(self.sample, self.n_bins, density=True)
        centers = 0.5 * (interval_boundaries[1:] + interval_boundaries[:-1])
        interval_widths = interval_boundaries[1:] - interval_boundaries[:-1]
        range_x = np.linspace(interval_boundaries[0], interval_boundaries[-1], num=self.n_bins)
        pdf_expected = np.array([Gaussian(x, self.mu, self.sigma) for x in centers])

        error = TestDEMInlet.Error(empirical_pdf, pdf_expected, interval_widths)

        TestDEMInlet.Say('Error = ' + '{:.4f}'.format(error), '( self.tolerance = ' + '{:.4f}'.format(self.tolerance), ')')

        if TestDEMInlet.debug_mode:
            import matplotlib.pyplot as plt
            plt.plot(centers, empirical_pdf)
            plt.hist(self.sample, bins = self.n_bins, density=True, label='obtained (' + str(i) + ' samples)')
            plt.plot(centers, pdf_expected, label='desired')
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
        model = KratosMultiphysics.Model()

        # Test parallel computation.
        with open(parameters_file_name,'r') as parameter_file:
            project_parameters = KratosMultiphysics.Parameters(parameter_file.read())
        TestDEMInlet(model, project_parameters).Run()


if __name__ == "__main__":
    Logger.GetDefaultOutput().SetSeverity(Logger.Severity.DETAIL)
    KratosUnittest.main()
