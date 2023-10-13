import os
import numpy as np
import KratosMultiphysics as Kratos
from Kratos import Logger
import KratosMultiphysics.KratosUnittest as KratosUnittest

import KratosMultiphysics.DEMApplication.DEM_analysis_stage as dem_analysis

import random_variable_tests_files.random_variable as rv

debug_mode = False

class TestDEMPiecewiseInletAnalysis(dem_analysis.DEMAnalysisStage):

    @staticmethod
    def Say(*args):
        Logger.PrintInfo("DEM", *args)
        Logger.Flush()

    @staticmethod # This error is bounded by 1
    def Error(values, values_ref, interval_widths):
        return 0.5 * sum(abs((values_ref - values)) * interval_widths)

    @staticmethod
    def StaticGetMainPath():
        return os.path.join(os.path.dirname(os.path.realpath(__file__)), "random_variable_based_inlet_tests_files")

    def GetMainPath(self):
        return self.StaticGetMainPath()

    def __init__(self, model, DEM_parameters):
        self.tolerance = 0.05
        self.error_norm = self.tolerance + 1
        self.n_bins = 50 # for the histogram

        super().__init__(model, DEM_parameters)

    def Initialize(self):
        self.samples = dict()
        return super().Initialize()

    def InitializeSolutionStep(self):
        return super().InitializeSolutionStep()

    def HasConverged(self):
        return self.error_norm < self.tolerance

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

            self.error_norm = TestDEMPiecewiseInletAnalysis.Error(self.empirical_pdf, self.pdf_expected, interval_widths)
            TestDEMPiecewiseInletAnalysis.Say('Error = ' + '{:.4f}'.format(self.error_norm), '( self.tolerance = ' + '{:.4f}'.format(self.tolerance), ')')


            inlet_sub_model_part = self.dem_inlet_model_part.GetSubModelPart('Inlet_inlet')
            self.max_radius = self.DEM_inlet.GetMaxRadius(inlet_sub_model_part)
            self.max_radius_expected = breakpoints[np.where(pdf_values)[0].max()]

        super().FinalizeSolutionStep()

    def Finalize(self):
        self.samples = np.array(list(self.samples.values()))

        self.procedures.RemoveFoldersWithResults(str(self.main_path), str(self.problem_name), '')

        if debug_mode:
            import matplotlib.pyplot as plt
            plt.plot(self.centers, self.pdf_expected, label='desired')
            plt.plot(self.centers, self.empirical_pdf, label='obtained (' + str(len(self.samples)) + ' samples)')
            plt.hist(self.samples, bins = self.n_bins, density=True, color = "skyblue", lw=0, alpha=0.25)
            plt.xlabel("radius")
            plt.ylabel("probability density")
            plt.legend()
            plt.savefig('histogram_' + self.problem_name + '.pdf')
            plt.close()

        super().Finalize()

class TestDEMDiscreteInletAnalysis(TestDEMPiecewiseInletAnalysis):

    def __init__(self, model, DEM_parameters):
        super().__init__(model, DEM_parameters)
        rv_settings = DEM_parameters['dem_inlets_settings']['Inlet_inlet']['random_variable_settings']
        self.possible_values = np.array(rv_settings["possible_values"].GetVector())
        self.probabilities = np.array(rv_settings["relative_frequencies"].GetVector())
        self.relative_closeness_tolerance = rv_settings["relative_closeness_tolerance"].GetDouble()
        self.histogram_expected = self.probabilities / sum(self.probabilities)
        self.histogram = np.zeros(len(self.possible_values))

    def AreClose(self, a, b):
        return (a < b + self.relative_closeness_tolerance
            and a > b - self.relative_closeness_tolerance)

    def FinalizeSolutionStep(self):
        particle_nodes = [node for node in self.spheres_model_part.Nodes if node.IsNot(Kratos.BLOCKED)]

        for node in particle_nodes:
            sample = node.GetSolutionStepValue(Kratos.RADIUS)
            not_found = True
            for i, b in enumerate(list(self.possible_values)):

                if self.AreClose(sample, b):
                    # Mapping to Id so it works even if some particles are removed before the end of the simulation
                    self.samples[node.Id] = sample
                    self.histogram[i] += 1
                    not_found = False
                    break

            if not_found:
                raise ValueError('The DEM inlet has generated a particle with a radius that is not within' +
                                 'the chosen tolerance (' + '{:.2e}'.format(self.relative_closeness_tolerance) + ') of any' +
                                 'of the possible values in the imposed discrete distribution.')

        if sum(self.histogram):
            error = self.histogram_expected - self.histogram / sum(self.histogram)
            self.error_norm = sum(abs(p) for p in error)

            inlet_sub_model_part = self.dem_inlet_model_part.GetSubModelPart('Inlet_inlet')
            self.max_radius = self.DEM_inlet.GetMaxRadius(inlet_sub_model_part)
            self.max_radius_expected = self.possible_values[np.where(self.probabilities)[0].max()]


        TestDEMPiecewiseInletAnalysis.Say('Error = ' + '{:.2e}'.format(self.error_norm) + ' (Error tolerance = ' + '{:.2e}'.format(self.tolerance) + ')')

        dem_analysis.DEMAnalysisStage.FinalizeSolutionStep(self)

    def Finalize(self):
        self.samples = np.array(list(self.samples.values()))

        self.procedures.RemoveFoldersWithResults(str(self.main_path), str(self.problem_name), '')

        if debug_mode:
            import matplotlib.pyplot as plt
            normalization_factor = sum(self.histogram) / sum(self.probabilities)
            support_diameter = self.possible_values[-1] - self.possible_values[0]

            min_distance_between_possible_values = min(self.possible_values[i] - self.possible_values[i-1]
                                                       for i in range(1, len(self.possible_values)))

            column_width = 0.4 * max(min_distance_between_possible_values, 0.05 * support_diameter / len(self.possible_values))
            plt.bar(self.possible_values - 0.5 * column_width, self.probabilities * normalization_factor, width=column_width, label='desired')
            plt.bar(self.possible_values + 0.5 * column_width, self.histogram, width=column_width,
                    label='obtained' + ' (error = ' + '{:.2e}'.format(self.error_norm) + ')')
            plt.xlabel('radius')
            plt.ylabel('quantity')
            plt.legend()

            plt.savefig('histogram_' + self.problem_name + '.pdf')
            plt.close()

        dem_analysis.DEMAnalysisStage.Finalize(self)

class TestPieceWiseLinearDEMInlet(KratosUnittest.TestCase):
    def setUp(self):
        self.parameters_file_name = 'PiecewiseLinearProjectParametersDEM.json'
        self.path = TestDEMPiecewiseInletAnalysis.StaticGetMainPath()
        self.analysis = TestDEMPiecewiseInletAnalysis

    def test_piecewise_linear_inlet(self):
        parameters_file_name = os.path.join(self.path, self.parameters_file_name)
        model = Kratos.Model()

        with open(parameters_file_name, 'r') as parameter_file:
            project_parameters = Kratos.Parameters(parameter_file.read())

        analysis = self.analysis(model, project_parameters)
        analysis.Run()
        self.assertAlmostEqual(analysis.max_radius_expected, analysis.max_radius)


class TestDiscreteDEMInlet(TestPieceWiseLinearDEMInlet):
    def setUp(self):
        self.parameters_file_name = 'DiscreteProjectParametersDEM.json'
        self.path = TestDEMDiscreteInletAnalysis.StaticGetMainPath()
        self.analysis = TestDEMDiscreteInletAnalysis


if __name__ == "__main__":
    if debug_mode:
        Logger.GetDefaultOutput().SetSeverity(Logger.Severity.INFO)
    else:
        Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)

    KratosUnittest.main()
