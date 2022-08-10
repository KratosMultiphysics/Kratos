from math import exp, sin, pi

import KratosMultiphysics as Kratos
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
import KratosMultiphysics.KratosUnittest as UnitTest

from KratosMultiphysics.kratos_utilities import DeleteFileIfExisting
from KratosMultiphysics.FluidDynamicsApplication.adjoint_stabilization_utilities import ComputeStabilizationCoefficient

@UnitTest.skipIfApplicationsNotAvailable("StatisticsApplication")
class AdjointStabilizationUtilitiesTests(UnitTest.TestCase):
    def setUp(self):
        self.write_json_output = False

    def testComputeStabilizationCoefficient(self):
        stabilization_parameters = Kratos.Parameters("""{
            "initial_coefficient_bounds"  : [0.0, 1.0],
            "tolerance"                   : 1e-3,
            "plateau_time_range"          : [5.0, 10.0],
            "plateau_max_slope"           : 0.5,
            "max_iterations"              : 20,
            "adjoint_parameters_file_name": "AdjointQSVMSSensitivity2DTest/one_element_test_adjoint_parameters.json"
        }""")

        class DummySolver:
            def __init__(self, model_part):
                self.main_model_part = model_part

        class DummyAnalysisClass:
            def __init__(self, model, parameters):
                self.model = model
                self.model_part = self.model.CreateModelPart("test")
                self.model_part.AddNodalSolutionStepVariable(Kratos.SHAPE_SENSITIVITY)
                self.model_part.CreateNewNode(1, 0, 0, 0)
                self.model_part.CreateNewNode(2, 0, 1, 0)
                self.model_part.CreateNewNode(3, 0, 0, 1)
                self.model_part.CreateNewNode(4, 0, 1, 1)
                self.model_part.CreateNewNode(5, 1, 1, 1)
                self.__solver = DummySolver(self.model_part)
                self.__stabilization_coefficient = parameters["problem_data"]["echo_level"].GetDouble()
                self.time = 1.0

            def OutputSolutionStep(self):
                for node in self.model_part.Nodes:
                    current_value = node.GetSolutionStepValue(Kratos.SHAPE_SENSITIVITY)
                    increment_rate = exp(-self.time * 2.0 * self.__stabilization_coefficient)
                    increment = Kratos.Array3([(node.Id + 1.0) * increment_rate, (node.Id * 2.0) * increment_rate, (node.Id / 3.0) * increment_rate])
                    node.SetSolutionStepValue(Kratos.SHAPE_SENSITIVITY,  current_value + increment)
                self.time += 1.0

            def Finalize(self):
                pass

            def _GetSolver(self):
                return self.__solver

        value = ComputeStabilizationCoefficient(DummyAnalysisClass, stabilization_parameters, execution_method=AdjointStabilizationUtilitiesTests.__ExecuteAnalysis)
        self.assertAlmostEqual(value, 6.860352e-01, 7)
        DeleteFileIfExisting("adjoint_stabilization_data.dat")

    def testFluidFFTUtilities(self):
        # generating signal data
        N = 1000
        delta_time = 0.001
        windowing_length = N*delta_time / 2

        x = 2.0
        amplitude_1 = (x**2 + x)
        amplitude_1_derivative = (2*x + 1)
        amplitude_2 = (x**3 + 3*x)
        amplitude_2_derivative = (3*x**2 + 3)

        time_steps = [0.0] * N
        drag_values = [0.0] * N
        drag_value_sensitivities = [0.0] * N
        for i in range(N):
            current_time = delta_time * (i+1)
            time_steps[i] = current_time

            drag_values[i] += amplitude_1 * sin(2 * pi * 100.0 * current_time)
            drag_values[i] += amplitude_2 * sin(2 * pi * 50.0 * current_time)

            drag_value_sensitivities[i] += amplitude_1_derivative * sin(2 * pi * 100.0 * current_time)
            drag_value_sensitivities[i] += amplitude_2_derivative * sin(2 * pi * 50.0 * current_time)

        # calculating from primal
        fft_utilities = KratosCFD.FluidFFTUtilities(time_steps[-1], windowing_length, delta_time)
        frequency_list, frequency_real_components, frequency_imag_components, frequency_amplitudes_square = fft_utilities.CalculateFFTFrequencyDistribution(drag_values)

        bin_index = int(frequency_list.index(50))

        self.assertAlmostEqual(frequency_amplitudes_square[bin_index], amplitude_2 ** 2, 9)

        # calculating from adjoint
        fft_utilities = KratosCFD.FluidFFTUtilities(time_steps[-1], windowing_length, delta_time)

        adjoint_drag_bin_real_value = 0.0
        adjoint_drag_bin_imag_value = 0.0
        adjoint_drag_bin_real_value_sensitivities = 0.0
        adjoint_drag_bin_imag_value_sensitivities = 0.0
        for i in range(N):
            time = time_steps[i]
            hann_windowing_coefficient = fft_utilities.CalculateHannWindowCoefficient(time)
            real_coefficient = fft_utilities.CalculateFFTRealCoefficient(bin_index, time)
            imag_coefficient = fft_utilities.CalculateFFTImagCoefficient(bin_index, time)

            adjoint_drag_bin_real_value += hann_windowing_coefficient * drag_values[i] * real_coefficient
            adjoint_drag_bin_imag_value += hann_windowing_coefficient * drag_values[i] * imag_coefficient

            adjoint_drag_bin_real_value_sensitivities += hann_windowing_coefficient * drag_value_sensitivities[i] * real_coefficient
            adjoint_drag_bin_imag_value_sensitivities += hann_windowing_coefficient * drag_value_sensitivities[i] * imag_coefficient

        self.assertAlmostEqual(frequency_real_components[bin_index], adjoint_drag_bin_real_value, 9)
        self.assertAlmostEqual(frequency_imag_components[bin_index], adjoint_drag_bin_imag_value, 9)

        adjoint_drag_amplitude_square_value = fft_utilities.CalculateFFTAmplitudeSquare(adjoint_drag_bin_real_value, adjoint_drag_bin_imag_value)
        self.assertAlmostEqual(adjoint_drag_amplitude_square_value, amplitude_2 ** 2, 9)

        adjoint_drag_amplitude_square_value_sensitivities = fft_utilities.CalculateFFTAmplitudeSquareDerivative(
            adjoint_drag_bin_real_value,
            adjoint_drag_bin_real_value_sensitivities,
            adjoint_drag_bin_imag_value,
            adjoint_drag_bin_imag_value_sensitivities)

        self.assertAlmostEqual(adjoint_drag_amplitude_square_value_sensitivities, 2 * amplitude_2 * amplitude_2_derivative, 9)


    @staticmethod
    def __ExecuteAnalysis(analysis_class_type, adjoint_parameters, stabilization_coefficient, solve_id):
        model = Kratos.Model()
        adjoint_parameters["problem_data"]["echo_level"].SetDouble(stabilization_coefficient)
        analysis = analysis_class_type(model, adjoint_parameters)
        for i in range(0, 10):
            analysis.OutputSolutionStep()
        analysis.Finalize()
        return analysis.time_series_data



if __name__ == '__main__':
    UnitTest.main()