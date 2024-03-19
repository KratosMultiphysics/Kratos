from math import sin, pi

import KratosMultiphysics.OptimizationApplication as KratosOA
import KratosMultiphysics.KratosUnittest as UnitTest

class TestFFTUtils(UnitTest.TestCase):
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
        fft_utilities = KratosOA.FFTUtils(time_steps[-1], windowing_length, delta_time)
        frequency_list, frequency_real_components, frequency_imag_components, frequency_amplitudes_square = fft_utilities.CalculateFFTFrequencyDistribution(drag_values)

        bin_index = int(frequency_list.index(50))

        self.assertAlmostEqual(frequency_amplitudes_square[bin_index], amplitude_2 ** 2, 9)

        # calculating from adjoint
        fft_utilities = KratosOA.FFTUtils(time_steps[-1], windowing_length, delta_time)

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

if __name__ == '__main__':
    UnitTest.main()