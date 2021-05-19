import os
import numpy as np
import KratosMultiphysics
from KratosMultiphysics import Logger
import KratosMultiphysics.DEMApplication as DEM
import KratosMultiphysics.KratosUnittest as KratosUnittest

class TestRandomVariable(KratosUnittest.TestCase):
    debug_mode = False

    @staticmethod
    def Error(values, values_ref, interval_widths):
        return np.sqrt(sum((values_ref - values)**2 * interval_widths) / sum(interval_widths))

    @staticmethod
    def Say(*args):
        Logger.PrintInfo("DEM", *args)
        Logger.Flush()

    @staticmethod
    def GetMean(pdf_breakpoints, pdf_values):
        mean = 0.0
        total_area = 0.0
        for i in range(len(pdf_breakpoints[:-1])):
            x1 = pdf_breakpoints[i]
            x2 = pdf_breakpoints[i+1]
            v1 = pdf_values[i]
            v2 = pdf_values[i+1]
            h = x2 - x1
            trapezoid_area = 0.5 * (v1 + v2) * h
            min_v = min(v1, v2)
            max_v = max(v1, v2)
            square_sub_area = min_v * h
            triangle_sub_area = 0.5 * (max_v - min_v)  * h
            delta_v = v2 - v1
            x_triangle = (0.5 + 1.0/6 * np.sign(delta_v)) * h
            x_square = 0.5 * h
            x_centroid_within_trapezoid = (triangle_sub_area * x_triangle + square_sub_area * x_square) / trapezoid_area
            mean += (x1 + x_centroid_within_trapezoid) * trapezoid_area
            total_area += trapezoid_area
        return mean / total_area

    def setUp(self):
        if TestRandomVariable.debug_mode:
            Logger.GetDefaultOutput().SetSeverity(Logger.Severity.DETAIL)
        else:
            Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)

    def test_random_variable(self):
        settings = KratosMultiphysics.Parameters("""
        {
            "pdf_breakpoints" : [0, 1, 2, 3, 3.5, 4, 4.5, 5, 5.5, 6],
            "pdf_values"      : [1, 2, 3, 4, 3.5, 3, 2.5, 2, 1.5, 1],
            "do_use_seed"     : true,
            "seed"            : 1,
            "relative_closeness_tolerance" : 1e-6
        }
        """)

        seed = settings["seed"].GetInt()

        n_bins = 25 # for the histogram
        tolerance = 0.005
        n_max_experiments = int(1e6) # should be enough with the tolerance=0.005

        if settings["do_use_seed"].GetBool():
            random_variable = DEM.PiecewiseLinearRandomVariable(settings, seed)
        else:
            random_variable = DEM.PiecewiseLinearRandomVariable(settings)
        mean_expected = TestRandomVariable.GetMean(settings["pdf_breakpoints"].GetVector(), settings["pdf_values"].GetVector())
        mean = random_variable.GetMean()
        self.assertAlmostEqual(mean_expected, mean, None, 'The calculated PDF mean differs from the expected one.')

        TestRandomVariable.Say('Mean:', mean)

        n_experiments = 2 * n_bins
        error = tolerance + 1
        i_iteration = 1

        # Increase number of experiments until the 'L2-norm' of the
        # error of the recovered pdf is smaller than the tolerance
        while error > tolerance:
            sample = np.zeros(n_experiments)

            for i in range(n_experiments):
                sample[i] = random_variable.Sample()

            empirical_pdf, interval_boundaries = np.histogram(sample, n_bins, density=True)
            centers = 0.5 * (interval_boundaries[1:] + interval_boundaries[:-1])
            interval_widths = interval_boundaries[1:] - interval_boundaries[:-1]
            pdf_expected = np.array([random_variable.ProbabilityDensity(x) for x in centers])

            error = TestRandomVariable.Error(empirical_pdf, pdf_expected, interval_widths)

            if i_iteration > 1:
                n_experiments *= int(max(2, (error / tolerance) ** 2))
            else:
                n_experiments *= 2
            if n_experiments > n_max_experiments:
                raise ValueError('\nThe requested tolerance (' + '{:.4f}'.format(tolerance) + ')'
                                + ' requires a greater number of experiments than the maximum allowed (' + '{:d}'.format(n_max_experiments)
                                + ')! Please, increase it (reduce desired precision).')

            TestRandomVariable.Say('Sample size', n_experiments)
            TestRandomVariable.Say('Error = ' + '{:.4f}'.format(error), '( tolerance = ' + '{:.4f}'.format(tolerance), ')')
            i_iteration += 1

        if TestRandomVariable.debug_mode:
            import matplotlib.pyplot as plt
            plt.plot(centers, empirical_pdf)
            plt.hist(sample, bins = n_bins, density=True, label='obtained (' + str(n_experiments) + ' samples)')
            plt.plot(centers, pdf_expected, label='desired')
            plt.legend()
            plt.savefig('piecewise_pdf.pdf')

    def tearDown(self):
        pass

if __name__ == "__main__":
    Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
    KratosUnittest.main()
