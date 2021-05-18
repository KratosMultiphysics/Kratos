import os
import numpy as np
import KratosMultiphysics
from KratosMultiphysics import Logger
import KratosMultiphysics.DEMApplication as DEM
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.DEMApplication.DEM_analysis_stage

import KratosMultiphysics.kratos_utilities as kratos_utils

Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)

def Say(*args):
    Logger.PrintInfo("DEM", *args)
    Logger.Flush()

def Error(values, values_ref, interval_widths):
    return np.sqrt(sum((values_ref - values)**2 * interval_widths) / sum(interval_widths))
class TestRandomVariable(KratosUnittest.TestCase):

    def setUp(self):
        pass

    @classmethod
    def test_random_variable(self):
        settings = KratosMultiphysics.Parameters("""
        {
            "pdf_breakpoints" : [0, 1, 2, 3, 4, 5, 6],
            "pdf_values"      : [1, 2, 2, 4, 2, 5, 0],
            "do_use_seed"     : true,
            "seed"            : 1
        }
        """)

        seed = settings["seed"].GetInt()

        n_bins = 25 # for the histogram
        tolerance = 0.005
        n_max_experiments = int(1e6) # should be enough with the tolerance=0.005
        do_print_graph = True

        if settings["do_use_seed"].GetBool():
            random_variable = DEM.PiecewiseLinearRandomVariable(settings, seed)
        else:
            random_variable = DEM.PiecewiseLinearRandomVariable(settings)

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

            error = Error(empirical_pdf, pdf_expected, interval_widths)

            if i_iteration > 1:
                n_experiments *= int(max(2, (error / tolerance) ** 2))
            else:
                n_experiments *= 2
            if n_experiments > n_max_experiments:
                raise ValueError('\nThe requested tolerance (' + '{:.4f}'.format(tolerance) + ')'
                                 + ' requires a greater number of experiments than the maximum allowed (' + '{:d}'.format(n_max_experiments)
                                 + ')! Please, increase it (reduce desired precision).')

            Say('Sample size', n_experiments)
            Say('Error = ' + '{:.4f}'.format(error), '(tolerance = ' + '{:.4f}'.format(tolerance), ')')
            i_iteration += 1

        if do_print_graph:
            import matplotlib.pyplot as plt
            plt.plot(centers, empirical_pdf)
            plt.hist(sample, bins = n_bins, density=True, label='obtained (' + str(n_experiments) + ' samples)')
            plt.plot(centers, pdf_expected, label='desired')
            plt.legend()
            plt.savefig('piecewise_pdf.pdf')

    def tearDown(self):
        pass

    if __name__ == "__main__":
        KratosUnittest.main()




