import os
import numpy as np
import random
import KratosMultiphysics
from KratosMultiphysics import Logger
import KratosMultiphysics.DEMApplication as DEM
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.DEMApplication.DEM_analysis_stage

import KratosMultiphysics.kratos_utilities as kratos_utils
import auxiliary_functions_for_tests

Logger.GetDefaultOutput().SetSeverity(Logger.Severity.DETAIL)


def Say(*args):
    Logger.PrintInfo("SwimmingDEM", *args)
    Logger.Flush()

this_working_dir_backup = os.getcwd()

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)


class Parameters:
    def __init__(self, *args, **kwargs):
        pass

class RandomVariable:
    def __init__(self, parameters):
        pass

    def Sample(self,  n_to_pick=1):
        return 0.0

def Gaussian(x, mu=0, sigma=1):
    return 1/(sigma * np.sqrt(2 * np.pi)) * np.exp( - (x - mu)**2 / (2 * sigma**2))

def Error(values, values_ref, interval_widths):
    return np.sqrt(sum((values_ref - values)**2 * interval_widths) / sum(interval_widths))

class PiecewiseLinearRV(RandomVariable):
    def __init__(self, parameters):
        self.boundaries = np.array(parameters.boundaries, dtype=np.float32)
        self.heights = np.array(parameters.heights, dtype=np.float32)
        self.ranges = np.array([(self.boundaries[i+1] - self.boundaries[i] for i in range(len(self.boundaries)))])
        self.areas = np.zeros(len(self.heights) - 1)
        self.trapezoid_indices = range(len(self.areas))
        self.Normalize()
        super().__init__(parameters)

    def InterpolateHeight(self, x):
        if x < self.boundaries[0] or x > self.boundaries[-1]:
            raise Exception('x =' + str(x) + 'is out of bounds (' + str([self.boundaries[0], self.boundaries[-1]]) + ')')

        for i, x_boundary in enumerate(self.boundaries):
            if x < x_boundary:
                break

        b1 = self.heights[i - 1]
        b2 = self.heights[i]
        h = self.boundaries[i] - self.boundaries[i - 1]
        x1 = self.boundaries[i - 1]
        x2 = self.boundaries[i]
        return (b1 * (x2 - x) + b2 * (x - x1)) / h

    def Normalize(self):
        self.total_area = 0.0
        for i in range(len(self.heights[:-1])):
            area = 0.5 * (self.boundaries[i+1] - self.boundaries[i]) * (self.heights[i+1] + self.heights[i])
            self.areas[i] = area
            self.total_area += area
        self.heights /= self.total_area
        self.discrete_probabilities = self.areas / self.total_area
        print('area = ', self.total_area)
        print('trapezoid_indices: ', self.trapezoid_indices)
        print('discrete_probabilities = ', self.discrete_probabilities)

    def Sample(self):
        i_trapezoid = self.SampleTrapezoidChoice()[0]
        x0 = self.boundaries[i_trapezoid]
        H = self.boundaries[i_trapezoid + 1] - x0
        B1 = self.heights[i_trapezoid]
        B2 = self.heights[i_trapezoid + 1]
        x_within = self.SampleWithinTrapezoid(H, B1, B2)
        return x0 + x_within

    def SampleTrapezoidChoice(self, n_to_pick=1):
        draw = np.random.choice(len(self.areas), n_to_pick, p=self.discrete_probabilities)
        return draw

    def SampleWithinTrapezoid(self, H, B1, B2):
        if B1 == 0:
            x = self.SamplePositiveSlopingStandardTriangle()
        else:
            beta = B2/B1
            b = 2.0 / (1 + beta)
            x = self.SampleXWithinStandardTrapezoid(b)
        return H * x

    def SampleXWithinStandardTrapezoid(self, b):
        alpha = random.uniform(0, 1)
        if alpha < b/2:
            y = self.SampleNegativeSlopingStandardTriangle()
        else:
            y = self.SamplePositiveSlopingStandardTriangle()
        return y

    def SamplePositiveSlopingStandardTriangle(self):
        return np.sqrt(random.uniform(0, 1))

    def SampleNegativeSlopingStandardTriangle(self):
        return 1.0 - self.SamplePositiveSlopingStandardTriangle()

class TestRandomVariable(KratosUnittest.TestCase):

    def setUp(self):
        pass

    @classmethod
    def test_random_variable(self):
        params = Parameters()
        settings = KratosMultiphysics.Parameters("""
        {
            "pdf_breakpoints" : [0, 1, 2, 3, 4, 5, 6],
            "pdf_values"      : [1, 2, 2, 4, 2, 5, 0]
        }
        """)
        params.boundaries = settings['pdf_breakpoints'].GetVector()
        params.heights = settings['pdf_values'].GetVector()

        rvar = PiecewiseLinearRV(params) # python alternative
        random_variable = DEM.PiecewiseLinearRandomVariable(settings, 1)

        # print('sample trapezoid:', rvar.SampleTrapezoid())

        n_bins = 25
        tolerance = 0.005
        max_experiments = int(1e6)

        import matplotlib.pyplot as plt
        n_experiments = 2 * n_bins
        e = tolerance + 1
        i_iteration = 1

        while e > tolerance:
            pdf = np.zeros(n_experiments)
            piecewise_linear_pdf = np.zeros(n_experiments)

            for i in range(n_experiments):
                # value = rvar.Sample()
                value = random_variable.Sample()
                piecewise_linear_pdf[i] = value

            y, edges = np.histogram(piecewise_linear_pdf, n_bins, density=True)
            centers = 0.5 * (edges[1:] + edges[:-1])
            widths = edges[1:] - edges[:-1]
            pdf_expected = np.array([rvar.InterpolateHeight(x) for x in centers])

            e = Error(y, pdf_expected, widths)

            if i_iteration > 1:
                n_experiments *= int(max(2, (e / tolerance) ** 2))
            else:
                n_experiments *= 2
            if n_experiments > max_experiments:
                raise ValueError('\nThe requested tolerance (' + '{:.4f}'.format(tolerance) + ')'
                                 + ' requires a greater number of experiments than the maximum allowed (' + '{:d}'.format(max_experiments)
                                 + ')! Please, increase it (reduce desired precision).')

            Say('Sample size', n_experiments)
            Say('Error = ' + '{:.4f}'.format(e), '(tolerance = ' + '{:.4f}'.format(tolerance), ')')
            i_iteration += 1

        plt.plot(centers, y)
        plt.hist(piecewise_linear_pdf, bins = n_bins, density=True, label='obtained (' + str(n_experiments) + ' samples)')
        plt.plot(centers, pdf_expected, label='desired')
        plt.legend()
        plt.savefig('piecewise_pdf.pdf')

    def tearDown(self):
        pass


    if __name__ == "__main__":
        KratosUnittest.main()




