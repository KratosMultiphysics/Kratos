# Importing the Kratos Library
import KratosMultiphysics as KM

from KratosMultiphysics.ShallowWaterApplication.benchmarks.base_benchmark_process import BaseBenchmarkProcess

# Other imports
import numpy as np

def Factory(settings, model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return PlanarSurfaceInParabolaBenchmark(model, settings["Parameters"])

class PlanarSurfaceInParabolaBenchmark(BaseBenchmarkProcess):
    def __init__(self, model, settings):
        # The base class sets the model_part, variables and benchmark_settings
        super(PlanarSurfaceInParabolaBenchmark, self).__init__(model, settings)

        benchmark_default_settings = KM.Parameters("""
            {
                "depth"     : 1.0,
                "amplitude" : 1.0
            }
            """
            )
        self.benchmark_settings.ValidateAndAssignDefaults(benchmark_default_settings)

        self.h0 = self.benchmark_settings["depth"].GetDouble()
        self.a = self.benchmark_settings["amplitude"].GetDouble()
        self.g = self.model_part.ProcessInfo[KM.GRAVITY_Z]
        self.B = self.__B()
        self.C = self.__C()
        self.L = self.__L()

    def Topography(self, coordinates):
        x = coordinates.X
        return self.h0 * (1/self.a**2 * (x - 0.5*self.L)**2 - 1.0)

    def Height(self, coordinates, time):
        x = coordinates.X
        x0 = self.__x0(time)
        x1 = self.__x1(time)
        if x0 < x < x1:
            return -self.h0*(((x - 0.5*self.L)/self.a + 0.5/self.a*np.cos(self.C*time/self.a))**2 - 1)
        else:
            return 0.0

    def Velocity(self, coordinates, time):
        x = coordinates.X
        x0 = self.__x0(time)
        x1 = self.__x1(time)
        if x0 < x < x1:
            return [self.B * np.sin(self.C*time/self.a), 0.0, 0.0]
        else:
            return [0.0, 0.0, 0.0]

    def Check(self):
        super(PlanarSurfaceInParabolaBenchmark, self).Check()
        label = self.__class__.__name__
        if self.g <= 0:
            msg = label + "Gravity must be a positive value. Please, check the definition of GRAVITY_Z compontent in the ProcessInfo."
            raise Exception(msg)
        elif self.L <= 0:
            msg = label + "The length must be a positive value. Please, check the Parameters."
            raise Exception(msg)
        elif self.h0 <= 0:
            msg = label + "The depth must be a positive value. Please, check the Parameters."
            raise Exception(msg)
        elif self.a <= 0:
            msg = label + "The amplitude must be a positive value. Please, check the Parameters."
            raise Exception(msg)

    def __B(self):
        return self.__C() / 2.0 / self.a

    def __C(self):
        return np.sqrt(2*self.g*self.h0)

    def __L(self):
        x0 = 1.0
        x1 = -1.0
        for node in self.model_part.Nodes:
            x0 = min(x0, node.X)
            x1 = max(x1, node.X)

        tolerance = 1e-6
        if abs(x0) > tolerance:
            KM.Logger.PrintWarning(self.__class__.__name__, "This benchmark expects an x-aligned model part starting at x=0")

        if x1 <= 0.0:
            KM.Logger.PrintWarning(self.__class__.__name__, "This benchmark expects a model part with x>0")

        return x1 - x0

    def __x0(self, time):
        return -0.5*np.cos(self.C/self.a*time) - self.a + 0.5*self.L

    def __x1(self, time):
        return -0.5*np.cos(self.C/self.a*time) + self.a + 0.5*self.L
