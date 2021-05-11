# Importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.ShallowWaterApplication as SW

from KratosMultiphysics.ShallowWaterApplication.benchmarks.base_benchmark_process import BaseBenchmarkProcess
from KratosMultiphysics.process_factory import Factory as ProcessFactory

# Other imports
import numpy as np
from scipy.integrate import odeint

def Factory(settings, model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return MacDonaldShockBenchmark(model, settings["Parameters"])

class MacDonaldShockBenchmark(BaseBenchmarkProcess):
    """Mac Donald's shock benchmark.

    This is a Mac Donald's type solution with a smooth transition and a shock in a
    short domain, with Manning's friction coefficient.

    The length of the channel is 100m and the discharge at steady state is q=2m^2/s. The flow
    is fluvial both upstream and downstream, the boundary conditions are fixed as follows:
        - upstream: q=2m^2/s
        - downstream: h=h_ex(100)

    This process sets the upstream and downstream boundary conditions.

    O. Delestre, C. Lucas, P.-A. Ksinant, F. Darboux, C. Laguerre, T.N.T. Vo, F. James, S. Cordier
    SWASHES: a compilation of Shallow Water Analytic Solutions for Hydraulic and Environmental Studies
    International Journal for Numerical Methods in Fluids, Wiley, 2013, 72 (3), pp.269-300.
    """

    def __init__(self, model, settings):
        """Constructor of the benchmark.

        The base class validates the settings and sets the model_part, the variables and the benchmark_settings
        """

        super().__init__(model, settings)

        self.n = self.benchmark_settings["manning"].GetDouble()
        self.q = self.benchmark_settings["discharge"].GetDouble()
        self.g = self.model_part.ProcessInfo[KM.GRAVITY_Z]
        self.x0 = 0
        self.x100 = 100
        self.h0 = self._H(self.x0)
        self.h100 = self._H(self.x100)

        self.__PreComputeTopography()

    @staticmethod
    def _GetBenchmarkDefaultSettings():
        return KM.Parameters("""
            {
                "discharge"             : 2,
                "manning"               : 0.0328,
                "upstream_model_part"   : "model_part.upstream",
                "downstream_model_part" : "model_part.downstream"
            }
            """
            )

    def _Topography(self, coordinates):
        x = coordinates.X
        return self._Z(x)

    def _Height(self, coordinates, time):
        x = coordinates.X
        if time > 0:
            return self._H(x)
        else:
            return self._InitialH(x)

    def _Momentum(self, coordinates, time):
        if time > 0:
            return [self.q, 0.0, 0.0]
        else:
            return [0.0, 0.0, 0.0]

    def _Velocity(self, coordinates, time):
        return [q / self._Height(coordinates, time) for q in self._Momentum(coordinates, time)]

    def Check(self):
        """This method checks if the input values have physical sense."""

        super().Check()
        label = self.__class__.__name__
        if self.g <= 0:
            msg = label + "Gravity must be a positive value. Please, check the definition of GRAVITY_Z component in the ProcessInfo."
            raise Exception(msg)
        elif self.n < 0:
            msg = label + "The manning coefficient must be a positive value. Please, check the Parameters."
            raise Exception(msg)
        elif self.q <= 0:
            msg = label + "The discharge must be a positive value. Please, check the Parameters."
            raise Exception(msg)
        self._CheckDomain()

    def ExecuteInitialize(self):
        """This method sets the topography, the initial conditions and the upstream/downstream boundary conditions"""

        super().ExecuteInitialize()
        for process in self._GetListOfBoundaryConditionsProcesses():
            process.ExecuteInitialize()
        KM.VariableUtils().SetVariable(SW.MANNING, self.n, self.model_part.Nodes)

    def _CheckDomain(self):
        x_min = 1.0
        x_max = -1.0
        for node in self.model_part.Nodes:
            x_min = min(x_min, node.X)
            x_max = max(x_max, node.X)

        tolerance = 1e-6
        if abs(x_min - self.x0) > tolerance:
            KM.Logger.PrintWarning(self.__class__.__name__, "This benchmark expects an x-aligned model part starting at x=0")

        if abs(x_max - self.x100) > tolerance:
            KM.Logger.PrintWarning(self.__class__.__name__, "This benchmark expects an x-aligned model part ending at x=100")

    def __PreComputeTopography(self):
        X = np.linspace(self.x100, 0)
        z100 = 0
        Z = odeint(self._dZ, z100, X)
        Z = np.ndarray.flatten(Z)
        self.__X = X[::-1]
        self.__Z = Z[::-1]

    def _Z(self, x):
        return np.interp(x, self.__X, self.__Z)

    def _H1(self, x):
        g = self.g
        return (4/g)**(1/3) * (4/3 - x/100) - 9*x/1000 * (x/100 - 2/3)

    def _H2(self, x):
        g = self.g
        a1 = 0.674202
        a2 = 21.7112
        a3 = 14.492
        a4 = 1.4305
        return (4/g)**(1/3) * (a1*(x/100 - 2/3)**4 + a1*(x/100 - 2/3)**3 - a2*(x/100 - 2/3)**2 + a3*(x/100 - 2/3) + a4)

    def _dH1(self, x):
        g = self.g
        return -9*x/50000 - (4/g)**(1/3)/100 + 0.006

    def _dH2(self, x):
        g = self.g
        return (4/g)**(1/3)*(-0.00434224*x + 0.02696808*(x/100 - 0.666666666666667)**3 + 0.02022606*(x/100 - 0.666666666666667)**2 + 0.434402666666667)

    def _H(self, x):
        if x < 200/3:
            return self._H1(x)
        else:
            return self._H2(x)

    def _dH(self, x):
        if x < 200/3:
            return self._dH1(x)
        else:
            return self._dH2(x)

    def _Sf(self, h):
        return self.n**2 * self.q**2 / h**(10/3)

    def _dZ(self, z, x):
        q = self.q
        g = self.g
        return (q**2 / (g * self._H(x)**3) - 1) * self._dH(x) - self._Sf(self._H(x))

    def _InitialH(self, x):
        return np.maximum(self.h100 - self._Z(x), self.h0)

    def _GetListOfBoundaryConditionsProcesses(self):
        if not hasattr(self, 'list_of_bc_processes'):
            self.list_of_bc_processes = self._CreateListOfBoundaryConditionsProcesses()
        return self.list_of_bc_processes

    def _CreateListOfBoundaryConditionsProcesses(self):
        self.upstream_settings = KM.Parameters("""{
            "process_name" : "ApplyConstantVectorValueProcess",
            "Parameters"   : {
                "variable_name"   : "MOMENTUM",
                "is_fixed_x"      : true,
                "is_fixed_y"      : true,
                "direction"       : [1.0, 0.0, 0.0]}
        }""")
        self.upstream_settings["Parameters"].AddValue("model_part_name", self.benchmark_settings["upstream_model_part"])
        self.upstream_settings["Parameters"].AddDouble("modulus", self.q)

        self.downstream_settings = KM.Parameters("""{
            "process_name" : "ApplyConstantScalarValueProcess",
            "Parameters"   : {
                "variable_name"   : "HEIGHT",
                "is_fixed"        : true
            }
        }""")
        self.downstream_settings["Parameters"].AddValue("model_part_name", self.benchmark_settings["downstream_model_part"])
        self.downstream_settings["Parameters"].AddDouble("value", self.h100)

        list_of_bc_processes = []
        list_of_bc_processes.append(ProcessFactory(self.upstream_settings, self.model))
        list_of_bc_processes.append(ProcessFactory(self.downstream_settings, self.model))

        return list_of_bc_processes
