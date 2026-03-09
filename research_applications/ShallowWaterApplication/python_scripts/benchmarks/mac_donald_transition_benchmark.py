# Importing the Kratos Library
import KratosMultiphysics as KM

from KratosMultiphysics.ShallowWaterApplication.benchmarks.mac_donald_shock_benchmark import MacDonaldShockBenchmark
from KratosMultiphysics.process_factory import Factory as ProcessFactory

def Factory(settings, model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return MacDonaldTransitionBenchmark(model, settings["Parameters"])

class MacDonaldTransitionBenchmark(MacDonaldShockBenchmark):
    """Mac Donald's transcritical flow.

    This is a Mac Donald's type solution with a smooth transition to supercritical flow in a
    short domain, with Manning's friction coefficient.

    The length of the channel is 100m and the discharge at steady state is q=2m^2/s. The flow
    is fluvial upstream and torrential downstream, the boundary conditions are fixed as follows:
        - upstream: q=2m^2/s
        - downstream: free

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

    @staticmethod
    def _GetBenchmarkDefaultSettings():
        return KM.Parameters("""
            {
                "discharge"             : 2,
                "manning"               : 0.0328,
                "upstream_model_part"   : "model_part.upstream"
            }
            """)

    def _H(self, x):
        g = self.g
        return (4/g)**(1/3) * (1 - (x-50)/200 + (x-50)**2/30000)

    def _dH(self, x):
        g = self.g
        return (4/g)**(1/3) * (x/15000 - 5/600)

    def _InitialH(self, x):
        if x > self.x0:
            return 1e-6
        else:
            return self.h0

    def _CreateListOfBoundaryConditionsProcesses(self):
        benchmark_settings = self.settings["benchmark_settings"]

        self.upstream_settings = KM.Parameters("""{
            "process_name" : "ApplyConstantVectorValueProcess",
            "Parameters"   : {
                "variable_name"   : "MOMENTUM",
                "is_fixed_x"      : true,
                "is_fixed_y"      : true,
                "direction"       : [1.0, 0.0, 0.0]}
        }""")
        self.upstream_settings["Parameters"].AddValue("model_part_name", benchmark_settings["upstream_model_part"])
        self.upstream_settings["Parameters"].AddDouble("modulus", self.q)

        list_of_bc_processes = []
        list_of_bc_processes.append(ProcessFactory(self.upstream_settings, self.model))

        return list_of_bc_processes
