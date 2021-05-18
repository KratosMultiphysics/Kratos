# Importing the Kratos Library
import KratosMultiphysics as KM

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_convergence_accelerator import CoSimulationConvergenceAccelerator

# CoSimulation imports
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools

def Create(settings):
    cs_tools.SettingsTypeCheck(settings)
    return ConstantRelaxationConvergenceAccelerator(settings)

class ConstantRelaxationConvergenceAccelerator(CoSimulationConvergenceAccelerator):
    ## The constructor.
    # @param alpha relaxation factor.
    def __init__(self, settings):
        super().__init__(settings)
        self.alpha = self.settings["alpha"].GetDouble()

    ## UpdateSolution(r, x)
    # @param r residual r_k
    # @param x solution x_k
    # Computes the approximated update in each iteration.
    def UpdateSolution( self, r, x ):
        if self.echo_level > 3:
            cs_tools.cs_print_info(self._ClassName(), "Doing relaxation with factor = ", "{0:.1g}".format(self.alpha))
        delta_x = self.alpha * r
        return delta_x

    @classmethod
    def SupportsDistributedData(cls):
        return True

    @classmethod
    def _GetDefaultParameters(cls):
        this_defaults = KM.Parameters("""{
            "alpha" : 0.125
        }""")
        this_defaults.AddMissingParameters(super()._GetDefaultParameters())
        return this_defaults

