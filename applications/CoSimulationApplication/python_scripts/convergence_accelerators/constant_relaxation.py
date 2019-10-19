from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_convergence_accelerator import CoSimulationConvergenceAccelerator

# CoSimulation imports
from KratosMultiphysics.CoSimulationApplication.co_simulation_tools import classprint

def Create(settings, solver_wrapper):
    cs_tools.SettingsTypeCheck(settings)
    raise NotImplementedError("This class needs some updates see MVQN and Aitken")
    return ConstantRelaxationConvergenceAccelerator(settings, solver_wrapper)

class ConstantRelaxationConvergenceAccelerator(CoSimulationConvergenceAccelerator):
    ## The constructor.
    # @param alpha relaxation factor.
    def __init__( self, settings, solvers, cosim_solver_details ):
        super(ConstantRelaxationConvergenceAccelerator, self).__init__(settings, solvers, cosim_solver_details)
        if "alpha" in self.settings:
            self.alpha = self.settings["alpha"]
        else:
            self.alpha = 0.125

    ## UpdateSolution(r, x)
    # @param r residual r_k
    # @param x solution x_k
    # Computes the approximated update in each iteration.
    def UpdateSolution( self, r, x ):
        if self.echo_level > 3:
            classprint(self._ClassName(), "Doing relaxation with factor = ", "{0:.1g}".format(self.alpha))
        delta_x = self.alpha * r
        return delta_x
