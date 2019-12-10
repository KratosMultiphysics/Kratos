from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.CoSimulationApplication as KMC

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_solver_wrapper import CoSimulationSolverWrapper

# Other imports
from .sdof_static_solver import SDoFStaticSolver
from .sdof_solver_wrapper import SdofSolverWrapper

def Create(settings, solver_name):
    return SdofStaticSolverWrapper(settings, solver_name)

class SdofStaticSolverWrapper(SdofSolverWrapper):
    """ This class implements a wrapper for an SDof solver to be used in CoSimulation
    """
    def __init__(self, settings, solver_name):
        CoSimulationSolverWrapper.__init__(self,settings, solver_name)

        input_file_name = self.settings["solver_wrapper_settings"]["input_file"].GetString()

        self.mp = self.model.CreateModelPart("Sdof_Static")
        self.mp.ProcessInfo[KM.DOMAIN_SIZE] = 1
        self._sdof_solver = SDoFStaticSolver(input_file_name)
