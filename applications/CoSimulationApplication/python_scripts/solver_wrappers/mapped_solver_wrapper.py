from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics as KM

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_solver_wrapper import CoSimulationSolverWrapper
import KratosMultiphysics.CoSimulationApplication.factories.solver_wrapper_factory as solver_wrapper_factory

# Other imports
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools

def Create(settings, solver_name):
    return MappedSolverWrapper(settings, solver_name)

class MappedSolverWrapper(CoSimulationSolverWrapper):
    """This class wraps another solver-wrapper and  in addition takes care of the mapping
    """
    def __init__(self, settings, solver_name):
        super(MappedSolverWrapper, self).__init__(settings, solver_name)

        self.wrapped_solver = solvers_wrapper_factory.CreateSolverWrapper(self.settings["solver_wrapper_settings"]["wrapped_solver_settings"], solver_name)
        self.input_mapper = None
        self.output_mapper = None


    def Initialize(self):
        super(MappedSolverWrapper, self).Initialize()

    def Finalize(self):
        pass

    def AdvanceInTime(self, current_time):
        pass

    def InitializeSolutionStep(self):
        pass

    def Predict(self):
        pass

    def SolveSolutionStep(self):
        pass

    def FinalizeSolutionStep(self):
        pass

    def OutputSolutionStep(self):
        pass