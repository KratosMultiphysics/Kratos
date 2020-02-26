from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics as KM

# CoSimulation imports
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
import KratosMultiphysics.CoSimulationApplication.colors as colors

class CoSimulationConvergenceAccelerator(object):
    """Baseclass for the convergence acceleratos used for CoSimulation
    Relaxes the solution to increase the speed of convergence in a (strongly) coupled simulation

    Note that the interface matches the convergence accelerators in the FSIApplication such that they can be used interchangeable
    ("FSIApplication/custom_utilities/convergence_accelerator.hpp")
    """
    def __init__(self, settings):
        self.settings = settings
        self.settings.RecursivelyValidateAndAssignDefaults(self._GetDefaultSettings())

        self.echo_level = self.settings["echo_level"].GetInt()

    def Initialize(self):
        pass

    def Finalize(self):
        pass

    def InitializeSolutionStep(self):
        pass

    def FinalizeSolutionStep(self):
        pass

    def InitializeNonLinearIteration(self):
        pass

    def FinalizeNonLinearIteration(self):
        pass

    def PrintInfo(self):
        '''Function to print Info abt the Object
        Can be overridden in derived classes to print more information
        '''
        cs_tools.cs_print_info("Convergence Accelerator", colors.bold(self._ClassName()))

    def Check(self):
        print("ConvAcc does not yet implement Check")

    def UpdateSolution(self, residual, iteration_guess):
        # TODO this should update the solution in place, otherwise not compatible with the conv-acc in the FSI-App
        # => would probably not be compatible with any C++ Conv-Acc
        raise NotImplementedError('"UpdateSolution" has to be implemented in the derived class!')

    @classmethod
    def SupportsDistributedData(cls):
        return False

    @classmethod
    def _ClassName(cls):
        return cls.__name__

    @classmethod
    def _GetDefaultSettings(cls):
        return KM.Parameters("""{
            "type"       : "UNSPECIFIED",
            "echo_level" : 0
        }""")
