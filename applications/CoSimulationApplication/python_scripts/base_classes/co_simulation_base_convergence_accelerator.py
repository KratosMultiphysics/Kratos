from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Other imports
import numpy as np
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools

class CoSimulationBaseConvergenceAccelerator(object):
    def __init__(self, settings, solver):
        self.settings = settings
        self.settings.RecursivelyValidateAndAssignDefaults(self._GetDefaultSettings())

        self.solver = solver
        self.interface_data = self.solver.GetInterfaceData(self.settings["data_name"].GetString())

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
        # Saving the previous data for the computation of the residual
        # and the computation of the solution update
        self.input_data = self.interface_data.GetNumpyArray()

    def FinalizeNonLinearIteration(self):
        pass

    def ComputeUpdate(self):
        current_data = self.interface_data.GetNumpyArray()
        residual = current_data - self.input_data
        updated_data = self.input_data + self._ComputeUpdate(residual, self.input_data)
        self.interface_data.ApplyUpdateToData(updated_data)

    def PrintInfo(self):
        '''Function to print Info abt the Object
        Can be overridden in derived classes to print more information
        '''
        cs_tools.classprint("Convergence Accelerator", cs_tools.bold(self._Name()))

    def SetEchoLevel(self, level):
        self.echo_level = level

    def _Name(self):
        return self.__class__.__name__

    def Check(self):
        print("ConvAcc does not yet implement Check")

    def _ComputeUpdate( self, residual, previous_data ):
        raise Exception('"_ComputeUpdate" has to be implemented in the derived class!')

    @classmethod
    def _GetDefaultSettings(cls):
        return cs_tools.cs_data_structure.Parameters("""{
            "type"       : "UNSPECIFIED",
            "solver"     : "UNSPECIFIED",
            "data_name"  : "UNSPECIFIED",
            "echo_level" : 0
        }""")
