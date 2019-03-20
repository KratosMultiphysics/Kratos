from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Other imports
import numpy as np
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools

class CoSimulationBaseConvergenceAccelerator(object):
    def __init__(self, settings, solvers):
        self.settings = settings
        self.solvers = solvers
        self.echo_level = 0

        ## Here we preallocate the arrays that will be used to exchange data
        self.data_array = np.array([])
        self.prev_data = np.array([])

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
        solver = self.solvers[self.settings["solver"]]
        data_name = self.settings["data_name"]
        cs_tools.ImportArrayFromSolver(solver, data_name, self.prev_data)

    def FinalizeNonLinearIteration(self):
        pass

    def ComputeUpdate(self):
        solver = self.solvers[self.settings["solver"]]
        data_name = self.settings["data_name"]

        cs_tools.ImportArrayFromSolver(solver, data_name, self.data_array)

        residual = self.data_array - self.prev_data

        updated_data = self.prev_data + self._ComputeUpdate(residual, self.prev_data)

        cs_tools.ExportArrayToSolver(solver, data_name, updated_data)

    def PrintInfo(self):
        '''Function to print Info abt the Object
        Can be overridden in derived classes to print more information
        '''
        cs_tools.classprint("Convergence Accelerator", cs_tools.bold(self._Name()))

    def SetEchoLevel(self, level):
        self.echo_level = level

    def _Name(self):
        raise Exception('"_Name" has to be implemented in the derived class!')

    def Check(self):
        print("ConvAcc does not yet implement Check")

    def _ComputeUpdate( self, residual, previous_data ):
        raise Exception('"_ComputeUpdate" has to be implemented in the derived class!')
