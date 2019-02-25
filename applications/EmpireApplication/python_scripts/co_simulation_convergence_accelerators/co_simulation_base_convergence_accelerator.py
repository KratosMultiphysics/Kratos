from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Other imports
import numpy as np
import co_simulation_tools as cs_tools

class CoSimulationBaseConvergenceAccelerator(object):
    def __init__(self, settings, solvers, level):
        self.settings = settings
        self.solvers = solvers
        self.lvl = level
        self.echo_level = 0

        ## Here we preallocate the arrays that will be used to exchange data
        num_data = len(self.settings["data_list"])
        self.arrays = [np.array([]) for e in range(num_data)]

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
        self.data_sizes = [] # saving the sizes of the data to later split them again
        size_counter = 0
        for i, data_entry in enumerate(self.settings["data_list"]):
            solver = self.solvers[data_entry["solver"]]
            data_name = data_entry["data_name"]
            cs_tools.ImportArrayFromSolver(solver, data_name, self.arrays[i])

            size_counter += self.arrays[i].size
            self.data_sizes.append(size_counter)

        self.combined_prev_data = np.concatenate(self.arrays)

    def FinalizeNonLinearIteration(self):
        pass

    def ComputeUpdate(self):
        for i, data_entry in enumerate(self.settings["data_list"]):
            solver = self.solvers[data_entry["solver"]]
            data_name = data_entry["data_name"]
            cs_tools.ImportArrayFromSolver(solver, data_name, self.arrays[i])

        combined_residuals = np.concatenate(self.arrays) - self.combined_prev_data

        combined_new_data = self.combined_prev_data + self._ComputeUpdate(combined_residuals, self.combined_prev_data)

        updated_data = np.split(combined_new_data, self.data_sizes)

        for data_entry, data_update in zip(self.settings["data_list"], updated_data):
            solver = self.solvers[data_entry["solver"]]
            data_name = data_entry["data_name"]
            cs_tools.ExportArrayToSolver(solver, data_name, data_update)

    def PrintInfo(self):
        '''Function to print Info abt the Object
        Can be overridden in derived classes to print more information
        '''
        cs_tools.classprint(self.lvl, "Convergence Accelerator", cs_tools.bold(self._Name()))

    def SetEchoLevel(self, level):
        self.echo_level = level

    def _Name(self):
        raise Exception('"_Name" has to be implemented in the derived class!')

    def Check(self):
        print("ConvAcc does not yet implement Check")

    def _ComputeUpdate( self, residual, previous_data ):
        raise Exception('"_ComputeUpdate" has to be implemented in the derived class!')
