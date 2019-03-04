from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Other imports
import numpy as np
import co_simulation_tools as cs_tools

class CosimulationBasePredictor(object):
    def __init__(self, settings, solvers, level):
        self.settings = settings
        self.solvers = solvers
        self.lvl = level
        self.echo_level = 0

    def Initialize(self):
        pass

    def Finalize(self):
        pass

    def InitializeSolutionStep(self):
        pass

    def Predict(self):
        raise Exception('"Predict" has to be implemented in the derived class!')

    def FinalizeSolutionStep(self):
        pass

    def SetDeltaTime(self, delta_time):
        self.delta_time = delta_time

    def PrintInfo(self):
        '''Function to print Info abt the Object
        Can be overridden in derived classes to print more information
        '''
        cs_tools.classprint(self.lvl, "Predictor", cs_tools.bold(self._Name()))

    def Check(self):
        print("The predictors do not yet implement Check!")

    def SetEchoLevel(self, level):
        self.echo_level = level

    def _Name(self):
        raise Exception('"_Name" has to be implemented in the derived class!')

    def _UpdateData(self, updated_data):
        for data_entry, data_update in zip(self.settings["data_list"], updated_data):
            solver = self.solvers[data_entry["solver"]]
            data_name = data_entry["data_name"]
            cs_tools.ExportArrayToSolver(solver, data_name, data_update)

        if self.echo_level > 3:
            cs_tools.classprint(self.lvl, self._Name(), "Computed prediction")