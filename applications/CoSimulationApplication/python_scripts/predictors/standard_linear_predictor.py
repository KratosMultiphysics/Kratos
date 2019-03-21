from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_base_predictor import CosimulationBasePredictor

# Other imports
import numpy as np
from KratosMultiphysics.CoSimulationApplication.co_simulation_tools import classprint

def Create(predictor_settings, solvers):
    return StandardLinearPredictor(predictor_settings, solvers)

class StandardLinearPredictor(CosimulationBasePredictor):
    def __init__(self, settings, solvers):
        super(StandardLinearPredictor, self).__init__(settings, solvers)
        # TODO add comment why we do this
        self.data_array_t0 = np.array([])
        self.data_array_t1 = np.array([])

        # TODO check buffer size!

    def Predict(self):
        solver = self.solvers[self.settings["solver"].GetString()]
        data_name = self.settings["data_name"].GetString()
        cs_tools.ImportArrayFromSolver(solver, data_name, self.data_array_t0, 0)
        cs_tools.ImportArrayFromSolver(solver, data_name, self.data_array_t1, 1)

        self.data_array_t0 = 2*self.data_array_t0 - data_array_t1

        self._UpdateData(self.data_array_t0)

    def _Name(self):
        return self.__class__.__name__

