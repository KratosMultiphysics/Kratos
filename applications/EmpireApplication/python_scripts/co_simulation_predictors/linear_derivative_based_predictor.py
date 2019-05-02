from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Importing the base class
from co_simulation_base_predictor import CosimulationBasePredictor

# Other imports
import numpy as np
import co_simulation_tools as cs_tools

def Create(predictor_settings, solvers, level):
    return LinearDerivativeBasedPredictor(predictor_settings, solvers, level)

class LinearDerivativeBasedPredictor(CosimulationBasePredictor):
    def __init__(self, settings, solvers, level):
        super(LinearDerivativeBasedPredictor, self).__init__(settings, solvers, level)
        # TODO add comment why we do this
        num_data = len(self.settings["data_list"])
        self.data_arrays            = [np.array([]) for e in range(num_data)]
        self.derivative_data_arrays = [np.array([]) for e in range(num_data)]

        # TODO check buffer size!

    def Predict(self):
        for i, data_entry in enumerate(self.settings["data_list"]):
            solver = self.solvers[data_entry["solver"]]
            data_name = data_entry["data_name"]
            deriv_data_name = data_entry["derivative_data_name"]
            cs_tools.ImportArrayFromSolver(solver, deriv_data_name, self.derivative_data_arrays[i], 1)
            cs_tools.ImportArrayFromSolver(solver, data_name, self.data_arrays[i], 1)

            self.data_arrays[i] += solver.GetDeltaTime() * self.derivative_data_arrays[i]

        self._UpdateData(self.data_arrays)

    def _Name(self):
        return self.__class__.__name__

