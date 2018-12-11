from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Importing the base class
from co_simulation_base_predictor import CosimulationBasePredictor

# Other imports
import numpy as np
import co_simulation_tools as cs_tools

# Predictor implemented according to:
# "A new staggered scheme for fluid-structure interaction"; W.G. Dettmer and D. Peric
# Numerical Methods in Engineering 2013; 93; 1-22

def Create(predictor_settings, solvers, level):
    return AverageValuePredictor(predictor_settings, solvers, level)

class AverageValuePredictor(CosimulationBasePredictor):
    # @param beta factor for weighting last and current value of the predicted values. Can be set in interval: [0, 1.0]
    def __init__(self, settings, solvers, level):
        super(AverageValuePredictor, self).__init__(settings, solvers, level)
        if "beta" in self.settings:
            self.beta = self.settings["beta"]
            if self.beta > 1.0 or self.beta < 0:
                raise Exception("Wrong value for beta. Admissible interval [0.0, 1.0]")
        else:
            self.beta = 0.5
        # TODO check buffer_size

        # TODO add comment why we do this
        num_data = len(self.settings["data_list"])
        self.data_arrays_prediction = [np.array([]) for e in range(num_data)]
        self.data_arrays_aux = [np.array([]) for e in range(num_data)]

    def Predict(self):
        for i, data_entry in enumerate(self.settings["data_list"]):
            solver = self.solvers[data_entry["solver"]]
            data_name = data_entry["data_name"]
            cs_tools.ImportArrayFromSolver(solver, data_name, self.data_arrays_prediction[i], 0) 
            cs_tools.ImportArrayFromSolver(solver, data_name, self.data_arrays_aux[i], 1) 

            self.data_arrays_prediction[i] = 2*self.data_arrays_prediction[i] - self.data_arrays_aux[i]

        self._UpdateData(self.data_arrays_prediction)

    def FinalizeSolutionStep(self):
        for i, data_entry in enumerate(self.settings["data_list"]):
            solver = self.solvers[data_entry["solver"]]
            data_name = data_entry["data_name"]
            cs_tools.ImportArrayFromSolver(solver, data_name, self.data_arrays_aux[i], 0)

            self.data_arrays_prediction[i] = self.beta * self.data_arrays_aux[i] + (1-self.beta) * self.data_arrays_prediction[i]

        self._UpdateData(self.data_arrays_prediction)


    def _Name(self):
        return self.__class__.__name__

