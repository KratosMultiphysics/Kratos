from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Importing the base class
from co_simulation_base_predictor import CosimulationBasePredictor

# Other imports
import numpy as np
import co_simulation_tools as cs_tools

def Create(predictor_settings, solvers, level):
    return LinearDerivativeBasedStaticPredictor(predictor_settings, solvers, level)

class LinearDerivativeBasedStaticPredictor(CosimulationBasePredictor):

    def Predict(self):
        pass

    def _Name(self):
        return self.__class__.__name__