## @module aitken
# This module contains the class Aitken
# Author: Wei He
# Updated : Aditya Ghantasala
# Date: Feb. 20, 2017

from base_classes.co_simulation_base_predictor import CoSimulationBasePredictor
import co_simulation_tools as cs_tools
data_structure = cs_tools.cs_data_structure

try :
    import numpy as np
except ModuleNotFoundError:
    cs_tools.PrintWarning(cs_tools.bcolors.FAIL + 'Numpy is not available', 'Using python default lists for computation'+ cs_tools.bcolors.ENDC)

def Create(settings, solver):
    predictor = LinearPredictor(settings, solver)
    return predictor


## Class LinearPredictor.
# This class contains the implementation for a linear predictor

class LinearPredictor(CoSimulationBasePredictor):

    def __init__( self, settings, solver ):
        super(LinearPredictor, self).__init__(settings, solver)
        default_settings = self._GetDefaultSettings()
        self.settings.RecursivelyValidateAndAssignDefaults(default_settings)
        self.data_prev_iter = []
        self.data_current_iter = []
        self.ResidualStorage = deque( maxlen = 2 ) # 0 is the latest , 1 is the previous


    def _GetDefaultSettings(self):
        default_setting = data_structure.Parameters("""
            {
                "type"          : "linear",
                "data" : {
                        "solver"   : "",
                        "data_name"     : ""
                },
                "settings":{
                }
            }
        """)
        return default_setting


    def InitializeSolutionStep(self):
        self.iteration = 0
        self.alpha_old = self.initial_alpha
        self.data_prev_iter = self.data_current_iter

    def FinalizeSolutionStep(self):
        pass

    ## InitializeNonLinearIteration : Function initializes the non linear iteration (coupling iteration)
    #                               Called at the beginning of the nonlinear iteration (coupling iteration)
    #
    def InitializeNonLinearIteration(self):
        pass

    ## FinalizeNonLinearIteration : Function finalizes the non linear iteration (coupling iteration)
    #                               Called at the end of the nonlinear iteration (coupling iteration)
    #
    def FinalizeNonLinearIteration(self):
        pass

    ## PrintInfo : Function to print the information of the convergence accelerator
    #
    def PrintInfo(self):
        cs_tools.PrintInfo(cs_tools.bcolors.HEADER + "This is an object of linear predictor." +cs_tools.bcolors.ENDC )

    ## Check : Function to Check the setup of the convergence accelerator
    #
    def Check(self):
        pass

    ## _Name : Function to get the name of the convergence accelerator
    #
    def _Name(self):
        return "linear_predictor"