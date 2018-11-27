## @module aitken
# This module contains the class Aitken
# Author: Wei He
# Updated : Aditya Ghantasala
# Date: Feb. 20, 2017
try :
    import numpy as np
except ModuleNotFoundError:
    print(cs_tools.bcolors.FAIL + 'Numpy is not available ! Using python default lists for computation'+ cs_tools.bcolors.ENDC)

from base_co_simulation_classes.co_simulation_base_extrapolator import CoSimulationBaseExtrapolator
import co_simulation_tools as cs_tools
data_structure = cs_tools.cs_data_structure


def Create(settings, solver):
    extrapolator = LinearExtrapolator(settings, solver)
    return extrapolator


## Class LinearExtrapolator.
# This class contains the implementation for a linear extrapolator

class LinearExtrapolator(CoSimulationBaseExtrapolator):

    def __init__( self, settings, solver ):
        super(LinearExtrapolator, self).__init__(settings, solver)
        default_settings = self._GetDefaultSettings()
        self.settings.RecursivelyValidateAndAssignDefaults(default_settings)
        self.data_prev_iter = []
        self.data_current_iter = []


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
        self.data_current_iter = cs_tools.GetDataAsList(self.solver, self.data_name)
        self._CalculateResidual()
        self._ComputeUpdatedRelaxFactor()
        self._ApplyRelaxationToData()

    ## FinalizeNonLinearIteration : Function finalizes the non linear iteration (coupling iteration)
    #                               Called at the end of the nonlinear iteration (coupling iteration)
    #
    def FinalizeNonLinearIteration(self):
        self.data_prev_iter = self.data_current_iter
        self.iteration = self.iteration + 1

    ## PrintInfo : Function to print the information of the convergence accelerator
    #
    def PrintInfo(self):
        print(cs_tools.bcolors.HEADER + "This is an object of linear extrapolator." +cs_tools.bcolors.ENDC )

    ## Check : Function to Check the setup of the convergence accelerator
    #
    def Check(self):
        pass

    ## _Name : Function to get the name of the convergence accelerator
    #
    def _Name(self):
        return "linear_extrapolator"