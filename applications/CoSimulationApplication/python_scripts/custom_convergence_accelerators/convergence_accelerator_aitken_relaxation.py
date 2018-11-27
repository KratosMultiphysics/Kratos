## @module aitken
# This module contains the class Aitken
# Author: Wei He
# Updated : Aditya Ghantasala
# Date: Feb. 20, 2017
try :
    import numpy as np
except ModuleNotFoundError:
    print(cs_tools.bcolors.FAIL + 'Numpy is not available ! Using python default lists for computation'+ cs_tools.bcolors.ENDC)

from base_co_simulation_classes.co_simulation_base_convergence_accelerator import CoSimulationBaseConvergenceAccelerator
import co_simulation_tools as cs_tools
data_structure = cs_tools.cs_data_structure

from copy import deepcopy
from collections import deque
import math

def Create(settings, solver):
    accelerator = AitkenAccelerator(settings, solver)
    return accelerator


## Class Aitken.
# This class contains the implementation of Aitken relaxation and helper functions.
# Reference: Ulrich Küttler et al., "Fixed-point fluid–structurce interaction solvers with dynamic relaxation"
class AitkenAccelerator(CoSimulationBaseConvergenceAccelerator):
    ## The constructor.
    # @param init_alpha Initial relaxation factor in the first time step.
    # @param init_alpha_max Maximum relaxation factor for the first iteration in each time step
    def __init__( self, settings, solver ):
        super(AitkenAccelerator, self).__init__(settings, solver)
        default_settings = self._GetDefaultSettings()
        self.settings.RecursivelyValidateAndAssignDefaults(default_settings)
        self.R = deque( maxlen = 2 )
        self.alpha_max = self.settings["settings"]["alpha_max"].GetDouble()
        self.initial_alpha = self.settings["settings"]["initial_alpha"].GetDouble()
        self.alpha_old = self.initial_alpha
        self.current_alpha = self.initial_alpha
        self.data_name = self.settings["data"]["data_name"].GetString()

        self.iteration = 0

        self.data_prev_iter = []
        self.data_current_iter = []

        self.residual = []
        self.update = []

    def _GetDefaultSettings(self):
        default_setting = data_structure.Parameters("""
            {
                "type"          : "aitken",
                "data" : {
                        "solver"   : "fluid",
                        "data_name"     : ""
                },
                "settings":{
                    "initial_alpha" : 0.2,
                    "alpha_max" : 2.0
                }
            }
        """)
        return default_setting


    ## ComputeUpdate(r, x)
    # @param r residual r_k
    # @param x solution x_k
    # Computes the approximated update in each iteration.
    def _ComputeUpdatedRelaxFactor( self ):
        self.R.appendleft( deepcopy(self.residual) )
        ## For the first iteration, do relaxation only
        if self.iteration == 0:
            alpha = self.initial_alpha
            print( cs_tools.bcolors.BLUE + "\tAitken: Doing relaxation in the first iteration with initial factor = " , alpha, cs_tools.bcolors.ENDC)
            self.current_alpha = alpha
            self.update = [data * self.current_alpha for data in self.R[0]]
        else:
            r_diff = self._Difference(self.R[0] , self.R[1])
            numerator = cs_tools.InnerProduct( self.residual, r_diff )
            denominator = cs_tools.InnerProduct( r_diff, r_diff )
            print("#############################")
            print("Numerator :: ", numerator)
            print("Denominator :: ", denominator)
            print("#############################")
            if(abs(denominator)<1E-15):
                denominator = 1.0
            self.current_alpha = -self.alpha_old * numerator/denominator
            print( cs_tools.bcolors.BLUE + "\tAitken: Doing relaxation with factor = " + cs_tools.bcolors.ENDC, self.current_alpha )
            if self.current_alpha > self.alpha_max:
                self.current_alpha = self.alpha_max
                print(cs_tools.bcolors.WARNING + "WARNING: dynamic relaxation factor reaches upper bound: "+ self.alpha_max + cs_tools.bcolors.ENDC)
            elif self.current_alpha < -self.alpha_max:
                self.current_alpha = -self.alpha_max
                print(cs_tools.bcolors.WARNING + "WARNING: dynamic relaxation factor reaches lower bound: "+ self.alpha_max + cs_tools.bcolors.ENDC)
            self.update = [data * self.current_alpha for data in self.R[0]]
            self.alpha_old = self.current_alpha

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
        print(cs_tools.bcolors.HEADER + "This is an object of Aitken relaxation accelerator. Initial alpha is ", self.init_alpha_max, ", current alpha is : ", self.alpha_old,""+cs_tools.bcolors.ENDC )

    ## Check : Function to Check the setup of the convergence accelerator
    #

    def Check(self):
        pass

    ## _Name : Function to get the name of the convergence accelerator
    #

    def _Name(self):
        return "aitken"

    ## _CalculateResidual : Calculates residual of the data specified in the settings
    #                       Numpy can be used in the variants of this class.
    #                       residual = data_in_current_iter - data_in_previous_iter
    #
    #  @param self            The object pointer
    def _CalculateResidual(self):
        self.residual = []
        if(self.iteration == 0):
            self.residual = self.data_current_iter
            return

        self.residual = self._Difference(self.data_current_iter , self.data_prev_iter)

    def _ApplyRelaxationToData(self):
        cs_tools.ApplyUpdateToData(self.solver, self.data_name, self.update)

    def _Difference(self, list_one, list_two):
        diff = []
        if(len(list_one) == len(list_two)):
            for i in range(len(list_one)):
                diff.append(list_one[i] - list_two[i])

        return diff

