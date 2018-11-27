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


## Class AitkenAccelerator.
# This class contains the implementation of Aitken relaxation.
# Reference: Ulrich Küttler et al., "Fixed-point fluid–structurce interaction solvers with dynamic relaxation"
class AitkenAccelerator(CoSimulationBaseConvergenceAccelerator):
    ## The constructor.
    # @param settings Settings for the relaxation
    # @param solver The solver which the relaxation should work with
    def __init__( self, settings, solver ):
        super(AitkenAccelerator, self).__init__(settings, solver)
        default_settings = self._GetDefaultSettings()
        self.settings.RecursivelyValidateAndAssignDefaults(default_settings)
        self.ResidualStorage = deque( maxlen = 2 ) # 0 is the latest , 1 is the previous
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

    ## _GetDefaultSettings :  Function to define the default parameters for this class
    #
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


    ## _ComputeUpdatedRelaxFactor: Function computes the new relaxation factor for the current non-linear iteration
    #
    #
    def _ComputeUpdatedRelaxFactor( self ):
        ## For the first iteration, do relaxation only
        if self.iteration == 0:
            alpha = self.initial_alpha
            print( cs_tools.bcolors.BLUE + "\tAitken: Doing relaxation in the first iteration with initial factor = " , alpha, cs_tools.bcolors.ENDC)
            self.current_alpha = alpha
            self.update = [data * self.current_alpha for data in self.ResidualStorage[0]]
        else:
            r_diff = self._Difference(self.ResidualStorage[1] , self.ResidualStorage[0])
            numerator = cs_tools.InnerProduct( self.ResidualStorage[1], r_diff )
            denominator = math.sqrt( cs_tools.InnerProduct( r_diff, r_diff ) )
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
            self.alpha_old = self.current_alpha

    ## InitializeSolutionStep : Called once at the beginning of the solution step
    #
    def _CalculateUpdate(self):
        self.update = [data * self.current_alpha for data in self.ResidualStorage[0]]

    ## InitializeSolutionStep : Called once at the beginning of the solution step
    #
    def InitializeSolutionStep(self):
        self.iteration = 0
        self.alpha_old = self.initial_alpha

    ## FinalizeSolutionStep : Called once at the end of the solution step
    #
    def FinalizeSolutionStep(self):
        pass

    ## InitializeNonLinearIteration : Function initializes the non linear iteration (coupling iteration)
    #                               Called at the beginning of the nonlinear iteration (coupling iteration)
    #
    def InitializeNonLinearIteration(self):
        self.input_data_current_iter = cs_tools.GetDataAsList(self.solver, self.data_name)


    ## FinalizeNonLinearIteration : Function finalizes the non linear iteration (coupling iteration)
    #                               Called at the end of the nonlinear iteration (coupling iteration)
    #
    def FinalizeNonLinearIteration(self):
        self.output_data_current_iter = cs_tools.GetDataAsList(self.solver, self.data_name)
        self._CalculateResidual()
        self._ComputeUpdatedRelaxFactor()
        self._CalculateUpdate()
        self._ApplyRelaxationToData()
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
    #                       residual = output_data_current_iter - input_data_current_iter
    def _CalculateResidual(self):
        if(self.iteration == 0):
            self.ResidualStorage.appendleft( self.output_data_current_iter )
            return

        self.ResidualStorage.appendleft( self._Difference(self.output_data_current_iter , self.input_data_current_iter) )

    ## _ApplyRelaxationToData : updates the data with the update calculated
    #
    def _ApplyRelaxationToData(self):
        updated_data = [ input_data + update  for input_data, update in zip(self.input_data_current_iter, self.update)]
        cs_tools.ApplyUpdateToData(self.solver, self.data_name, updated_data)


    ## _Difference : Calculates the difference of two vectors provided in terms of python lists
    #
    def _Difference(self, list_one, list_two):
        diff = [ i-j for i,j in zip(list_one, list_two)]

        return diff

