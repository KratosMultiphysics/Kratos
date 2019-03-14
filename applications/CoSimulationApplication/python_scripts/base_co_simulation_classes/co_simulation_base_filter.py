# co simulation imports
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
from KratosMultiphysics.CoSimulationApplication.base_co_simulation_classes.co_simulation_base_solver import CoSimulationBaseSolver
# Other imports
import KratosMultiphysics.CoSimulationApplication.co_simulation_data_structure as data_str
cs_data_structure = data_str.__DATA_STRUCTURE__
import collections


##
#  IMPORTANT : This is a BASE CLASS
#               Please do not change any thing in this class.
#
#  This class is intended to serve as the base class for all the filters (predictors, accelerators, etc) applied on the data field.
class CoSimulationBaseDataFilter(object):
    ## Constructor :  The base class for the base filter class
    #  @param cosim_solver_settings     parameter object with the filter settings
    def __init__(self, filter_settings):
        default_setting = cs_data_structure.Parameters("""
        {
            "type" : "",
            "settings":{}
        }
        """)
        self.name = solver_name
        self.filter_settings = filter_settings
        self.filter_settings.ValidateAndAssignDefaults(default_setting)

    ## Initialize : Initialize function. Called only once
    #               all member variables are initialized here.
    def Initialize(self):
        pass

    ## Finalize : Finalize function. Called only once
    #               all member variables are finalized here.
    def Finalize(self):
        pass

    ## InitializeSolutionStep : Called once in the beginning of the solution step
    #
    #  @param self            The object pointer.
    def InitializeSolutionStep(self):
        raise NotImplementedError(tools.bcolors.FAIL + "CoSimulationBaseDataFilter : Calling InitializeSolutionStep function !" + tools.bcolors.ENDC)

    ## FinalizeSolutionStep : Called once at the end of the solution step
    #
    #  @param self            The object pointer.
    def InitializeSolutionStep(self):
        raise NotImplementedError(tools.bcolors.FAIL + "CoSimulationBaseDataFilter : Calling InitializeSolutionStep function !" + tools.bcolors.ENDC)

    ## FinalizeNonLinearIteration : Function initializes the non linear iteration (coupling iteration)
    #                               Called at the beginning of the nonlinear iteration (coupling iteration)
    #
    def InitializeNonLinearIteration(self):
        raise NotImplementedError(tools.bcolors.FAIL + "CoSimulationBaseDataFilter : Calling InitializeNonLinearIteration function !" + tools.bcolors.ENDC)

    ## FinalizeNonLinearIteration : Function finalizes the non linear iteration (coupling iteration)
    #                               Called at the end of the nonlinear iteration (coupling iteration)
    #
    def FinalizeNonLinearIteration(self):
        raise NotImplementedError(tools.bcolors.FAIL + "CoSimulationBaseDataFilter : Calling FinalizeNonLinearIteration function !" + tools.bcolors.ENDC)

    ## Apply : Apply function for the accelerator class. Necessary
    #               Calculation of the update and applying update to the data field will happen here.
    #  @param self                      The object pointer.
    def Apply(self):
        raise NotImplementedError(tools.bcolors.FAIL + "CoSimulationBaseDataFilter : Calling Apply function from base accelerator !" + tools.bcolors.ENDC)