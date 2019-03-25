from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Other imports
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as tools


##
#  IMPORTANT : This is a BASE CLASS
#               Please do not change any thing in this class.
#
# This Class servers as a base class for all the filters which are to be applied on the co simulation data
# used in co simulation.
class CoSimulationBaseFilter(object):
    def __init__(self, settings, data):
        self.data = data
        self.settings = settings
        self.echo_level = 0

    ## Initialize : Initialize function for the accelerator class. Necessary
    #               initialization of the variables and objects to be done here.
    #  @param self                      The object pointer.
    def Initialize(self):
        pass

    ## Finalize :  Initialize function for the accelerator class.
    #               finalization of the variables and objects to be done here.
    #  @param self                      The object pointer.
    def Finalize(self):
        pass

    ## Apply : Apply function for the accelerator class. Necessary
    #               Calculation of the update and applying update to the data field will happen here.
    #  @param self                      The object pointer.
    def Apply(self):
        raise NotImplementedError(tools.bcolors.FAIL + "CoSimulationBaseFilter : Calling Apply function from base filter !" + tools.bcolors.ENDC)

    ## InitializeSolutionStep : Called once in the beginning of the solution step
    #

    def InitializeSolutionStep(self):
        pass

    ## FinalizeSolutionStep : Called once at the end of the solution step
    #

    def FinalizeSolutionStep(self):
        pass

    ## InitializeCouplingIteration : Called once in the beginning of the coupled iteration
    #
    #  @param self            The object pointer.
    def InitializeCouplingIteration(self):
        raise NotImplementedError(tools.bcolors.FAIL + "CoSimulationBaseFilter : Calling InitializeCouplingIteration function !" + tools.bcolors.ENDC)
    ## FinalizeCouplingIteration : Called once at the end of the coupled iteration
    #
    #  @param self            The object pointer.
    def FinalizeCouplingIteration(self):
        raise NotImplementedError(tools.bcolors.FAIL + "CoSimulationBaseFilter : Calling FinalizeCouplingIteration function !" + tools.bcolors.ENDC)

    ## ComputeUpdate : Function to compute the update for the fields. Should be called during the nonlinear (coupling iteration).
    #
    def ComputeUpdate(self):
        raise NotImplementedError(tools.bcolors.FAIL + "CoSimulationBaseFilter : Calling ComputeUpdate function !" + tools.bcolors.ENDC)

    ## PrintInfo : Function to print the information of the convergence accelerator
    #

    def PrintInfo(self):
        raise NotImplementedError(tools.bcolors.FAIL + "CoSimulationBaseFilter : Calling PrintInfo function !" + tools.bcolors.ENDC)

    ## SetEchoLevel : Function to set the echo_level of the convergence accelerator
    #

    #  @param level           int : echo level to be set
    def SetEchoLevel(self, level):
        self.echo_level = level

    ## Check : Function to Check the setup of the convergence accelerator
    #

    def Check(self):
        raise NotImplementedError(tools.bcolors.FAIL + "CoSimulationBaseFilter : Calling Check function !" + tools.bcolors.ENDC)

    ## _Name : Function to get the name of the convergence accelerator
    #

    def _Name(self):
        raise NotImplementedError(tools.bcolors.FAIL + "CoSimulationBaseFilter : Calling _Name function !" + tools.bcolors.ENDC)
