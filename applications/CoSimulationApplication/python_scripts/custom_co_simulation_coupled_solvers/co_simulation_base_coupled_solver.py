# co simulation imports
import tools
# Importing the CoSimulation application 
import KratosMultiphysics.CoSimulationApplication as CoSimulationApplication
from co_simulation_solver_interfaces.co_simulation_base_solver import CoSimulationBaseSolver

# Other imports
import os

def CreateSolver(custom_settings):
    return GaussSeidelIterativeStrongCouplingSolver(custom_settings)

class CoSimulationBaseCoupledSolver(CoSimulationBaseSolver):
    '''
        This class is intended to server as the base class for all the
        coupled solvers.
    '''
    ## The constructor
    #
    #  @param self            The object pointer.
    #  @param parameters     parameters for configuring the CoSimulationBaseCoupledSolver
    def __init__(self, custom_settings):

        ##settings string in json format
        # default values for all available settings
        # for mandatory settings, the type is defined
        defaultSettings = {}
        defaultSettings["echo_level"] = 1
        defaultSettings["convergence_acceleration"] = dict    #MANDATORY
        defaultSettings["convergence_criteria"] = dict    #MANDATORY
        defaultSettings["max_iteration_per_step"] = 10
        defaultSettings["participants"] = list