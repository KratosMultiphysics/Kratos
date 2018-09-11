# co simulation imports
import co_simulation_tools as tools
# Importing the CoSimulation application 
import KratosMultiphysics.CoSimulationApplication as CoSimulationApplication
from custom_co_simulation_solver_interfaces.co_simulation_base_solver import CoSimulationBaseSolver

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
        self.full_settings = custom_settings
        defaultSettings = {}
        defaultSettings["echo_level"] = 1
        defaultSettings["max_coupling_iterations"] = 10
        defaultSettings["participants"] = list
        defaultSettings["start_coupling_time"] = float #MANDATORY
        self.settings = tools.ValidateAndAssignInputParameters(defaultSettings, custom_settings["coupled_solver_settings"], False)
        self.number_of_participants = len( self.settings['participants'] )

        # Get the participating solvers a map with their names and objects
        self.participating_solvers = tools.GetSolvers(self.full_settings['solvers'])

        # With this setting the coupling can start later
        self.start_coupling_time = 0.0
        if "start_coupling_time" in self.settings:
            self.start_coupling_time = self.settings["start_coupling_time"]
        if self.start_coupling_time > 0.0:
            self.coupling_started = False
        else:
            self.coupling_started = True

    def Initialize(self):
        for solver_name, solver in self.participating_solvers.items():
            solver.Initialize()

    def Finalize(self):
        for solver_name, solver in self.participating_solvers.items():
            solver.Finalize()


    def InitializeSolutionStep(self):
        for solver_name, solver in self.participating_solvers.items():
            solver.InitializeTimeStep()

    def FinalizeSolutionStep(self):
        for solver_name, solver in self.participating_solvers.items():
            solver.FinalizeTimeStep()