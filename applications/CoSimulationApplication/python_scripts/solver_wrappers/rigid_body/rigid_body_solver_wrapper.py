# Importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.CoSimulationApplication as KMC

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_solver_wrapper import CoSimulationSolverWrapper

# Other imports
from .rigid_body_solver import RigidBodySolver
from KratosMultiphysics.CoSimulationApplication.utilities.data_communicator_utilities import GetRankZeroDataCommunicator

def Create(settings, model, solver_name):
    return RigidBodySolverWrapper(settings, model, solver_name)

class RigidBodySolverWrapper(CoSimulationSolverWrapper):
    """ This class implements a wrapper for an Rigid Body solver to be used in CoSimulation
    """
    def __init__(self, settings, model, solver_name):
        super().__init__(settings, model, solver_name)

        input_file_name = settings["solver_wrapper_settings"]["input_file"].GetString()
        if not input_file_name.endswith(".json"):
            input_file_name += ".json"
        with open(input_file_name,'r') as parameter_file:
            self.project_parameters = KM.Parameters(parameter_file.read())

        input_file_name = self.settings["solver_wrapper_settings"]["input_file"].GetString()
        self._rigid_body_solver = self._CreateRigidBodySolver(self.model, self.project_parameters)

    @classmethod
    def _CreateRigidBodySolver(cls, model, parameters):
        return RigidBodySolver(model, parameters)

    def Initialize(self):
        super().Initialize()
        self._rigid_body_solver.Initialize()
    
    def Finalize(self):
        super().Finalize()
        self._rigid_body_solver.Finalize()

    def OutputSolutionStep(self):
        self._rigid_body_solver.OutputSolutionStep()

    def AdvanceInTime(self, current_time):
        return self._rigid_body_solver.AdvanceInTime(current_time)

    def Predict(self):
        self._rigid_body_solver.Predict()

    def InitializeSolutionStep(self):
        self._rigid_body_solver.InitializeSolutionStep()

    def FinalizeSolutionStep(self):
        self._rigid_body_solver.FinalizeSolutionStep()

    def SolveSolutionStep(self):
        self._rigid_body_solver.SolveSolutionStep()

    def Check(self):
        for process in self._rigid_body_solver._list_of_processes:
            process.Check()
    
    def _GetDataCommunicator(self):
        # this solver does not support MPI
        return GetRankZeroDataCommunicator()
