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

    # TODO: Do we need this?
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
        #self._rigid_body_solver.SetSolutionStepValue("ROOT_POINT_DISPLACEMENT", self._rigid_body_solver.root_point_model_part[KM.DISPLACEMENT], 0)
        #self._rigid_body_solver.SetSolutionStepValue("ROOT_POINT_ROTATION",     self._rigid_body_solver.root_point_model_part[KM.ROTATION], 0)
        #self._rigid_body_solver.SetSolutionStepValue("FORCE",                   self._rigid_body_solver.rigid_body_model_part[KM.FORCE], 0)
        #self._rigid_body_solver.SetSolutionStepValue("MOMENT",                  self._rigid_body_solver.rigid_body_model_part[KM.MOMENT], 0)

        self._rigid_body_solver.SolveSolutionStep()

        #self._rigid_body_solver.rigid_body_model_part[KM.DISPLACEMENT] =    self._rigid_body_solver.GetSolutionStepValue("DISPLACEMENT", 0)
        #self._rigid_body_solver.rigid_body_model_part[KM.ROTATION] =        self._rigid_body_solver.GetSolutionStepValue("ROTATION", 0)
        #self._rigid_body_solver.root_point_model_part[KM.REACTION] =        self._rigid_body_solver.GetSolutionStepValue("REACTION", 0)
        #self._rigid_body_solver.root_point_model_part[KM.REACTION_MOMENT] = self._rigid_body_solver.GetSolutionStepValue("REACTION_MOMENT", 0)

    def Check(self):
        '''
        # making sure only a set of vaiables can be used
        admissible_variables = [
            "ROOT_POINT_DISPLACEMENT",
            "ROOT_POINT_ROTATION",
            "FORCE",
            "MOMENT",
            "DISPLACEMENT",
            "ROTATION",
            "REACTION",
            "REACTION_MOMENT"
        ]
        for data in self.data_dict.values():
            if data.variable.Name() not in admissible_variables:
                raise Exception('Variable "{}" of interface data "{}" of solver "{}" cannot be used for the Rigid Body Solver!\nOnly the following variables are allowed: {}'.format(data.variable.Name(), data.name, data.solver_name, admissible_variables))
        '''
        for process in self._rigid_body_solver._list_of_processes:
            process.Check()
    
    def _GetDataCommunicator(self):
        # this solver does not support MPI
        return GetRankZeroDataCommunicator()
