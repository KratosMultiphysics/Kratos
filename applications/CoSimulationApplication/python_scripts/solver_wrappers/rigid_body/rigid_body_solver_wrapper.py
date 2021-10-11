# Importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.CoSimulationApplication as KMC

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_solver_wrapper import CoSimulationSolverWrapper

# Other imports
from .rigid_body_solver import RigidBodySolver

def Create(settings, model, solver_name):
    return RigidBodySolverWrapper(settings, model, solver_name)

class RigidBodySolverWrapper(CoSimulationSolverWrapper):
    """ This class implements a wrapper for an Rigid Body solver to be used in CoSimulation
    """
    def __init__(self, settings, model, solver_name, model_part_name="RigidBody"):
        super().__init__(settings, model, solver_name)

        self.mp = self.model.CreateModelPart(model_part_name)
        self.mp.ProcessInfo[KM.DOMAIN_SIZE] = 3

        input_file_name = self.settings["solver_wrapper_settings"]["input_file"].GetString()
        self._rigid_body_solver = self._CreateRigidBodySolver(input_file_name, self.mp)

    @classmethod
    def _CreateRigidBodySolver(cls, input_file_name, model_part):
        return RigidBodySolver(input_file_name, model_part)

    def Initialize(self):
        super().Initialize()
        self._rigid_body_solver.Initialize()

    def OutputSolutionStep(self):
        self._rigid_body_solver.OutputSolutionStep()

    def AdvanceInTime(self, current_time):
        return self._rigid_body_solver.AdvanceInTime(current_time)

    def SolveSolutionStep(self):
        self._rigid_body_solver.SetSolutionStepValue("ROOT_POINT_DISPLACEMENT", self.mp[KMC.ROOT_POINT_DISPLACEMENT], 0)
        self._rigid_body_solver.SetSolutionStepValue("ROOT_POINT_ROTATION",     self.mp[KMC.ROOT_POINT_ROTATION], 0)
        self._rigid_body_solver.SetSolutionStepValue("FORCE",                   self.mp[KM.FORCE], 0)
        self._rigid_body_solver.SetSolutionStepValue("MOMENT",                  self.mp[KM.MOMENT], 0)

        self._rigid_body_solver.SolveSolutionStep()

        self.mp[KM.DISPLACEMENT] =    self._rigid_body_solver.GetSolutionStepValue("DISPLACEMENT", 0)
        self.mp[KM.ROTATION] =        self._rigid_body_solver.GetSolutionStepValue("ROTATION", 0)
        self.mp[KM.REACTION] =        self._rigid_body_solver.GetSolutionStepValue("REACTION", 0)
        self.mp[KM.REACTION_MOMENT] = self._rigid_body_solver.GetSolutionStepValue("REACTION_MOMENT", 0)

    def Check(self):
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
