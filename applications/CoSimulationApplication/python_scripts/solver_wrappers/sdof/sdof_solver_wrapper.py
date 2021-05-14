# Importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.CoSimulationApplication as KMC

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_solver_wrapper import CoSimulationSolverWrapper

# Other imports
from .sdof_solver import SDoFSolver

def Create(settings, model, solver_name):
    return SdofSolverWrapper(settings, model, solver_name)

class SdofSolverWrapper(CoSimulationSolverWrapper):
    """ This class implements a wrapper for an SDof solver to be used in CoSimulation
    """
    def __init__(self, settings, model, solver_name, model_part_name="Sdof"):
        super().__init__(settings, model, solver_name)

        self.mp = self.model.CreateModelPart(model_part_name)
        self.mp.ProcessInfo[KM.DOMAIN_SIZE] = 1

        input_file_name = self.settings["solver_wrapper_settings"]["input_file"].GetString()
        self._sdof_solver = self._CreateSDofSolver(input_file_name)

    @classmethod
    def _CreateSDofSolver(cls, input_file_name):
        return SDoFSolver(input_file_name)

    def Initialize(self):
        super().Initialize()
        self._sdof_solver.Initialize()

    def OutputSolutionStep(self):
        self._sdof_solver.OutputSolutionStep()

    def AdvanceInTime(self, current_time):
        return self._sdof_solver.AdvanceInTime(current_time)

    def SolveSolutionStep(self):
        self._sdof_solver.SetSolutionStepValue("ROOT_POINT_DISPLACEMENT", self.mp[KMC.SCALAR_ROOT_POINT_DISPLACEMENT], 0)
        self._sdof_solver.SetSolutionStepValue("LOAD",                    self.mp[KMC.SCALAR_FORCE], 0)

        self._sdof_solver.SolveSolutionStep()

        self.mp[KMC.SCALAR_DISPLACEMENT]        = self._sdof_solver.GetSolutionStepValue("DISPLACEMENT", 0)
        self.mp[KMC.SCALAR_REACTION]            = self._sdof_solver.GetSolutionStepValue("REACTION", 0)
        self.mp[KMC.SCALAR_VOLUME_ACCELERATION] = self._sdof_solver.GetSolutionStepValue("VOLUME_ACCELERATION", 0)

    def Check(self):
        # making sure only a set of vaiables can be used
        admissible_variables = [
            "SCALAR_ROOT_POINT_DISPLACEMENT",
            "SCALAR_FORCE",
            "SCALAR_DISPLACEMENT",
            "SCALAR_REACTION",
            "SCALAR_VOLUME_ACCELERATION",
        ]
        for data in self.data_dict.values():
            if data.variable.Name() not in admissible_variables:
                raise Exception('Variable "{}" of interface data "{}" of solver "{}" cannot be used for the SDof Solver!\nOnly the following variables are allowed: {}'.format(data.variable.Name(), data.name, data.solver_name, admissible_variables))
