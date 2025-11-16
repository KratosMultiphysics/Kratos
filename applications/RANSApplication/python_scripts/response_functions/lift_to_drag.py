"""
This module contains an interface to the available response functions
"""
import time as timer
from pathlib import Path

import KratosMultiphysics as Kratos
from KratosMultiphysics.response_functions.response_function_interface import ResponseFunctionInterface

from KratosMultiphysics.RANSApplication.response_functions.utilities import SolvePrimalProblem
from KratosMultiphysics.RANSApplication.response_functions.utilities import SolveAdjointProblem
from KratosMultiphysics.RANSApplication.response_functions.utilities import CalculateTimeAveragedDrag
from KratosMultiphysics.RANSApplication.response_functions.utilities import RecursiveCopy

class LiftToDrag(ResponseFunctionInterface):
    def __init__(self, identifier, response_settings, model):
        self.identifier = identifier
        self.response_settings = response_settings

        default_parameters = Kratos.Parameters("""
        {
            "response_type"       : "lift_to_drag_response",
            "problem_setup_folder": "PLEASE_SPECIFY_PROBLEM_SETUP_FOLDER",
            "problem_setup_files" : {
                "primal_project_parameters_file"      : "PLEASE_SPECIFY_PRIMAL_PROJECT_PARAMETERS_FILE",
                "lift_adjoint_project_parameters_file": "PLEASE_SPECIFY_LIFT_ADJOINT_PROJECT_PARAMETERS_FILE",
                "drag_adjoint_project_parameters_file": "PLEASE_SPECIFY_DRAG_ADJOINT_PROJECT_PARAMETERS_FILE"
            }
        }
        """)

        self.response_settings.ValidateAndAssignDefaults(default_parameters)

        self.problem_setup_folder = Path(self.response_settings["problem_setup_folder"].GetString()).absolute()
        if (not self.problem_setup_folder.is_dir()):
            raise Exception("The provided \"problem_setup_folder\" is not a directory [ \"problem_setup_folder\" = {:s} ].".format(str(self.problem_setup_folder)))

        self.problem_setup_file_settings = self.response_settings["problem_setup_files"]

        # checks for lift response settings
        lift_adjoint_project_parameters_file = self.problem_setup_folder / self.problem_setup_file_settings["lift_adjoint_project_parameters_file"].GetString()
        with open(lift_adjoint_project_parameters_file, "r") as file_input:
            lift_adjoint_settings = Kratos.Parameters(file_input.read())

        lift_adjoint_response_function_settings = lift_adjoint_settings["solver_settings"]["response_function_settings"]
        lift_response_type = lift_adjoint_response_function_settings["response_type"].GetString()
        if (lift_response_type != "drag"):
            raise Exception("Lift to drag needs \"response_type\" in \"response_function_settings\" of \"lift_adjoint_project_parameters_file\" file to be drag. [ \"lift_adjoint_project_parameters_file\" = {:s}, \"response_type\" = {:s} ].".format(
                lift_adjoint_project_parameters_file, lift_response_type))

        # check for drag response settings
        drag_adjoint_project_parameters_file = self.problem_setup_folder / self.problem_setup_file_settings["drag_adjoint_project_parameters_file"].GetString()
        with open(drag_adjoint_project_parameters_file, "r") as file_input:
            drag_adjoint_settings = Kratos.Parameters(file_input.read())

        drag_adjoint_response_function_settings = drag_adjoint_settings["solver_settings"]["response_function_settings"]
        drag_response_type = drag_adjoint_response_function_settings["response_type"].GetString()
        if (drag_response_type != "drag"):
            raise Exception("Lift to drag needs \"response_type\" in \"response_function_settings\" of \"drag_adjoint_project_parameters_file\" file to be drag. [ \"drag_adjoint_project_parameters_file\" = {:s}, \"response_type\" = {:s} ].".format(
                drag_adjoint_project_parameters_file, drag_response_type))

        self.lift_model_part_name = lift_adjoint_response_function_settings["custom_settings"]["structure_model_part_name"].GetString()
        self.lift_direction = lift_adjoint_response_function_settings["custom_settings"]["drag_direction"].GetVector()

        self.drag_model_part_name = drag_adjoint_response_function_settings["custom_settings"]["structure_model_part_name"].GetString()
        self.drag_direction = drag_adjoint_response_function_settings["custom_settings"]["drag_direction"].GetVector()

        # check for model part name settings
        primal_project_parameters_file = self.problem_setup_folder / self.problem_setup_file_settings["primal_project_parameters_file"].GetString()
        with open(primal_project_parameters_file, "r") as file_input:
            primal_settings = Kratos.Parameters(file_input.read())

        self.mdpa_name = primal_settings["solver_settings"]["model_import_settings"]["input_filename"].GetString()

        self.main_model_part_name = primal_settings["solver_settings"]["model_part_name"].GetString()
        self.lift_model_part_name =  self.main_model_part_name + "." + self.lift_model_part_name
        self.drag_model_part_name =  self.main_model_part_name + "." + self.drag_model_part_name

    def Initialize(self):
        pass

    def UpdateDesign(self, updated_model_part, variable):
        self.updated_model_part = updated_model_part

    def InitializeSolutionStep(self):
        self.lift = None
        self.drag = None

        startTime = timer.time()

        # copy data to the working directory
        RecursiveCopy(str(self.problem_setup_folder), ".")

        # write the new shape mdpa
        Kratos.ModelPartIO(self.mdpa_name, Kratos.IO.WRITE | Kratos.IO.MESH_ONLY).WriteModelPart(self.updated_model_part)

        Kratos.Logger.PrintInfo(self._GetLabel(), "Time needed to copy data = ",round(timer.time() - startTime,2),"s")

    def CalculateValue(self):
        startTime = timer.time()

        # open primal settings
        primal_settings_file_name = self.problem_setup_file_settings["primal_project_parameters_file"].GetString()
        primal_parameters = SolvePrimalProblem(primal_settings_file_name)

        # calculate lift and drag
        self.lift = CalculateTimeAveragedDrag(primal_parameters, self.lift_model_part_name, self.lift_direction)
        self.drag = CalculateTimeAveragedDrag(primal_parameters, self.drag_model_part_name, self.drag_direction)

        Kratos.Logger.PrintInfo(self._GetLabel(), "Time needed for calculating the response value = ",round(timer.time() - startTime,2),"s")

    def CalculateGradient(self):
        # solve adjoint lift problem
        lift_start_time = timer.time()

        self.lift_adjoint_model = Kratos.Model()
        self.lift_adjoint_simulation = SolveAdjointProblem(
                self.lift_adjoint_model,
                self.problem_setup_file_settings["lift_adjoint_project_parameters_file"].GetString(),
                "lift_adjoint_evaluation.log")

        Kratos.Logger.PrintInfo(self._GetLabel(), "Time needed for solving the lift adjoint analysis = ",round(timer.time() - lift_start_time,2),"s")

        # solve adjoint drag problem
        drag_start_time = timer.time()

        self.drag_adjoint_model = Kratos.Model()
        self.drag_adjoint_simulation = SolveAdjointProblem(
                self.drag_adjoint_model,
                self.problem_setup_file_settings["drag_adjoint_project_parameters_file"].GetString(),
                "drag_adjoint_evaluation.log")

        Kratos.Logger.PrintInfo(self._GetLabel(), "Time needed for solving the drag adjoint analysis = ",round(timer.time() - drag_start_time,2),"s")
        Kratos.Logger.PrintInfo(self._GetLabel(), "Time needed for solving the total adjoint analysis = ",round(timer.time() - lift_start_time,2),"s")

    def GetValue(self):
        if (self.lift is None or self.drag is None):
            raise RuntimeError("GetValue: Lift or Drag is not computed. Please use CalculateValue method first.")
        return self.lift / self.drag

    def GetNodalGradient(self, variable):
        if (self.lift is None or self.drag is None):
            raise RuntimeError("GetNodalGradient: Lift or Drag is not computed. Please use CalculateValue method first.")
        if variable != Kratos.SHAPE_SENSITIVITY:
            raise RuntimeError("GetNodalGradient: No gradient for {}!".format(variable.Name))

        gradient = {}
        lift_adjoint_model_part = self.lift_adjoint_model[self.main_model_part_name]
        drag_adjoint_model_part = self.drag_adjoint_model[self.main_model_part_name]
        for lift_node, drag_node in zip(lift_adjoint_model_part.Nodes, drag_adjoint_model_part.Nodes):
            node_id = lift_node.Id
            if (lift_node.Id != drag_node.Id):
                raise RuntimeError("Node mismatch")
            gradient[node_id] = \
                lift_node.GetSolutionStepValue(variable) / self.drag - \
                self.lift * drag_node.GetSolutionStepValue(variable) / (self.drag ** 2)
        return gradient

    def IsEvaluatedInFolder(self):
        return True

    @staticmethod
    def _GetLabel():
        return "AdjointLiftToDrag"

