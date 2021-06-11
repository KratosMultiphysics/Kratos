"""
This module contains an interface to the available response functions
"""
import time as timer
import os, shutil
from pathlib import Path

import KratosMultiphysics as Kratos
from KratosMultiphysics.response_functions.response_function_interface import ResponseFunctionInterface

from KratosMultiphysics.RANSApplication.response_functions.utilities import SolvePrimalProblem
from KratosMultiphysics.RANSApplication.response_functions.utilities import SolveAdjointProblem
from KratosMultiphysics.RANSApplication.response_functions.utilities import CalculateTimeAveragedDrag
from KratosMultiphysics.RANSApplication.response_functions.utilities import RecursiveCopy

class DragConfiguration:
    def __init__(self, parameters):
        default_parameters = Kratos.Parameters("""{
            "primal_problem_solving_type"     : "solved",
            "problem_setup_folder"            : "PLEASE_SPECIFY_PROBLEM_SETUP_FOLDER",
            "primal_project_parameters_file"  : "PLEASE_SPECIFY_PRIMAL_PROJECT_PARAMETERS_FILE",
            "adjoint_project_parameters_file" : "PLEASE_SPECIFY_LIFT_ADJOINT_PROJECT_PARAMETERS_FILE",
            "primal_problem_location_variable": "RANS_PRIMAL_SOLUTION_LOCATION_1"
        }""")

        self.parameters = parameters
        self.parameters.ValidateAndAssignDefaults(default_parameters)

        self.problem_setup_folder = Path(self.parameters["problem_setup_folder"].GetString()).absolute()
        if (not self.problem_setup_folder.is_dir()):
            raise Exception("The provided \"problem_setup_folder\" is not a directory [ \"problem_setup_folder\" = {:s} ].".format(str(self.problem_setup_folder)))

        # check for drag response settings
        adjoint_project_parameters_file = self.problem_setup_folder / self.parameters["adjoint_project_parameters_file"].GetString()
        with open(adjoint_project_parameters_file, "r") as file_input:
            drag_adjoint_settings = Kratos.Parameters(file_input.read())

        drag_adjoint_response_function_settings = drag_adjoint_settings["solver_settings"]["response_function_settings"]
        drag_response_type = drag_adjoint_response_function_settings["response_type"].GetString()
        if (drag_response_type != "drag"):
            raise Exception("Lift to drag needs \"response_type\" in \"response_function_settings\" of \"adjoint_project_parameters_file\" file to be drag. [ \"adjoint_project_parameters_file\" = {:s}, \"response_type\" = {:s} ].".format(
                adjoint_project_parameters_file, drag_response_type))

        self.drag_model_part_name = drag_adjoint_response_function_settings["custom_settings"]["structure_model_part_name"].GetString()
        self.drag_direction = drag_adjoint_response_function_settings["custom_settings"]["drag_direction"].GetVector()

        # check for model part name settings
        primal_project_parameters_file = self.problem_setup_folder / self.parameters["primal_project_parameters_file"].GetString()
        with open(primal_project_parameters_file, "r") as file_input:
            primal_settings = Kratos.Parameters(file_input.read())

        self.mdpa_name = primal_settings["solver_settings"]["model_import_settings"]["input_filename"].GetString()

        self.main_model_part_name = primal_settings["solver_settings"]["model_part_name"].GetString()
        self.drag_model_part_name =  self.main_model_part_name + "." + self.drag_model_part_name
        self.primal_problem_location_variable = Kratos.KratosGlobals.GetVariable(self.parameters["primal_problem_location_variable"].GetString())

    def GetPrimalProjectParametersFileName(self):
        return self.parameters["primal_project_parameters_file"].GetString()

    def GetAdjointProjectParametersFileName(self):
        return self.parameters["adjoint_project_parameters_file"].GetString()

    def GetMainModelPartName(self):
        return self.main_model_part_name

    def GetDragModelPartName(self):
        return self.drag_model_part_name

    def GetDragDirection(self):
        return self.drag_direction

class DragConfigurationPrimalFromFile(DragConfiguration):
    def __init__(self, parameters):
        super().__init__(parameters)


    def InitializeSolutionStep(self, updated_model_part):
        # copy data to the working directory
        RecursiveCopy(str(self.problem_setup_folder), ".")

        # write the new shape mdpa
        Kratos.ModelPartIO(self.mdpa_name, Kratos.IO.WRITE | Kratos.IO.MESH_ONLY).WriteModelPart(updated_model_part)

        updated_model_part.SetValue(self.primal_problem_location_variable, str(Path.cwd().absolute()))

        Kratos.Logger.PrintInfo("DragConfigurationPrimalFromFile", "Copied data successfully.")

    def FinalizeSolutionStep(self):
        pass

class DragConfigurationPrimalFromExistingSolution(DragConfiguration):
    def __init__(self, parameters):
        super().__init__(parameters)

    def InitializeSolutionStep(self, updated_model_part):
        self.current_path = Path.cwd().absolute()
        if (updated_model_part.Has(self.primal_problem_location_variable)):
            os.chdir(updated_model_part.GetValue(self.primal_problem_location_variable))
        else:
            raise RuntimeError("Primal evaluation is not found.")

        Kratos.Logger.PrintInfo("DragConfigurationPrimalFromExistingSolution", "Initialized successfully.")

    def FinalizeSolutionStep(self):
        os.chdir(self.current_path)

class Drag(ResponseFunctionInterface):
    def __init__(self, identifier, response_settings, model):
        self.identifier = identifier
        self.response_settings = response_settings

        default_parameters = Kratos.Parameters("""
        {
            "response_type"       : "drag",
            "problem_setup_settings" : {
                "primal_problem_solving_type"     : "solved"
            },
            "clean_primal_solution": false,
            "primal_solution_folder_name": "PLEASE_SPECIFY_PRIMAL_SOLUTION_FOLDER_NAME",
            "start_time": 0.0
        }
        """)
        self.identifier = identifier
        self.response_settings.ValidateAndAssignDefaults(default_parameters)

        self.problem_setup_settings = self.response_settings["problem_setup_settings"]
        self.primal_problem_solving_type = self.problem_setup_settings["primal_problem_solving_type"].GetString()

        if (self.primal_problem_solving_type == "solved"):
            self.drag_configuration = DragConfigurationPrimalFromFile(self.problem_setup_settings)
            self.is_evaluated_in_folder = True
        elif (self.primal_problem_solving_type == "read"):
            self.drag_configuration = DragConfigurationPrimalFromExistingSolution(self.problem_setup_settings)
            self.is_evaluated_in_folder = False
        else:
            raise RuntimeError("Unsupported \"primal_problem_solving_type\" requested. Requested type = \"{}\". Supported types are:\n\t \"solved\"\n\t\"read\"\n")

    def Initialize(self):
        pass

    def UpdateDesign(self, updated_model_part, variable):
        self.updated_model_part = updated_model_part

    def InitializeSolutionStep(self):
        self.drag = None

        startTime = timer.time()

        self.drag_configuration.InitializeSolutionStep(self.updated_model_part)

        Kratos.Logger.PrintInfo(self._GetLabel(), "Time needed for initialization = ",round(timer.time() - startTime,2),"s")

    def CalculateValue(self):
        startTime = timer.time()

        # open primal settings
        primal_parameters = SolvePrimalProblem(self.drag_configuration.GetPrimalProjectParametersFileName())

        # calculate drag
        self.drag = CalculateTimeAveragedDrag(primal_parameters, self.drag_configuration.GetDragModelPartName(), self.drag_configuration.GetDragDirection(), self.response_settings["start_time"].GetDouble())

        Kratos.Logger.PrintInfo(self._GetLabel(), "Time needed for calculating the response value = ",round(timer.time() - startTime,2),"s")

    def CalculateGradient(self):
        # solve adjoint drag problem
        drag_start_time = timer.time()

        self.gradient = {}

        drag_adjoint_model = Kratos.Model()
        _ = SolveAdjointProblem(
                drag_adjoint_model,
                self.drag_configuration.GetAdjointProjectParametersFileName(),
                self.identifier + ".log")

        drag_adjoint_model_part = drag_adjoint_model[self.drag_configuration.GetMainModelPartName()]
        for drag_node in drag_adjoint_model_part.Nodes:
            self.gradient[drag_node.Id] = drag_node.GetSolutionStepValue(Kratos.SHAPE_SENSITIVITY)

        if self.response_settings["clean_primal_solution"].GetBool():
            shutil.rmtree(self.response_settings["primal_solution_folder_name"].GetString())

        Kratos.Logger.PrintInfo(self._GetLabel(), "Time needed for solving the drag adjoint analysis = ",round(timer.time() - drag_start_time,2),"s")

    def GetValue(self):
        if (self.drag is None):
            raise RuntimeError("GetValue: Drag is not computed. Please use CalculateValue method first.")
        return self.drag

    def GetNodalGradient(self, variable):
        if (self.drag is None):
            raise RuntimeError("GetNodalGradient: Drag is not computed. Please use CalculateValue method first.")
        if variable != Kratos.SHAPE_SENSITIVITY:
            raise RuntimeError("GetNodalGradient: No gradient for {}!".format(variable.Name))

        return self.gradient

    def FinalizeSolutionStep(self):
        self.drag_configuration.FinalizeSolutionStep()

    def IsEvaluatedInFolder(self):
        return self.is_evaluated_in_folder

    @staticmethod
    def _GetLabel():
        return "AdjointDrag"

