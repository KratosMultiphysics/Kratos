"""
This module contains an interface to the available response functions
"""
import time as timer
import shutil
from pathlib import Path

import KratosMultiphysics as Kratos
from KratosMultiphysics.response_functions.response_function_interface import ResponseFunctionInterface

from KratosMultiphysics.RANSApplication.response_functions.utilities import SolvePrimalProblem
from KratosMultiphysics.RANSApplication.response_functions.utilities import SolveAdjointProblem
from KratosMultiphysics.RANSApplication.response_functions.utilities import CalculateTimeAveragedResponseValue
from KratosMultiphysics.RANSApplication.response_functions.utilities import RecursiveCopy

class DomainIntegrated3DVectorMagnitudeSquarePowerMean(ResponseFunctionInterface):
    def __init__(self, identifier, response_settings, model):
        self.identifier = identifier
        self.response_settings = response_settings

        default_parameters = Kratos.Parameters("""
        {
            "response_type"       : "domain_integrated_3d_vector_magnitude_square_power_mean",
            "problem_setup_folder": "PLEASE_SPECIFY_PROBLEM_SETUP_FOLDER",
            "problem_setup_files" : {
                "primal_project_parameters_file" : "PLEASE_SPECIFY_PRIMAL_PROJECT_PARAMETERS_FILE",
                "adjoint_project_parameters_file": "PLEASE_SPECIFY_ADJOINT_PROJECT_PARAMETERS_FILE"
            },
            "clean_primal_solution": false,
            "primal_solution_folder_name": "PLEASE_SPECIFY_PRIMAL_SOLUTION_FOLDER_NAME"
        }
        """)

        self.response_settings.ValidateAndAssignDefaults(default_parameters)

        self.problem_setup_folder = Path(self.response_settings["problem_setup_folder"].GetString()).absolute()
        if (not self.problem_setup_folder.is_dir()):
            raise Exception("The provided \"problem_setup_folder\" is not a directory [ \"problem_setup_folder\" = {:s} ].".format(str(self.problem_setup_folder)))

        self.problem_setup_file_settings = self.response_settings["problem_setup_files"]

        # checks for lift response settings
        adjoint_project_parameters_file = self.problem_setup_folder / self.problem_setup_file_settings["adjoint_project_parameters_file"].GetString()
        with open(adjoint_project_parameters_file, "r") as file_input:
            adjoint_settings = Kratos.Parameters(file_input.read())

        self.adjoint_response_function_settings = adjoint_settings["solver_settings"]["response_function_settings"]
        response_type = self.adjoint_response_function_settings["response_type"].GetString()
        if (response_type != "domain_integrated_3d_vector_magnitude_square_power_mean"):
            raise Exception("DomainIntegrated3DVectorMagnitudeSquarePowerMeanneeds \"response_type\" in \"response_function_settings\" of \"adjoint_project_parameters_file\" file to be domain_integrated_3d_vector_magnitude_square_power_mean. [ \"adjoint_project_parameters_file\" = {:s}, \"response_type\" = {:s} ].".format(
                adjoint_project_parameters_file, response_type))

        self.model_part_name = self.adjoint_response_function_settings["custom_settings"]["model_part_name"].GetString()

        # check for model part name settings
        primal_project_parameters_file = self.problem_setup_folder / self.problem_setup_file_settings["primal_project_parameters_file"].GetString()
        with open(primal_project_parameters_file, "r") as file_input:
            primal_settings = Kratos.Parameters(file_input.read())

        self.mdpa_name = primal_settings["solver_settings"]["model_import_settings"]["input_filename"].GetString()

        self.main_model_part_name = primal_settings["solver_settings"]["model_part_name"].GetString()
        self.model_part_name =  self.main_model_part_name + "." + self.model_part_name

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
        self.value = CalculateTimeAveragedResponseValue(primal_parameters, self.model_part_name, self.adjoint_response_function_settings)

        Kratos.Logger.PrintInfo(self._GetLabel(), "Time needed for calculating the response value = ",round(timer.time() - startTime,2),"s")

    def CalculateGradient(self):
        # solve adjoint lift problem
        start_time = timer.time()

        self.adjoint_model = Kratos.Model()
        self.adjoint_simulation = SolveAdjointProblem(
                self.adjoint_model,
                self.problem_setup_file_settings["adjoint_project_parameters_file"].GetString(),
                "adjoint_evaluation.log")

        if self.response_settings["clean_primal_solution"].GetBool():
            shutil.rmtree(self.response_settings["primal_solution_folder_name"].GetString())

        Kratos.Logger.PrintInfo(self._GetLabel(), "Time needed for solving the adjoint analysis = ",round(timer.time() - start_time,2),"s")

    def GetValue(self):
        if (self.value is None):
            raise RuntimeError("GetValue: Response value is not computed. Please use CalculateValue method first.")
        return self.value

    def GetNodalGradient(self, variable):
        if variable != Kratos.SHAPE_SENSITIVITY:
            raise RuntimeError("GetNodalGradient: No gradient for {}!".format(variable.Name))

        gradient = {}
        adjoint_model_part = self.adjoint_model[self.main_model_part_name]
        for node in adjoint_model_part.Nodes:
            gradient[node.Id] = node.GetSolutionStepValue(variable)
        return gradient

    def IsEvaluatedInFolder(self):
        return True

    @staticmethod
    def _GetLabel():
        return "AdjointDomainIntegrated3DVectorMagnitudeSquarePowerMean"

