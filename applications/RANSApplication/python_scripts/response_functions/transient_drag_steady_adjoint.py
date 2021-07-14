"""
This module contains an interface to the available response functions
"""
import time as timer

import KratosMultiphysics as Kratos

from KratosMultiphysics.RANSApplication.response_functions.utilities import SolvePrimalProblem
from KratosMultiphysics.RANSApplication.response_functions.utilities import SolveAdjointProblem
from KratosMultiphysics.RANSApplication.response_functions.utilities import CalculateTimeAveragedDrag

from KratosMultiphysics.RANSApplication.response_functions.drag import Drag
from KratosMultiphysics.RANSApplication.response_functions.drag import DragConfigurationPrimalFromFile
from KratosMultiphysics.RANSApplication.response_functions.drag import DragConfigurationPrimalFromExistingSolution


class TransientDragSteadyAdjoint(Drag):
    def __init__(self, identifier, response_settings, model):
        self.identifier = identifier
        self.response_settings = response_settings

        default_parameters = Kratos.Parameters("""
        {
            "response_type"       : "drag",
            "problem_setup_settings" : {
                "primal_problem_solving_type"     : "solved"
            },
            "primal_transient_project_parameters_file": "PLEASE_SPECIFY_PRIMAL_TRANSIENT_PROJECT_PARAMETERS_FILE_NAME",
            "primal_steady_project_parameters_file"   : "PLEASE_SPECIFY_PRIMAL_steady_PROJECT_PARAMETERS_FILE_NAME",
            "statistics_start_point_control_value"    : 0.0
        }
        """)
        self.identifier = identifier
        self.response_settings.ValidateAndAssignDefaults(default_parameters)

        self.problem_setup_settings = self.response_settings["problem_setup_settings"]
        self.primal_problem_solving_type = self.problem_setup_settings["primal_problem_solving_type"].GetString()
        self.statistics_start_point_control_value = self.response_settings["statistics_start_point_control_value"].GetDouble()

        if (self.primal_problem_solving_type == "solved"):
            self.drag_configuration = DragConfigurationPrimalFromFile(self.problem_setup_settings)
            self.is_evaluated_in_folder = True
        elif (self.primal_problem_solving_type == "read"):
            self.drag_configuration = DragConfigurationPrimalFromExistingSolution(self.problem_setup_settings)
            self.is_evaluated_in_folder = False
        else:
            raise RuntimeError("Unsupported \"primal_problem_solving_type\" requested. Requested type = \"{}\". Supported types are:\n\t \"solved\"\n\t\"read\"\n")

    def CalculateValue(self):
        startTime = timer.time()

        # open primal settings
        primal_parameters = SolvePrimalProblem(self.response_settings["primal_transient_project_parameters_file"].GetString(), "primal_transient_evaluation.log")

        # calculate drag
        self.drag = CalculateTimeAveragedDrag(primal_parameters, self.drag_configuration.GetDragModelPartName(), self.drag_configuration.GetDragDirection(), self.statistics_start_point_control_value)

        Kratos.Logger.PrintInfo(self._GetLabel(), "Time needed for calculating the response value = ",round(timer.time() - startTime,2),"s")

    def CalculateGradient(self):
        # solve adjoint drag problem
        drag_start_time = timer.time()

        self.gradient = {}

        # now run the RANS analysis
        _ = SolvePrimalProblem(self.response_settings["primal_steady_project_parameters_file"].GetString(), "primal_steady_evaluation.log")

        drag_adjoint_model = Kratos.Model()
        _ = SolveAdjointProblem(
                drag_adjoint_model,
                self.drag_configuration.GetAdjointProjectParametersFileName(),
                self.identifier + ".log")

        drag_adjoint_model_part = drag_adjoint_model[self.drag_configuration.GetMainModelPartName()]
        for drag_node in drag_adjoint_model_part.Nodes:
            self.gradient[drag_node.Id] = drag_node.GetSolutionStepValue(Kratos.SHAPE_SENSITIVITY)

        Kratos.Logger.PrintInfo(self._GetLabel(), "Time needed for solving the drag adjoint analysis = ",round(timer.time() - drag_start_time,2),"s")
