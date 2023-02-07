import KratosMultiphysics as Kratos
from KratosMultiphysics.analysis_stage import AnalysisStage
from KratosMultiphysics.OptimizationApplication.algorithms.algorithm import Algorithm
from KratosMultiphysics.OptimizationApplication.model_part_controllers.model_part_controller import ModelPartController
from KratosMultiphysics.OptimizationApplication.execution_policies.execution_policy_wrapper import ExecutionPolicyWrapper
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.utilities.control_transformation_technique import ControlTransformationTechnique
from KratosMultiphysics.OptimizationApplication.optimization_info import OptimizationInfo
from KratosMultiphysics.OptimizationApplication.optimization_routine import OptimizationRoutine
from KratosMultiphysics.OptimizationApplication.utilities.helper_utils import CallOnAll
from KratosMultiphysics.OptimizationApplication.utilities.helper_utils import OptimizationRoutineFactory

class OptimizationAnalysis(AnalysisStage):
    def __init__(self, model: Kratos.Model, project_parameters: Kratos.Parameters):
        default_parametetrs = Kratos.Parameters("""{
            "model_parts"       : [],
            "analyses"          : [],
            "responses"         : [],
            "controls"          : [],
            "algorithm_settings": {}
        }""")

        self.model = model
        self.project_parameters = project_parameters
        self.project_parameters.AddMissingParameters(default_parametetrs)

        self.optimization_info = OptimizationInfo(self.project_parameters["problem_data"]["echo_level"].GetInt())

        self._CreateModelPartControllers()
        self._CreateAnalyses()
        self._CreateResponses()
        self._CreateControlTechniques()

        super().__init__(model, project_parameters)

    def Initialize(self):
        super().Initialize()
        CallOnAll(self._GetOrderedOptimizationRoutines(), OptimizationRoutine.Initialize)

    def InitializeSolutionStep(self):
        CallOnAll(self._GetOrderedOptimizationRoutines(), OptimizationRoutine.InitializeSolutionStep)
        super().InitializeSolutionStep()

    def FinalizeSolutionStep(self):
        CallOnAll(self._GetOrderedOptimizationRoutines(), OptimizationRoutine.FinalizeSolutionStep)
        super().FinalizeSolutionStep()

    def Finalize(self):
        CallOnAll(self._GetOrderedOptimizationRoutines(), OptimizationRoutine.Finalize)
        return super().Finalize()

    def _CreateSolver(self):
        default_algorithm_settings = Kratos.Parameters("""{
            "module"            : "KratosMultiphysics.OptimizationApplication.algorithms",
            "type"              : "PLEASE_PROVIDE_AN_ALGORITHM_CLASS_NAME"
        }""")
        algorithm_settings = self.project_parameters["algorithm_settings"]
        algorithm_settings.AddMissingParameters(default_algorithm_settings)

        return OptimizationRoutineFactory(algorithm_settings["module"].GetString(), algorithm_settings["type"].GetString(), self.model, algorithm_settings, self.optimization_info, Algorithm)

    def _GetSimulationName(self):
        return "OptimizationAnalysis"

    def PrintAnalysisStageProgressInformation(self):
        Kratos.Logger.PrintInfo(self._GetSimulationName(), "STEP: ", self.optimization_info["step"])

    def KeepAdvancingSolutionLoop(self):
        return self.optimization_info["step"] < self.end_time and not self._GetSolver().IsConverged()

    def _CreateModelPartControllers(self):
        for model_part_controller_settings in self.project_parameters["model_parts"]:
            default_model_part_controller_settings = Kratos.Parameters("""{
                "name"    : "",
                "module"  : "KratosMultiphysics.OptimizationApplication.model_part_controllers",
                "type"    : "MdpaModelPartController",
                "settings": {}
            }""")
            model_part_controller_settings.ValidateAndAssignDefaults(default_model_part_controller_settings)
            routine: ModelPartController = OptimizationRoutineFactory(model_part_controller_settings["module"].GetString(), model_part_controller_settings["type"].GetString(), self.model, model_part_controller_settings["settings"], self.optimization_info, ModelPartController)
            self.optimization_info.AddOptimizationRoutine(ModelPartController, model_part_controller_settings["name"].GetString(), routine)

    def _CreateAnalyses(self):
        for analyses_settings in self.project_parameters["analyses"]:
            routine = ExecutionPolicyWrapper(self.model, analyses_settings)
            self.optimization_info.AddOptimizationRoutine(ExecutionPolicyWrapper, routine.GetName(), routine)

    def _CreateResponses(self):
        default_settings = Kratos.Parameters("""{
            "name"     : "",
            "module"   : "KratosMultiphysics.OptimizationApplication.responses",
            "type"     : "",
            "settings" : {}
        }""")
        for response_settings in self.project_parameters["responses"]:
            response_settings.ValidateAndAssignDefaults(default_settings)
            routine: ResponseFunction = OptimizationRoutineFactory(response_settings["module"].GetString(), response_settings["type"].GetString(), self.model, response_settings["settings"], self.optimization_info, ResponseFunction)
            self.optimization_info.AddOptimizationRoutine(ResponseFunction, response_settings["name"].GetString(), routine)

    def _CreateControlTechniques(self):
        for control_settings in self.project_parameters["controls"]:
            routine = ControlTransformationTechnique(self.model, control_settings, self.optimization_info)
            self.optimization_info.AddOptimizationRoutine(ControlTransformationTechnique, routine.GetName(), routine)

    def _GetOrderedOptimizationRoutines(self):
        routines = []
        routines.extend(self.optimization_info.GetOptimizationRoutines(ModelPartController))
        routines.extend(self.optimization_info.GetOptimizationRoutines(ExecutionPolicyWrapper))
        routines.extend(self.optimization_info.GetOptimizationRoutines(ResponseFunction))
        routines.extend(self.optimization_info.GetOptimizationRoutines(ControlTransformationTechnique))
        return routines
