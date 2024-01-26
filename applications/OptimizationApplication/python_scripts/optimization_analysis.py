import KratosMultiphysics as Kratos
from KratosMultiphysics.process_factory import KratosProcessFactory
from KratosMultiphysics.OptimizationApplication.algorithms.algorithm import Algorithm
from KratosMultiphysics.OptimizationApplication.model_part_controllers.model_part_controller import ModelPartController
from KratosMultiphysics.OptimizationApplication.execution_policies.execution_policy_decorator import ExecutionPolicyDecorator
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.controls.control import Control
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import OptimizationComponentFactory
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import CallOnAll
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import OptimizationAnalysisTimeLogger

class OptimizationAnalysis:
    @classmethod
    def GetDefaultParameters(cls) -> Kratos.Parameters:
        return Kratos.Parameters("""{
            "problem_data"      : {},
            "model_parts"       : [],
            "analyses"          : [],
            "responses"         : [],
            "controls"          : [],
            "algorithm_settings": {},
            "processes"  : {
                "kratos_processes"           : {},
                "optimization_data_processes": {}
            }
        }""")

    def __init__(self, model: Kratos.Model, project_parameters: Kratos.Parameters):
        self.model = model
        self.project_parameters = project_parameters
        self.project_parameters.ValidateAndAssignDefaults(self.GetDefaultParameters())

        self.optimization_problem = OptimizationProblem(self.project_parameters["problem_data"]["echo_level"].GetInt())

        self.__list_of_model_part_controllers: 'list[ModelPartController]' = []
        self.__algorithm: Algorithm = None

        self._CreateModelPartControllers()
        self._CreateAnalyses()
        self._CreateControls()
        self._CreateResponses()
        self._CreateAlgorithm()
        self._CreateProcesses()

    def Initialize(self):
        CallOnAll(self.__list_of_model_part_controllers, ModelPartController.ImportModelPart)
        CallOnAll(self.__list_of_model_part_controllers, ModelPartController.Initialize)
        for process_type in self.__algorithm.GetProcessesOrder():
            CallOnAll(self.optimization_problem.GetListOfProcesses(process_type), Kratos.Process.ExecuteInitialize)
        CallOnAll(self.optimization_problem.GetListOfExecutionPolicies(), ExecutionPolicyDecorator.Initialize)

        self.__algorithm.Initialize()

    def Check(self):
        for process_type in self.__algorithm.GetProcessesOrder():
            CallOnAll(self.optimization_problem.GetListOfProcesses(process_type), Kratos.Process.Check)
        CallOnAll(self.optimization_problem.GetListOfExecutionPolicies(), ExecutionPolicyDecorator.Check)

        self.__algorithm.Check()

    def Finalize(self):
        self.__algorithm.Finalize()

        CallOnAll(self.__list_of_model_part_controllers, ModelPartController.Finalize)
        for process_type in self.__algorithm.GetProcessesOrder():
            CallOnAll(self.optimization_problem.GetListOfProcesses(process_type), Kratos.Process.ExecuteFinalize)
        CallOnAll(self.optimization_problem.GetListOfExecutionPolicies(), ExecutionPolicyDecorator.Finalize)

    def Run(self):
        with OptimizationAnalysisTimeLogger():
            self.Initialize()
            self.Check()
            self.__algorithm.Solve()
            self.Finalize()

    def _CreateModelPartControllers(self):
        default_settings = Kratos.Parameters("""{
            "type": "mdpa_model_part_controller",
            "module": "KratosMultiphysics.OptimizationApplication.model_part_controllers"
        }""")
        for model_part_controller_settings in self.project_parameters["model_parts"]:
            model_part_controller_settings.AddMissingParameters(default_settings)
            model_part_controller: ModelPartController = OptimizationComponentFactory(self.model, model_part_controller_settings, self.optimization_problem)
            self.__list_of_model_part_controllers.append(model_part_controller)

    def _CreateAnalyses(self):
        default_settings = Kratos.Parameters("""{
            "module": "KratosMultiphysics.OptimizationApplication.execution_policies"
        }""")
        for analyses_settings in self.project_parameters["analyses"]:
            analyses_settings.AddMissingParameters(default_settings)
            execution_policy = OptimizationComponentFactory(self.model, analyses_settings, self.optimization_problem)
            self.optimization_problem.AddComponent(execution_policy)

    def _CreateResponses(self):
        default_settings = Kratos.Parameters("""{
            "module": "KratosMultiphysics.OptimizationApplication.responses"
        }""")
        for response_settings in self.project_parameters["responses"]:
            response_settings.AddMissingParameters(default_settings)
            response_function: ResponseFunction = OptimizationComponentFactory(self.model, response_settings, self.optimization_problem)
            self.optimization_problem.AddComponent(response_function)

    def _CreateControls(self):
        default_settings = Kratos.Parameters("""{
            "module" : "KratosMultiphysics.OptimizationApplication.controls"
        }""")
        for control_settings in self.project_parameters["controls"]:
            control_settings.AddMissingParameters(default_settings)
            control = OptimizationComponentFactory(self.model, control_settings, self.optimization_problem)
            self.optimization_problem.AddComponent(control)

    def _CreateProcesses(self):
        process_settings = self.project_parameters["processes"]
        process_settings.ValidateAndAssignDefaults(self.GetDefaultParameters()["processes"])

        kratos_processes = process_settings["kratos_processes"]
        optimization_data_processes = process_settings["optimization_data_processes"]

        factory = KratosProcessFactory(self.model)

        optimization_data_process_default_settings = Kratos.Parameters("""{
            "module" : "KratosMultiphysics.OptimizationApplication.optimization_data_processes"
        }""")

        for process_type in self.__algorithm.GetProcessesOrder():
            self.optimization_problem.AddProcessType(process_type)
            if kratos_processes.Has(process_type):
                for process in factory.ConstructListOfProcesses(kratos_processes[process_type]):
                    self.optimization_problem.AddProcess(process_type, process)
            if optimization_data_processes.Has(process_type):
                for process_settings in optimization_data_processes[process_type]:
                    process_settings.AddMissingParameters(optimization_data_process_default_settings)
                    process = OptimizationComponentFactory(self.model, process_settings, self.optimization_problem)
                    self.optimization_problem.AddProcess(process_type, process)

    def _CreateAlgorithm(self):
        default_settings = Kratos.Parameters("""{
            "module" : "KratosMultiphysics.OptimizationApplication.algorithms"
        }""")
        algorithm_settings = self.project_parameters["algorithm_settings"]
        algorithm_settings.AddMissingParameters(default_settings)
        self.__algorithm = OptimizationComponentFactory(self.model, algorithm_settings, self.optimization_problem)
