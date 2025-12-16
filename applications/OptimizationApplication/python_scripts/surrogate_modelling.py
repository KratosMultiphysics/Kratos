import KratosMultiphysics as Kratos
from KratosMultiphysics.process_factory import KratosProcessFactory
from KratosMultiphysics.OptimizationApplication.algorithms.algorithm import Algorithm
from KratosMultiphysics.OptimizationApplication.model_part_controllers.model_part_controller import ModelPartController
from KratosMultiphysics.OptimizationApplication.execution_policies.execution_policy_decorator import ExecutionPolicyDecorator
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.controls.control import Control
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import SurrogateModellingComponentFactory
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import CallOnAll
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import OptimizationAnalysisTimeLogger

class FsiSurrogate:
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

    def __init__(self, models: 'dict[Kratos.Model]', project_parameters: Kratos.Parameters):
        self.models = models
        self.project_parameters = project_parameters
        self.project_parameters.ValidateAndAssignDefaults(self.GetDefaultParameters())
        self.optimization_problem = OptimizationProblem(self.project_parameters["problem_data"]["echo_level"].GetInt())
        self.__algorithm: Algorithm = None
        self._CreateAnalyses()
        self._CreateControls()
        self._CreateAlgorithm()
        self._CreateProcesses()

    def Initialize(self):
        CallOnAll(self.optimization_problem.GetListOfExecutionPolicies(), ExecutionPolicyDecorator.Initialize)
        CallOnAll(self.optimization_problem.GetListOfControls(), Control.Initialize)
        for process_type in self.__algorithm.GetProcessesOrder():
            CallOnAll(self.optimization_problem.GetListOfProcesses(process_type), Kratos.Process.ExecuteInitialize)
        self.__algorithm.Initialize()

    def Check(self):
        for process_type in self.__algorithm.GetProcessesOrder():
            CallOnAll(self.optimization_problem.GetListOfProcesses(process_type), Kratos.Process.Check)
        CallOnAll(self.optimization_problem.GetListOfExecutionPolicies(), ExecutionPolicyDecorator.Check)
        CallOnAll(self.optimization_problem.GetListOfControls(), Control.Check)
        CallOnAll(self.optimization_problem.GetListOfResponses(), ResponseFunction.Check)

        self.__algorithm.Check()

    def Finalize(self):
        self.__algorithm.Finalize()

        CallOnAll(self.__list_of_model_part_controllers, ModelPartController.Finalize)
        for process_type in self.__algorithm.GetProcessesOrder():
            CallOnAll(self.optimization_problem.GetListOfProcesses(process_type), Kratos.Process.ExecuteFinalize)
        CallOnAll(self.optimization_problem.GetListOfExecutionPolicies(), ExecutionPolicyDecorator.Finalize)
        CallOnAll(self.optimization_problem.GetListOfControls(), Control.Finalize)
        CallOnAll(self.optimization_problem.GetListOfResponses(), ResponseFunction.Finalize)

    def Run(self):
        with OptimizationAnalysisTimeLogger():
            self.Initialize()
            #self.Check()
            self.__algorithm.Solve()
            #self.Finalize()

    def _CreateAnalyses(self):
        default_settings = Kratos.Parameters("""{
            "module": "KratosMultiphysics.OptimizationApplication.execution_policies"
        }""")
        for analyses_settings in self.project_parameters["analyses"]:
            analyses_settings.AddMissingParameters(default_settings)
            execution_policy = SurrogateModellingComponentFactory(self.models, analyses_settings, self.optimization_problem)
            self.optimization_problem.AddComponent(execution_policy)

    def _CreateControls(self):
        default_settings = Kratos.Parameters("""{
            "module" : "KratosMultiphysics.OptimizationApplication.controls"
        }""")
        for control_settings in self.project_parameters["controls"].values():
            control_settings.AddMissingParameters(default_settings)
            control = SurrogateModellingComponentFactory(self.models, control_settings, self.optimization_problem)
            self.optimization_problem.AddComponent(control)

    def _CreateProcesses(self):
        process_settings = self.project_parameters["processes"]
        process_settings.ValidateAndAssignDefaults(self.GetDefaultParameters()["processes"])
        
        optimization_data_processes = process_settings["optimization_data_processes"]
        optimization_data_process_default_settings = Kratos.Parameters("""{
            "module" : "KratosMultiphysics.OptimizationApplication.optimization_data_processes"
        }""")

        for process_type in self.__algorithm.GetProcessesOrder():
            self.optimization_problem.AddProcessType(process_type)
            if optimization_data_processes.Has(process_type):
                for process_settings in optimization_data_processes[process_type]:
                    process_settings.AddMissingParameters(optimization_data_process_default_settings)
                    process = SurrogateModellingComponentFactory(self.models, process_settings, self.optimization_problem)
                    self.optimization_problem.AddProcess(process_type, process)

    def _CreateAlgorithm(self):
        default_settings = Kratos.Parameters("""{
            "module" : "KratosMultiphysics.OptimizationApplication.algorithms"
        }""")
        algorithm_settings = self.project_parameters["algorithm_settings"]
        algorithm_settings.AddMissingParameters(default_settings)
        self.__algorithm = SurrogateModellingComponentFactory(self.models, algorithm_settings, self.optimization_problem)
