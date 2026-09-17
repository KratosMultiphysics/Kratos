import pickle

import KratosMultiphysics as Kratos
from KratosMultiphysics.process_factory import KratosProcessFactory
from KratosMultiphysics.OptimizationApplication.algorithms.algorithm import Algorithm
from KratosMultiphysics.OptimizationApplication.model_part_controllers.model_part_controller import ModelPartController
from KratosMultiphysics.OptimizationApplication.execution_policies.execution_policy_decorator import ExecutionPolicyDecorator
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem_utilities import OptimizationComponentFactory
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import CallOnAll
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import OptimizationAnalysisTimeLogger
from KratosMultiphysics.OptimizationApplication.utilities.restart_file_naming import ResolveRestartLoadStep, SplitRestartFileName
from KratosMultiphysics.OptimizationApplication.processes.optimization_problem_restart_output_process import OptimizationProblemRestartOutputProcess
from KratosMultiphysics.OptimizationApplication.processes.optimization_problem_restart_input_process import OptimizationProblemRestartInputProcess

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
            },
            "restart_settings": {
                "save_restart"           : false,
                "load_restart"           : false,
                "restart_file_name"      : "Optimization_Restart/restart_<step>.pkl",
                "restart_save_frequency" : 1,
                "restart_load_step"      : "latest",
                "max_files_to_keep"      : -1,
                "echo_level"             : 0,
                "model_parts_settings"   : {
                    "save_model_parts"  : true,
                    "load_model_parts"  : true,
                    "serializer_trace"  : "no_trace",
                    "clean_before_save" : true
                }
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
        # _CreateRestart must run after _CreateProcesses: _CreateProcesses() unconditionally
        # calls optimization_problem.AddProcessType(process_type) for every process type, which
        # resets (rather than appends to) that type's process list.
        self._CreateRestart()

    def Initialize(self):
        CallOnAll(self.__list_of_model_part_controllers, ModelPartController.ImportModelPart)
        CallOnAll(self.__list_of_model_part_controllers, ModelPartController.Initialize)
        CallOnAll(self.optimization_problem.GetListOfExecutionPolicies(), ExecutionPolicyDecorator.Initialize)
        for process_type in self.__algorithm.GetProcessesOrder():
            CallOnAll(self.optimization_problem.GetListOfProcesses(process_type), Kratos.Process.ExecuteInitialize)

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
        for model_part_controller_settings in self.project_parameters["model_parts"].values():
            model_part_controller_settings.AddMissingParameters(default_settings)
            model_part_controller: ModelPartController = OptimizationComponentFactory(self.model, model_part_controller_settings, self.optimization_problem)
            self.__list_of_model_part_controllers.append(model_part_controller)

    def _CreateAnalyses(self):
        default_settings = Kratos.Parameters("""{
            "module": "KratosMultiphysics.OptimizationApplication.execution_policies"
        }""")
        for analyses_settings in self.project_parameters["analyses"].values():
            analyses_settings.AddMissingParameters(default_settings)
            execution_policy = OptimizationComponentFactory(self.model, analyses_settings, self.optimization_problem)
            self.optimization_problem.AddComponent(execution_policy)

    def _CreateResponses(self):
        default_settings = Kratos.Parameters("""{
            "module": "KratosMultiphysics.OptimizationApplication.responses"
        }""")
        for response_settings in self.project_parameters["responses"].values():
            response_settings.AddMissingParameters(default_settings)
            response_function: ResponseFunction = OptimizationComponentFactory(self.model, response_settings, self.optimization_problem)
            self.optimization_problem.AddComponent(response_function)

    def _CreateControls(self):
        default_settings = Kratos.Parameters("""{
            "module" : "KratosMultiphysics.OptimizationApplication.controls"
        }""")
        for control_settings in self.project_parameters["controls"].values():
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
                for process_settings in optimization_data_processes[process_type].values():
                    process_settings.AddMissingParameters(optimization_data_process_default_settings)
                    process = OptimizationComponentFactory(self.model, process_settings, self.optimization_problem)
                    self.optimization_problem.AddProcess(process_type, process)

    def _CreateRestart(self):
        restart_settings = self.project_parameters["restart_settings"]

        default_restart_settings = self.GetDefaultParameters()["restart_settings"]
        # "restart_load_step" is either the string "latest" or an explicit step (int).
        # ValidateAndAssignDefaults requires matching types, so the default's type is chosen to
        # match whatever the user provided (mirrors the same workaround already used by
        # OptimizationProblemRestartInputProcess.__init__).
        if restart_settings.Has("restart_load_step") and restart_settings["restart_load_step"].IsInt():
            default_restart_settings["restart_load_step"].SetInt(0)
        restart_settings.ValidateAndAssignDefaults(default_restart_settings)
        # ValidateAndAssignDefaults is not recursive (see _CreateProcesses' identical pattern for
        # "processes"), so a partially-specified "model_parts_settings" needs its own explicit
        # default-filling pass here.
        restart_settings["model_parts_settings"].ValidateAndAssignDefaults(default_restart_settings["model_parts_settings"])

        save_restart = restart_settings["save_restart"].GetBool()
        load_restart = restart_settings["load_restart"].GetBool()
        if not save_restart and not load_restart:
            return

        restart_files_path, restart_file_name = SplitRestartFileName(restart_settings["restart_file_name"].GetString())
        echo_level = restart_settings["echo_level"].GetInt()

        model_parts_settings = restart_settings["model_parts_settings"]

        restart_capable_controllers = [controller for controller in self.__list_of_model_part_controllers if controller.SupportsRestart()]

        if load_restart:
            payload = None
            resolved_step = ResolveRestartLoadStep(restart_files_path, restart_file_name, restart_settings["restart_load_step"])
            if resolved_step is not None:
                file_path = (restart_files_path / restart_file_name.replace("<step>", str(resolved_step))).resolve()
                if restart_files_path.resolve() not in file_path.parents:
                    raise RuntimeError(f"Resolved restart checkpoint path is outside restart_files_path. [ file_path = \"{file_path}\" ].")
                if not file_path.is_file():
                    raise FileNotFoundError(f"Restart checkpoint file not found: \"{file_path}\".")

                Kratos.Logger.PrintInfo("OptimizationAnalysis", f"Loading restart checkpoint from \"{file_path}\".")
                with open(file_path, "rb") as file_input:
                    payload = pickle.load(file_input)

                # Load restart-capable model parts directly here, synchronously, before
                # Initialize() ever calls ModelPartController.ImportModelPart() -- confirmed safe
                # even though _CreateModelPartControllers() already created these (empty) model
                # parts above: Kratos.ModelPart.load() correctly populates an already-created,
                # empty ModelPart of the same name in place. Reusing payload["serializer"] (rather
                # than each controller reading its own file) is required for the BufferedDict's
                # TensorAdaptor leaves (restored later, from this exact same Serializer instance,
                # by OptimizationProblemRestartInputProcess) to stay pointer-linked to these model
                # parts' Node/Element/Condition containers.
                if model_parts_settings["load_model_parts"].GetBool():
                    for controller in restart_capable_controllers:
                        model_part = controller.GetModelPart()
                        if model_part.Name in payload["model_part_names"]:
                            payload["serializer"].Load(model_part.Name, model_part)
                        elif echo_level > 0:
                            Kratos.Logger.PrintWarning("OptimizationAnalysis", f"No restart data found for model part \"{model_part.Name}\"; it will be imported normally.")

            input_process = OptimizationProblemRestartInputProcess(payload, self.optimization_problem, echo_level)
            self.optimization_problem.AddProcess("auxiliary_processes", input_process)

        if save_restart:
            output_settings = Kratos.Parameters("{}")
            for key in ("restart_file_name", "restart_save_frequency", "max_files_to_keep", "echo_level"):
                output_settings.AddValue(key, restart_settings[key])
            output_model_parts_settings = Kratos.Parameters("{}")
            for key in ("save_model_parts", "serializer_trace", "clean_before_save"):
                output_model_parts_settings.AddValue(key, model_parts_settings[key])
            output_settings.AddValue("model_parts_settings", output_model_parts_settings)

            list_of_model_parts = [controller.GetModelPart() for controller in restart_capable_controllers]
            output_process = OptimizationProblemRestartOutputProcess(output_settings, self.optimization_problem, list_of_model_parts)
            self.optimization_problem.AddProcess("output_processes", output_process)

    def _CreateAlgorithm(self):
        default_settings = Kratos.Parameters("""{
            "module" : "KratosMultiphysics.OptimizationApplication.algorithms"
        }""")
        algorithm_settings = self.project_parameters["algorithm_settings"]
        algorithm_settings.AddMissingParameters(default_settings)
        self.__algorithm = OptimizationComponentFactory(self.model, algorithm_settings, self.optimization_problem)

    def GetAlgorithm(self):
        return self.__algorithm
