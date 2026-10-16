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
from KratosMultiphysics.OptimizationApplication.utilities.buffered_dict import BufferedDict
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import OptimizationAnalysisTimeLogger
from KratosMultiphysics.OptimizationApplication.utilities.restart_file_naming import ResolveRestartLoadStep, SplitRestartFileName
from KratosMultiphysics.OptimizationApplication.utilities.optimization_restart import OptimizationProblemRestartOutput, CreateTensorAdaptorLoadTarget, CheckSerialRestart, ImportApplications

def _RestoreBufferedDictNodes(root: BufferedDict, node_metadata: dict) -> None:
    """Recreates every checkpointed BufferedDict node missing from root, parents first, with its saved
    (buffer_size, clear_buffer_when_advancing). Nodes that already exist (sized by their component)
    are kept, and so are nodes that hold no values (e.g. an empty history buffer)."""
    for path in sorted(node_metadata.keys(), key=lambda node_path: node_path.count("/")):
        parent_path, _, name = path.rpartition("/")
        parent = root.GetValue(parent_path) if parent_path else root
        if name not in parent.GetSubItems():
            buffer_size, clear_buffer_when_advancing = node_metadata[path]
            parent.SetValue(name, BufferedDict(buffer_size, clear_buffer_when_advancing))

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
                "serializer_trace"       : "no_trace",
                "echo_level"             : 0
            }
        }""")

    def __init__(self, model: Kratos.Model, project_parameters: Kratos.Parameters):
        self.model = model
        self.project_parameters = project_parameters
        self.project_parameters.ValidateAndAssignDefaults(self.GetDefaultParameters())

        self.optimization_problem = OptimizationProblem(self.project_parameters["problem_data"]["echo_level"].GetInt())

        self.__list_of_model_part_controllers: 'list[ModelPartController]' = []
        self.__algorithm: Algorithm = None

        # before any model part is created: Model.load() fills the still-empty Model in place
        self.__restart_payload = self._LoadRestartModel()

        self._CreateModelPartControllers()
        self._CreateAnalyses()
        self._CreateControls()
        self._CreateResponses()
        self._CreateAlgorithm()
        self._CheckRestartSupport()

    def Initialize(self):
        CallOnAll(self.__list_of_model_part_controllers, ModelPartController.ImportModelPart)
        CallOnAll(self.__list_of_model_part_controllers, ModelPartController.Initialize)
        self._CreateProcesses()
        # after all _Create*() calls (their SetDataBuffer() sizes the BufferedDict nodes) and before
        # any process or the algorithm is initialized, so that they all see the restored state.
        self._CreateRestart()
        for process_type in self.__algorithm.GetProcessesOrder():
            CallOnAll(self.optimization_problem.GetListOfProcesses(process_type), Kratos.Process.ExecuteInitialize)
            CallOnAll(self.optimization_problem.GetListOfProcesses(process_type), Kratos.Process.ExecuteBeforeSolutionLoop)
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

    def _GetValidatedRestartSettings(self) -> Kratos.Parameters:
        restart_settings = self.project_parameters["restart_settings"]

        default_restart_settings = self.GetDefaultParameters()["restart_settings"]
        # "restart_load_step" is "latest" or an int; match the default's type for the validation.
        if restart_settings.Has("restart_load_step") and restart_settings["restart_load_step"].IsInt():
            default_restart_settings["restart_load_step"].SetInt(0)
        restart_settings.ValidateAndAssignDefaults(default_restart_settings)
        return restart_settings

    def _LoadRestartModel(self) -> 'dict | None':
        """Loads the checkpoint's Model into self.model if load_restart is on.

        Must run before any model part is created, since Model.load() fills the empty Model in place.
        Imports the applications recorded in the checkpoint, reapplies the saved serializer flags and
        marks the problem as restarted.

        Returns:
            dict | None: the unpickled payload for _CreateRestart(), or None if load_restart is off
            or "restart_load_step" is "latest" and no checkpoint exists yet.
        """
        restart_settings = self._GetValidatedRestartSettings()
        if not restart_settings["load_restart"].GetBool():
            return None
        CheckSerialRestart()

        restart_files_path, restart_file_name = SplitRestartFileName(restart_settings["restart_file_name"].GetString())
        resolved_step = ResolveRestartLoadStep(restart_files_path, restart_file_name, restart_settings["restart_load_step"])
        if resolved_step is None:
            Kratos.Logger.PrintInfo("OptimizationAnalysis", f"\"load_restart\" is on but no checkpoint matching \"{restart_files_path / restart_file_name}\" was found; starting from scratch.")
            return None

        file_path = (restart_files_path / restart_file_name.replace("<step>", str(resolved_step))).resolve()
        if restart_files_path.resolve() not in file_path.parents:
            raise RuntimeError(f"Resolved restart checkpoint path is outside restart_files_path. [ file_path = \"{file_path}\" ].")
        if not file_path.is_file():
            raise FileNotFoundError(f"Restart checkpoint file not found: \"{file_path}\".")

        Kratos.Logger.PrintInfo("OptimizationAnalysis", f"Loading restart checkpoint from \"{file_path}\".")
        with open(file_path, "rb") as file_input:
            payload: dict = pickle.load(file_input)

        # a fresh interpreter has not imported the applications whose entities the Model holds yet
        ImportApplications(payload["applications"])

        serializer: Kratos.StreamSerializer = payload["serializer"]
        flags = Kratos.Flags()
        serializer.Load("SerializerFlags", flags)
        serializer.Set(flags)
        serializer.Load("Model", self.model)

        self.optimization_problem.SetRestarted(True)
        return payload

    def _CheckRestartSupport(self) -> None:
        """Raises if save_restart or load_restart is on for an algorithm without SUPPORTS_RESTART."""
        restart_settings = self._GetValidatedRestartSettings()
        if (restart_settings["save_restart"].GetBool() or restart_settings["load_restart"].GetBool()) and not self.__algorithm.SUPPORTS_RESTART:
            raise RuntimeError(
                f"\"{type(self.__algorithm).__name__}\" does not call SaveRestartFile() from its "
                "Solve() loop, so restart_settings.save_restart/load_restart cannot be honoured "
                "for it.")

    def _CreateRestart(self):
        """Restores the checkpointed BufferedDict and step, and installs the restart writer.

        Every saved value is restored unconditionally; missing nodes are recreated from the saved
        node metadata. Must run after all components are created (they size their data buffers) and
        before any of them is initialized.
        """
        restart_settings = self._GetValidatedRestartSettings()
        save_restart = restart_settings["save_restart"].GetBool()

        if self.__restart_payload is not None:
            payload = self.__restart_payload
            root = self.optimization_problem.GetProblemDataContainer()
            node_metadata = payload["node_metadata"]
            serializer: Kratos.StreamSerializer = payload["serializer"]

            # nodes of components absent from this run are recreated from node_metadata, then every
            # value is restored verbatim.
            _RestoreBufferedDictNodes(root, node_metadata)
            for (bufferdict_key, step_index), value in payload["scalars"].items():
                root.SetValue(bufferdict_key, value, step_index, overwrite=True)

            for (bufferdict_key, step_index), (serializer_tag, type_name) in payload["tensor_kinds"].items():
                tensor_adaptor = CreateTensorAdaptorLoadTarget(type_name)
                serializer.Load(serializer_tag, tensor_adaptor)
                root.SetValue(bufferdict_key, tensor_adaptor, step_index, overwrite=True)

            # the checkpoint is written after the step's UpdateControl(), so advance once, as the
            # checkpointed run would have done before its next iteration.
            self.optimization_problem.SetStep(payload["step"])
            self.optimization_problem.AdvanceStep()

            Kratos.Logger.PrintInfo("OptimizationAnalysis", f"Restored optimization problem state from step {payload['step']}.")

        if save_restart:
            output_settings = Kratos.Parameters("{}")
            for key in ("restart_file_name", "restart_save_frequency", "max_files_to_keep", "serializer_trace", "echo_level"):
                output_settings.AddValue(key, restart_settings[key])
            restart_output = OptimizationProblemRestartOutput(output_settings, self.optimization_problem, self.model)
            self.optimization_problem.SetRestartOutput(restart_output)

    def _CreateAlgorithm(self):
        default_settings = Kratos.Parameters("""{
            "module" : "KratosMultiphysics.OptimizationApplication.algorithms"
        }""")
        algorithm_settings = self.project_parameters["algorithm_settings"]
        algorithm_settings.AddMissingParameters(default_settings)
        self.__algorithm = OptimizationComponentFactory(self.model, algorithm_settings, self.optimization_problem)

    def GetAlgorithm(self):
        return self.__algorithm
