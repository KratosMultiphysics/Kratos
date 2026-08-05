import pickle
from pathlib import Path

import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.buffered_dict import BufferedDict
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView

def Factory(_: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> Kratos.Process:
    if not parameters.Has("settings"):
        raise RuntimeError(f"OptimizationProblemRestartInputProcess instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")
    return OptimizationProblemRestartInputProcess(parameters["settings"], optimization_problem)

def RestoreBufferedDict(node: BufferedDict, snapshot: dict, shape_owner, echo_level: int) -> None:
    """Replays a snapshot produced by optimization_problem_restart_output_process back into a
    live (freshly constructed, still empty) BufferedDict node.

    Tensor-valued leaves need a correctly-shaped placeholder Kratos.TensorAdaptors.* object to
    write the restored numpy array into (they cannot be unpickled directly). "shape_owner" is
    the component (Control/MasterControl) this subtree's tensors are shaped like -- see
    RestoreOptimizationProblemData for how it is resolved. If it is None, tensor-valued leaves
    are skipped (only scalars restored); this only happens for subtrees that get unconditionally
    recomputed before being read again, e.g. response value/gradient caches.
    """
    for step_index, slot in enumerate(snapshot["slots"]):
        for key, leaf in slot.items():
            kind = leaf["kind"]
            if kind == "scalar":
                node.SetValue(key, leaf["value"], step_index, overwrite=True)
                continue

            if shape_owner is None:
                if echo_level > 0:
                    Kratos.Logger.PrintWarning("OptimizationProblemRestartInputProcess", f"Skipping restore of tensor-valued key \"{key}\" (no known shape owner for this subtree).")
                continue

            placeholder = shape_owner.GetEmptyField()
            if kind == "tensor":
                placeholder.data[:] = leaf["array"]
            elif kind == "combined_tensor":
                parts = leaf["parts"]
                tensor_adaptors = placeholder.GetTensorAdaptors()
                if len(parts) != len(tensor_adaptors):
                    raise RuntimeError(
                        f"Restart data for key \"{key}\" has {len(parts)} parts, but the current placeholder has {len(tensor_adaptors)} parts."
                    )
                for part_array, part_ta in zip(parts, tensor_adaptors):
                    part_ta.data[:] = part_array
                # a combined tensor adaptor's own .data is a separate buffer from its parts';
                # CollectData() resyncs it from the parts we just wrote.
                placeholder.CollectData()
            else:
                raise RuntimeError(f"Unknown restart data kind \"{kind}\" for key \"{key}\".")

            if step_index == 0:
                # shape_owner.Update() writes straight to the shared mesh, so only the current
                # step's data (step_index 0) may be pushed there. Historical steps are restored
                # into the BufferedDict below for bookkeeping only (e.g. GetRelativeChange()'s
                # previous-step value), never onto the mesh.
                shape_owner.Update(placeholder)
            node.SetValue(key, placeholder, step_index, overwrite=True)

    for name, sub_snapshot in snapshot["sub_items"].items():
        sub_items = node.GetSubItems()
        if name not in sub_items:
            if echo_level > 0:
                Kratos.Logger.PrintWarning("OptimizationProblemRestartInputProcess", f"Skipping restore of \"{name}\" (not present yet in the live optimization problem).")
            continue
        RestoreBufferedDict(sub_items[name], sub_snapshot, shape_owner, echo_level)

def RestoreOptimizationProblemData(optimization_problem: OptimizationProblem, data_snapshot: dict, echo_level: int) -> None:
    root = optimization_problem.GetProblemDataContainer()

    # restore root-level scalar leaves too (currently just "step"; ExecuteInitialize also sets
    # it explicitly via SetStep, so this is redundant but harmless).
    for step_index, slot in enumerate(data_snapshot["slots"]):
        for key, leaf in slot.items():
            if leaf["kind"] == "scalar":
                root.SetValue(key, leaf["value"], step_index, overwrite=True)

    type_name_to_class = {cls.__name__: cls for cls in optimization_problem.GetComponentContainer().keys()}

    master_controls = list(optimization_problem.GetListOfMasterControls())
    algorithm_shape_owner = master_controls[0] if len(master_controls) > 0 else None

    root_sub_items = root.GetSubItems()
    for type_name, type_snapshot in data_snapshot["sub_items"].items():
        if type_name not in root_sub_items:
            continue
        type_node = root_sub_items[type_name]
        type_node_sub_items = type_node.GetSubItems()

        for component_name, component_snapshot in type_snapshot["sub_items"].items():
            if component_name not in type_node_sub_items:
                if echo_level > 0:
                    Kratos.Logger.PrintWarning("OptimizationProblemRestartInputProcess", f"Skipping restore of \"{type_name}/{component_name}\" (component not present in this run).")
                continue

            shape_owner = None
            if type_name == "object":
                # "object" holds ComponentDataView(<string>, ...) subtrees, e.g. the "algorithm"
                # buffered data every Algorithm implementation uses -- that one is mesh-shaped
                # like the master control. Everything else under "object" (response value/gradient
                # caches, projection helper bookkeeping) gets recomputed before being read again,
                # so it's fine to leave those tensors unresolved.
                if component_name == "algorithm":
                    shape_owner = algorithm_shape_owner
            elif type_name in type_name_to_class:
                component = optimization_problem.GetComponent(component_name, type_name_to_class[type_name])
                if hasattr(component, "GetEmptyField"):
                    shape_owner = component
                elif hasattr(component, "GetMasterControl"):
                    shape_owner = component.GetMasterControl()

            RestoreBufferedDict(type_node_sub_items[component_name], component_snapshot, shape_owner, echo_level)

class OptimizationProblemRestartInputProcess(Kratos.Process):
    def GetDefaultParameters(self) -> Kratos.Parameters:
        return Kratos.Parameters("""{
            "restart_files_path": "Optimization_Restart",
            "restart_file_name" : "restart_<step>.pkl",
            "restart_load_step" : "latest",
            "echo_level"         : 0
        }""")

    def __init__(self, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        Kratos.Process.__init__(self)

        # "restart_load_step" is either the string "latest" or an explicit step (int).
        # ValidateAndAssignDefaults requires matching types, so the default's type is chosen
        # to match whatever the user provided (mirrors the "scaled_ref_value" pattern used in
        # standardized_rgp_constraint.py).
        default_parameters = self.GetDefaultParameters()
        if parameters.Has("restart_load_step") and parameters["restart_load_step"].IsInt():
            default_parameters["restart_load_step"].SetInt(0)
        parameters.ValidateAndAssignDefaults(default_parameters)

        self.optimization_problem = optimization_problem
        self.restart_files_path = Path(parameters["restart_files_path"].GetString())
        self.restart_file_name = parameters["restart_file_name"].GetString()
        self.echo_level = parameters["echo_level"].GetInt()

        if "<step>" not in self.restart_file_name:
            raise RuntimeError(f"\"restart_file_name\" should contain the \"<step>\" placeholder [ restart_file_name = \"{self.restart_file_name}\" ].")

        load_step_param = parameters["restart_load_step"]
        self.restart_load_step = None if load_step_param.IsString() and load_step_param.GetString() == "latest" else load_step_param.GetInt()

    def ExecuteInitialize(self) -> None:
        file_path = self.__ResolveCheckpointFile()
        if file_path is None:
            Kratos.Logger.PrintInfo(self.__class__.__name__, f"No restart checkpoint found under \"{self.restart_files_path}\". Starting a fresh run.")
            return

        Kratos.Logger.PrintInfo(self.__class__.__name__, f"Restoring optimization problem state from \"{file_path}\".")
        with open(file_path, "rb") as file_input:
            payload = pickle.load(file_input)

        # Controls are normally initialized later, from the algorithm's own Initialize(). Restoring
        # tensor data needs GetEmptyField()/Update() on already-initialized controls, so initialize
        # the master control(s) here first; Control.Initialize() is idempotent, so the algorithm
        # initializing them again afterwards is harmless.
        for master_control in self.optimization_problem.GetListOfMasterControls():
            master_control.Initialize()

        self.optimization_problem.SetStep(payload["step"])
        RestoreOptimizationProblemData(self.optimization_problem, payload["data"], self.echo_level)

        # The checkpoint is written by Algorithm.Output(), which always runs before that step's own
        # UpdateControl() call. So the design just restored is the one that *produced* the
        # checkpointed step's results, not the design the live run had moved on to by the time it
        # stopped. Apply the checkpointed step's pending "control_field_update" now to get there.
        master_controls = list(self.optimization_problem.GetListOfMasterControls())
        if master_controls:
            algorithm_data = ComponentDataView("algorithm", self.optimization_problem).GetBufferedData()
            if algorithm_data.HasValue("control_field_update"):
                master_control = master_controls[0]
                resumed_field = master_control.GetControlField()
                resumed_field.data[:] += algorithm_data.GetValue("control_field_update").data
                # a combined tensor adaptor's .data is disconnected from its parts; Update() reads
                # the parts, so the change has to be pushed down via StoreData() first.
                Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor(resumed_field, perform_store_data_recursively=False, copy=False).StoreData()
                master_control.Update(resumed_field)

        # The checkpoint was taken before the checkpointed step's own AdvanceStep() call (which
        # only runs when the run doesn't converge there), so both restored buffer slots are still
        # "occupied" (current = checkpointed step, previous = the step before it). Advance now to
        # free up a fresh current slot for the next iteration, as AdvanceStep() would have done had
        # the checkpointed run continued.
        self.optimization_problem.AdvanceStep()

        # push the restored design through the normal Control.Update path so control-specific side
        # effects (filtering caches, etc.) run on it too, not just the raw writes above.
        for master_control in self.optimization_problem.GetListOfMasterControls():
            master_control.Update(master_control.GetControlField())

        Kratos.Logger.PrintInfo(self.__class__.__name__, f"Restored optimization problem at step {payload['step']}.")

    def __ResolveCheckpointFile(self) -> 'Path | None':
        if self.restart_load_step is not None:
            file_path = (self.restart_files_path / self.restart_file_name.replace("<step>", str(self.restart_load_step))).resolve()
            if self.restart_files_path.resolve() not in file_path.parents:
                raise RuntimeError(f"Resolved restart checkpoint path is outside restart_files_path. [ file_path = \"{file_path}\" ].")
            if not file_path.is_file():
                raise FileNotFoundError(f"Restart checkpoint file not found: \"{file_path}\".")
            return file_path

        checkpoint_glob = self.restart_file_name.replace("<step>", "*")
        checkpoints = list(self.restart_files_path.glob(checkpoint_glob)) if self.restart_files_path.is_dir() else []
        if len(checkpoints) == 0:
            return None

        return max(checkpoints, key=lambda p: p.stat().st_mtime)
