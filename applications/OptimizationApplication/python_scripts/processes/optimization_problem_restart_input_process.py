import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.buffered_dict import BufferedDict
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView

def Factory(*_args, **_kwargs):
    raise RuntimeError(
        "\"optimization_problem_restart_input_process\" can no longer be configured under "
        "\"processes\"; restart is now configured via the top-level \"restart_settings\" section "
        "of the optimization parameters (see OptimizationAnalysis.GetDefaultParameters())."
    )

def _GetExistingParentSubItem(root: BufferedDict, bufferdict_key: str) -> 'BufferedDict | None':
    """Walks bufferdict_key's parent path (all but its last "/"-separated segment) through
    already-existing sub items only, never creating new ones.

    Returns the parent BufferedDict node if the whole chain already exists live, or None if any
    segment is missing -- meaning the checkpoint refers to a component/subtree that isn't part of
    the current run (e.g. a Control that was removed from the optimization parameters since the
    checkpoint was written). Callers skip such leaves rather than letting BufferedDict.SetValue's
    own "/"-path auto-creation silently insert orphaned structure into the live tree.
    """
    node = root
    segments = bufferdict_key.split("/")
    for segment in segments[:-1]:
        sub_items = node.GetSubItems()
        if segment not in sub_items:
            return None
        node = sub_items[segment]
    return node

def RestoreScalars(root: BufferedDict, scalar_snapshot: dict, echo_level: int) -> None:
    for (bufferdict_key, step_index), value in scalar_snapshot.items():
        if _GetExistingParentSubItem(root, bufferdict_key) is None:
            if echo_level > 0:
                Kratos.Logger.PrintWarning("OptimizationProblemRestartInputProcess", f"Skipping restore of \"{bufferdict_key}\" (not present yet in the live optimization problem).")
            continue
        root.SetValue(bufferdict_key, value, step_index, overwrite=True)

def RestoreTensors(root: BufferedDict, serializer: Kratos.StreamSerializer, tensor_kinds: dict, echo_level: int) -> None:
    """Restores every TensorAdaptor leaf directly via the shared Serializer.

    Since `serializer` also holds (or, for BufferedDict tensors backed by a container outside any
    restart-capable model part, simply doesn't reference) whichever ModelPart(s) were saved
    alongside it, a natively-`Load()`ed tensor already comes back with the right shape *and* the
    right container reference -- no shape_owner/GetEmptyField() resolution needed here, unlike the
    numpy-array-based restore this replaces.
    """
    for (bufferdict_key, step_index), (serializer_tag, kind) in tensor_kinds.items():
        if _GetExistingParentSubItem(root, bufferdict_key) is None:
            if echo_level > 0:
                Kratos.Logger.PrintWarning("OptimizationProblemRestartInputProcess", f"Skipping restore of \"{bufferdict_key}\" (not present yet in the live optimization problem).")
            continue
        placeholder = getattr(Kratos.TensorAdaptors, kind)()
        serializer.Load(serializer_tag, placeholder)
        root.SetValue(bufferdict_key, placeholder, step_index, overwrite=True)

class OptimizationProblemRestartInputProcess(Kratos.Process):
    def __init__(self, payload: 'dict | None', optimization_problem: OptimizationProblem, echo_level: int = 0):
        """`payload` is the dict OptimizationAnalysis._CreateRestart() unpickled from the restart
        checkpoint file (None if no checkpoint was found, e.g. a first-ever run with load_restart
        enabled) -- see _CreateRestart() for why it's threaded in already-unpickled: it also has
        already used `payload["serializer"]` to load the restart-capable ModelPart(s) into
        `optimization_problem`'s model, and this process's ExecuteInitialize() must keep restoring
        BufferedDict tensors from that *exact same* Serializer instance, not a freshly re-read one,
        for the ModelPart<->tensor pointer relinking to hold.
        """
        Kratos.Process.__init__(self)
        self.payload = payload
        self.optimization_problem = optimization_problem
        self.echo_level = echo_level

    def ExecuteInitialize(self) -> None:
        if self.payload is None:
            Kratos.Logger.PrintInfo(self.__class__.__name__, "No restart checkpoint found. Starting a fresh run.")
            return

        Kratos.Logger.PrintInfo(self.__class__.__name__, f"Restoring optimization problem state from step {self.payload['step']}.")

        # Controls are normally initialized later, from the algorithm's own Initialize(). The
        # control_field_update reapplication below needs GetControlField()/Update() on
        # already-initialized controls, so initialize the master control(s) here first;
        # Control.Initialize() is idempotent, so the algorithm initializing them again afterwards
        # is harmless.
        for master_control in self.optimization_problem.GetListOfMasterControls():
            master_control.Initialize()

        self.optimization_problem.SetStep(self.payload["step"])
        root = self.optimization_problem.GetProblemDataContainer()
        RestoreScalars(root, self.payload["scalars"], self.echo_level)
        RestoreTensors(root, self.payload["serializer"], self.payload["tensor_kinds"], self.echo_level)

        # MasterControl.GetControlField()/Control.GetControlField() read live, current data
        # straight off each control's own storage (ultimately the mesh) -- they do NOT consult the
        # BufferedDict's "control_field" entry we just restored above. So restoring that entry into
        # the BufferedDict alone does not, by itself, make the design "reappear" on the mesh; it
        # has to be explicitly pushed through Control.Update().
        #
        # Route this through GetEmptyField() rather than pushing the restored "control_field"
        # tensor directly: when the model part(s) it's built over are *not themselves* also
        # restart-loaded (restart_settings.model_parts_settings.load_model_parts=false, or a
        # control's model part simply isn't restart-capable -- e.g. a
        # connectivity_preserving_model_part_controller-derived one), that tensor's own
        # Serializer-restored container is a fresh, disconnected reconstruction, not the live model
        # part Control.Update() expects (confirmed by test_system_identification_restart.py: this
        # raised "Updates for the required element container not found" before this fix).
        # GetEmptyField() is always correctly linked to the live mesh; copying just the restored
        # *values* into it sidesteps the container question entirely, matching how the
        # "control_field_update" reapplication below already safely only ever reads restored data
        # via .data (never relies on a restored tensor's own container).
        master_controls = list(self.optimization_problem.GetListOfMasterControls())
        if master_controls:
            algorithm_data = ComponentDataView("algorithm", self.optimization_problem).GetBufferedData()
            if algorithm_data.HasValue("control_field"):
                master_control = master_controls[0]
                live_field = master_control.GetEmptyField()
                live_field.data[:] = algorithm_data.GetValue("control_field").data
                # a combined tensor adaptor's .data is disconnected from its parts; Update() reads
                # the parts, so the change has to be pushed down via StoreData() first (same as the
                # control_field_update reapplication below).
                Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor(live_field, perform_store_data_recursively=False, copy=False).StoreData()
                master_control.Update(live_field)

            # The checkpoint is written by Algorithm.Output(), which always runs before that step's
            # own UpdateControl() call. So the design just restored is the one that *produced* the
            # checkpointed step's results, not the design the live run had moved on to by the time
            # it stopped. Apply the checkpointed step's pending "control_field_update" now to get
            # there.
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

        Kratos.Logger.PrintInfo(self.__class__.__name__, f"Restored optimization problem at step {self.payload['step']}.")
