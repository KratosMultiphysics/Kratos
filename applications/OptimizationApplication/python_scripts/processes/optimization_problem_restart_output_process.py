import os
import pickle

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.buffered_dict import BufferedDict
from KratosMultiphysics.OptimizationApplication.utilities.restart_file_naming import SplitRestartFileName, ParseStepFromFileName

def Factory(*_args, **_kwargs):
    raise RuntimeError(
        "\"optimization_problem_restart_output_process\" can no longer be configured under "
        "\"processes\"; restart is now configured via the top-level \"restart_settings\" section "
        "of the optimization parameters (see OptimizationAnalysis.GetDefaultParameters())."
    )

_SERIALIZER_TRACE_TYPES = {
    "no_trace":    Kratos.SerializerTraceType.SERIALIZER_NO_TRACE,
    "trace_error": Kratos.SerializerTraceType.SERIALIZER_TRACE_ERROR,
    "trace_all":   Kratos.SerializerTraceType.SERIALIZER_TRACE_ALL,
}

# (base type, matching combined type, matching NDData type) per dtype.
_DTYPE_TENSOR_ADAPTOR_TYPES = (
    (Kratos.TensorAdaptors.DoubleTensorAdaptor, Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor, Kratos.DoubleNDData),
    (Kratos.TensorAdaptors.IntTensorAdaptor, Kratos.TensorAdaptors.IntCombinedTensorAdaptor, Kratos.IntNDData),
    (Kratos.TensorAdaptors.BoolTensorAdaptor, Kratos.TensorAdaptors.BoolCombinedTensorAdaptor, Kratos.BoolNDData),
)

# isinstance against just the 3 base types also matches the corresponding combined subtype, since
# Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor (etc.) is bound as a Python subclass of
# Kratos.TensorAdaptors.DoubleTensorAdaptor (etc.) -- see kratos/python/add_tensor_adaptors_to_python.cpp.
_TENSOR_ADAPTOR_BASE_TYPES = tuple(base for base, _, _ in _DTYPE_TENSOR_ADAPTOR_TYPES)
_COMBINED_TENSOR_ADAPTOR_TYPES = tuple(combined for _, combined, _ in _DTYPE_TENSOR_ADAPTOR_TYPES)

# the only leaf (non-combined) types registered with Kratos::Serializer (kratos_application.cpp)
# and therefore safe to save/load directly, as themselves.
_EXACT_NATIVE_TENSOR_ADAPTOR_TYPES = _TENSOR_ADAPTOR_BASE_TYPES

_SCALAR_LEAF_TYPES = (int, float, bool, str, list, dict, type(None))

def _ToNativelySerializable(value):
    """Returns a TensorAdaptor safe to save/load directly through Kratos.Serializer.

    A combined adaptor's own C++ save() serializes each child through a polymorphic pointer, which
    throws for a further-derived, unregistered child type (e.g. SystemIdentificationApplication's
    PropertiesVariableTensorAdaptor) -- so recurse into every combined adaptor's children first and
    rebuild a plain combined instance from their laundered results, regardless of the combined
    adaptor's own exact type. The rebuilt combined instance's own top-level .data is freshly
    allocated and disconnected from its (already-correct) children until CollectData() pulls
    values in from them -- save() serializes that top-level buffer too, so skipping this would
    silently save garbage. perform_collect_data_recursively=False keeps that pull to "children's
    .data -> this .data" only: recursing further ("child's container -> child's .data") would be
    wrong for a child that doesn't represent something re-readable from its container at all (e.g.
    a momentum/moment estimate, whose container is only there for shape/identity).

    A leaf whose exact type is already registered returns unchanged; otherwise it is "laundered"
    into a freshly-constructed plain base-type instance wrapping the *same* container and a copy of
    the *same* data -- dropping whatever extra, subclass-specific state the derived type carries,
    since restart only needs to round-trip the current values and which container they belong to.

    Returns None if a leaf has no container to launder through (HasContainer() is False) and isn't
    itself one of the registered types -- the caller falls back to skipping it, with a warning, same
    as any other genuinely unsupported leaf.
    """
    if isinstance(value, _COMBINED_TENSOR_ADAPTOR_TYPES):
        children = [_ToNativelySerializable(child) for child in value.GetTensorAdaptors()]
        if any(child is None for child in children):
            return None
        for _, combined_type, _ in _DTYPE_TENSOR_ADAPTOR_TYPES:
            if isinstance(value, combined_type):
                native_value = combined_type(children, perform_collect_data_recursively=False, copy=False)
                native_value.CollectData()
                return native_value
        return None
    if type(value) in _EXACT_NATIVE_TENSOR_ADAPTOR_TYPES:
        return value
    if not value.HasContainer():
        return None
    for base_type, _, nd_data_type in _DTYPE_TENSOR_ADAPTOR_TYPES:
        if isinstance(value, base_type):
            return base_type(value.GetContainer(), nd_data_type(value.data, True), True)
    return None

def SaveBufferedDictToSerializer(node: BufferedDict, path: str, serializer: Kratos.StreamSerializer, scalar_snapshot: dict, tensor_kinds: dict, node_metadata: dict, keep_alive: list, echo_level: int = 0) -> None:
    """Recursively walks every buffer slot and sub item of a BufferedDict.

    Mirrors the traversal BufferedDict.PrintData/__Info already implements (iterate every buffer
    slot, then recurse into sub items). Kratos.TensorAdaptors.* leaves are saved natively into
    `serializer` (so they stay linked to the Model also saved into the same serializer -- see
    OptimizationProblemRestartOutput.SaveRestartFile); plain picklable leaves go into
    `scalar_snapshot` instead.

    "path" is the BufferedDict.SetValue()-style "/"-separated key of `node` itself (SetValue
    natively resolves such nested keys, see buffered_dict.py). Each leaf's own BufferedDict key is
    `f"{path}/{key}"` (or just `key` at the root); `tensor_kinds`/`scalar_snapshot` are keyed by
    `(bufferdict_key, step_index)` tuples, since the same key can hold different values across
    buffer slots. Tensor leaves additionally get a unique string tag (`tensor_kinds`' value)
    derived from that same pair, since Kratos.Serializer.Save/Load need a string tag.

    `node_metadata` records every sub item's own (buffer_size, clear_buffer_when_advancing), keyed
    by its path -- used on restore to recreate a sub item that no longer has a live counterpart
    (e.g. a control removed from this run's parameters since the checkpoint was written) with its
    original buffer size, instead of BufferedDict.SetValue's path auto-creation, which would default
    it to buffer_size=1.

    Leaves of an unsupported type (e.g. a raw C++-bound object such as a SystemIdentification
    Sensor stored via ComponentDataView's UnBuffered data) are skipped (with a warning) instead of
    aborting the checkpoint, since such objects are typically re-derived from their own inputs
    (e.g. read back from a settings/measurement file) on every run and are not cross-iteration
    state that needs restoring.

    `keep_alive` collects every laundered TensorAdaptor this call creates, and the caller must hold
    it until the whole save (not just this call) is done. Serializer's pointer-identity tracking
    (needed for the Model<->tensor relinking this whole mechanism relies on) keys off raw C++
    addresses; a laundered object that goes out of scope and gets garbage-collected mid-save frees
    its address for reuse by the next allocation, which the serializer would then wrongly treat as
    "the same object, already written" -- silently saving a stale reference instead of real data.
    """
    for step_index in range(node.GetBufferSize()):
        for key, value in node.GetValueItems(step_index).items():
            bufferdict_key = f"{path}/{key}" if path else key
            native_value = _ToNativelySerializable(value) if isinstance(value, _TENSOR_ADAPTOR_BASE_TYPES) else None
            if native_value is not None:
                keep_alive.append(native_value)
                serializer_tag = f"{bufferdict_key}:{step_index}"
                serializer.Save(serializer_tag, native_value)
                tensor_kinds[(bufferdict_key, step_index)] = (serializer_tag, type(native_value).__name__)
            elif isinstance(value, _SCALAR_LEAF_TYPES):
                scalar_snapshot[(bufferdict_key, step_index)] = value
            elif echo_level > 0:
                Kratos.Logger.PrintWarning(
                    "OptimizationProblemRestartOutput",
                    f"Skipping check-pointing of \"{bufferdict_key}\" [ value = {value} ]: does not know how to "
                    f"check-point a value of type \"{type(value).__name__}\". Only tensor adaptors and plain "
                    "int/float/bool/str/None/list/dict values are supported. If this is a new kind of "
                    "cross-iteration state (rather than something re-derived on every run, e.g. from a "
                    "settings/measurement file), either route it through a Kratos.TensorAdaptors.* "
                    "container or store it via ComponentDataView's Buffered/UnBuffered data as one of "
                    "the supported types.")

    for name, sub_item in node.GetSubItems().items():
        sub_path = f"{path}/{name}" if path else name
        node_metadata[sub_path] = (sub_item.GetBufferSize(), sub_item.GetClearBufferWhenAdvancing())
        SaveBufferedDictToSerializer(sub_item, sub_path, serializer, scalar_snapshot, tensor_kinds, node_metadata, keep_alive, echo_level)

class OptimizationProblemRestartOutput:
    """Writes restart checkpoints: SerializerFlags + the whole Kratos.Model + the whole
    OptimizationProblem BufferedDict, into one StreamSerializer, pickled to one file per step.

    Deliberately not a Kratos.OutputProcess: it has exactly one caller, Algorithm.SaveRestartFile()
    (called right after UpdateControl(), so the checkpoint always matches a design the mesh and
    BufferedDict agree on) -- registering this as an output process too would give it a second,
    live entry point and risk writing the same checkpoint twice per step.
    """
    def GetDefaultParameters(self) -> Kratos.Parameters:
        return Kratos.Parameters("""{
            "restart_file_name"     : "Optimization_Restart/restart_<step>.pkl",
            "restart_save_frequency": 1,
            "max_files_to_keep"     : -1,
            "serializer_trace"      : "no_trace",
            "echo_level"            : 0
        }""")

    def __init__(self, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem, model: Kratos.Model):
        parameters.ValidateAndAssignDefaults(self.GetDefaultParameters())

        self.optimization_problem = optimization_problem
        self.model = model
        self.restart_files_path, self.restart_file_name = SplitRestartFileName(parameters["restart_file_name"].GetString())
        self.restart_save_frequency = parameters["restart_save_frequency"].GetInt()
        if self.restart_save_frequency <= 0:
            raise RuntimeError(f"\"restart_save_frequency\" must be > 0. [ restart_save_frequency = {self.restart_save_frequency} ].")
        self.max_files_to_keep = parameters["max_files_to_keep"].GetInt()
        self.serializer_trace = _SERIALIZER_TRACE_TYPES[parameters["serializer_trace"].GetString()]
        self.echo_level = parameters["echo_level"].GetInt()

        if "<step>" not in self.restart_file_name:
            raise RuntimeError(f"\"restart_file_name\" should contain the \"<step>\" placeholder [ restart_file_name = \"{self.restart_file_name}\" ].")

        self.restart_files_path.mkdir(parents=True, exist_ok=True)

    def IsSaveStep(self) -> bool:
        return self.optimization_problem.GetStep() % self.restart_save_frequency == 0

    def SaveRestartFile(self) -> None:
        step = self.optimization_problem.GetStep()

        serializer = Kratos.StreamSerializer(self.serializer_trace)
        flags = Kratos.Serializer.SHALLOW_GLOBAL_POINTERS_SERIALIZATION
        serializer.Set(flags)
        # StreamSerializer's own pickle state does not carry serializer.Set() flags (only the
        # buffer + trace type), so save them as the first record and reapply on load -- same
        # pattern kratos/python_scripts/project.py uses for Model + stage checkpointing.
        serializer.Save("SerializerFlags", flags)

        root_model_parts = [self.model[name] for name in self.model.GetModelPartNames() if "." not in name]

        # NEIGHBOUR_ELEMENTS/etc. (e.g. from the Helmholtz filter) can't survive a Serializer
        # round-trip (nodes load before the entities they point to -- crash on load). They're
        # recomputable, so clear them for this save (covers both the Model save below and any
        # tensor, e.g. "control_field", whose container reaches the same nodes) and restore right
        # after so this run's own continuation is unaffected.
        KratosOA.OptAppModelPartUtils.ClearNeighbourEntitiesData(root_model_parts)
        try:
            # strips the dynamically-created "<OPTIMIZATION_APP_AUTO>"-prefixed sub-model-parts
            # (pure recomputable scratch/cache, e.g. for overhang-angle sensitivity) so they don't
            # bloat/complicate the restart file; they are recreated on demand when next needed.
            KratosOA.OptAppModelPartUtils.RemoveModelPartsWithCommonReferenceEntitiesBetweenReferenceListAndExaminedList(root_model_parts)

            serializer.Save("Model", self.model)

            scalar_snapshot = {}
            tensor_kinds = {}
            node_metadata = {}
            keep_alive = []  # see SaveBufferedDictToSerializer's own docstring
            SaveBufferedDictToSerializer(self.optimization_problem.GetProblemDataContainer(), "", serializer, scalar_snapshot, tensor_kinds, node_metadata, keep_alive, self.echo_level)
        finally:
            KratosOA.OptAppModelPartUtils.RestoreNeighbourEntitiesData()

        file_path = (self.restart_files_path / self.restart_file_name.replace("<step>", str(step))).resolve()
        if self.restart_files_path.resolve() not in file_path.parents:
            raise RuntimeError(f"Resolved restart checkpoint path is outside restart_files_path. [ file_path = \"{file_path}\" ].")

        # write to a temp file in the same directory and rename into place, so a run killed
        # mid-dump never leaves a truncated file that a later "latest" resolution could pick up.
        tmp_file_path = file_path.with_name(file_path.name + ".tmp")
        with open(tmp_file_path, "wb") as file_output:
            pickle.dump({
                "step": step,
                "serializer": serializer,
                "scalars": scalar_snapshot,
                "tensor_kinds": tensor_kinds,
                "node_metadata": node_metadata,
            }, file_output, protocol=pickle.HIGHEST_PROTOCOL)
        os.replace(tmp_file_path, file_path)

        if self.echo_level > 0:
            Kratos.Logger.PrintInfo(self.__class__.__name__, f"Wrote restart checkpoint to \"{file_path}\".")

        self.__PruneOldCheckpoints()

    def __PruneOldCheckpoints(self) -> None:
        if self.max_files_to_keep <= 0:
            return

        checkpoint_glob = self.restart_file_name.replace("<step>", "*")
        checkpoints = sorted(self.restart_files_path.glob(checkpoint_glob), key=lambda p: ParseStepFromFileName(self.restart_file_name, p.name))
        for stale_checkpoint in checkpoints[:-self.max_files_to_keep]:
            stale_checkpoint.unlink()
            if self.echo_level > 0:
                Kratos.Logger.PrintInfo(self.__class__.__name__, f"Removed old restart checkpoint \"{stale_checkpoint}\".")
