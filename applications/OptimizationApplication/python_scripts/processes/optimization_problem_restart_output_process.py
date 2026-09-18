import pickle

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.buffered_dict import BufferedDict
from KratosMultiphysics.OptimizationApplication.utilities.restart_file_naming import SplitRestartFileName

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

# isinstance against just the 3 base types also matches the corresponding combined subtype, since
# Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor (etc.) is bound as a Python subclass of
# Kratos.TensorAdaptors.DoubleTensorAdaptor (etc.) -- see kratos/python/add_tensor_adaptors_to_python.cpp.
# It also matches further-derived, application-specific subtypes not registered with the
# Serializer (e.g. SystemIdentificationApplication's PropertiesVariableTensorAdaptor) -- those get
# "laundered" into a plain base-type instance before saving; see _ToNativelySerializable.
_TENSOR_ADAPTOR_BASE_TYPES = (
    Kratos.TensorAdaptors.DoubleTensorAdaptor,
    Kratos.TensorAdaptors.IntTensorAdaptor,
    Kratos.TensorAdaptors.BoolTensorAdaptor,
)

# the only types actually registered with Kratos::Serializer (kratos/sources/kratos_application.cpp)
# and therefore safe to save/load directly, as themselves.
_EXACT_NATIVE_TENSOR_ADAPTOR_TYPES = (
    Kratos.TensorAdaptors.DoubleTensorAdaptor,
    Kratos.TensorAdaptors.IntTensorAdaptor,
    Kratos.TensorAdaptors.BoolTensorAdaptor,
    Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor,
    Kratos.TensorAdaptors.IntCombinedTensorAdaptor,
    Kratos.TensorAdaptors.BoolCombinedTensorAdaptor,
)

# (base tensor adaptor type, matching NDData type) per dtype, used by _ToNativelySerializable.
_DTYPE_TENSOR_ADAPTOR_TYPES = (
    (Kratos.TensorAdaptors.DoubleTensorAdaptor, Kratos.DoubleNDData),
    (Kratos.TensorAdaptors.IntTensorAdaptor, Kratos.IntNDData),
    (Kratos.TensorAdaptors.BoolTensorAdaptor, Kratos.BoolNDData),
)

_SCALAR_LEAF_TYPES = (int, float, bool, str, list, dict, type(None))

def _ToNativelySerializable(value):
    """Returns a TensorAdaptor safe to save/load directly through Kratos.Serializer.

    `value` itself if its exact type is already one of the few registered with the Serializer
    (kratos_application.cpp); otherwise "launders" it into a freshly-constructed plain base-type
    instance wrapping the *same* container and a copy of the *same* data -- dropping whatever
    extra, subclass-specific state a further-derived type carries (e.g.
    PropertiesVariableTensorAdaptor's mapped Properties Variable), since restart only needs to
    round-trip the current values and which container they belong to, not that subclass's own
    CollectData/StoreData behaviour.

    Returns None if `value` has no container to launder through (HasContainer() is False) and
    isn't itself one of the registered types -- the caller falls back to skipping it, with a
    warning, same as any other genuinely unsupported leaf.
    """
    if type(value) in _EXACT_NATIVE_TENSOR_ADAPTOR_TYPES:
        return value
    if not value.HasContainer():
        return None
    for base_type, nd_data_type in _DTYPE_TENSOR_ADAPTOR_TYPES:
        if isinstance(value, base_type):
            return base_type(value.GetContainer(), nd_data_type(value.data, True), True)
    return None

def SaveBufferedDictToSerializer(node: BufferedDict, path: str, serializer: Kratos.StreamSerializer, scalar_snapshot: dict, tensor_kinds: dict, echo_level: int = 0) -> None:
    """Recursively walks every buffer slot and sub item of a BufferedDict.

    Mirrors the traversal BufferedDict.PrintData/__Info already implements (iterate every buffer
    slot, then recurse into sub items). Kratos.TensorAdaptors.* leaves are saved natively into
    `serializer` (so they stay linked to whichever ModelPart(s) they were built from, also saved
    into the same serializer -- see OptimizationProblemRestartOutputProcess.PrintOutput); plain
    picklable leaves go into `scalar_snapshot` instead.

    "path" is the BufferedDict.SetValue()-style "/"-separated key of `node` itself (SetValue
    natively resolves such nested keys, see buffered_dict.py). Each leaf's own BufferedDict key is
    `f"{path}/{key}"` (or just `key` at the root); `tensor_kinds`/`scalar_snapshot` are keyed by
    `(bufferdict_key, step_index)` tuples, since the same key can hold different values across
    buffer slots. Tensor leaves additionally get a unique string tag (`tensor_kinds`' value)
    derived from that same pair, since Kratos.Serializer.Save/Load need a string tag.

    Leaves of an unsupported type (e.g. a raw C++-bound object such as a SystemIdentification
    Sensor stored via ComponentDataView's UnBuffered data) are skipped (with a warning) instead of
    aborting the checkpoint, since such objects are typically re-derived from their own inputs
    (e.g. read back from a settings/measurement file) on every run and are not cross-iteration
    state that needs restoring.
    """
    for step_index in range(node.GetBufferSize()):
        for key, value in node.GetValueItems(step_index).items():
            bufferdict_key = f"{path}/{key}" if path else key
            native_value = _ToNativelySerializable(value) if isinstance(value, _TENSOR_ADAPTOR_BASE_TYPES) else None
            if native_value is not None:
                serializer_tag = f"{bufferdict_key}:{step_index}"
                serializer.Save(serializer_tag, native_value)
                tensor_kinds[(bufferdict_key, step_index)] = (serializer_tag, type(native_value).__name__)
            elif isinstance(value, _SCALAR_LEAF_TYPES):
                scalar_snapshot[(bufferdict_key, step_index)] = value
            elif echo_level > 0:
                Kratos.Logger.PrintWarning(
                    "OptimizationProblemRestartOutputProcess",
                    f"Skipping check-pointing of \"{bufferdict_key}\" [ value = {value} ]: does not know how to "
                    f"check-point a value of type \"{type(value).__name__}\". Only tensor adaptors and plain "
                    "int/float/bool/str/None/list/dict values are supported. If this is a new kind of "
                    "cross-iteration state (rather than something re-derived on every run, e.g. from a "
                    "settings/measurement file), either route it through a Kratos.TensorAdaptors.* "
                    "container or store it via ComponentDataView's Buffered/UnBuffered data as one of "
                    "the supported types.")

    for name, sub_item in node.GetSubItems().items():
        SaveBufferedDictToSerializer(sub_item, f"{path}/{name}" if path else name, serializer, scalar_snapshot, tensor_kinds, echo_level)

class OptimizationProblemRestartOutputProcess(Kratos.OutputProcess):
    def GetDefaultParameters(self) -> Kratos.Parameters:
        return Kratos.Parameters("""{
            "restart_file_name"     : "Optimization_Restart/restart_<step>.pkl",
            "restart_save_frequency": 1,
            "max_files_to_keep"     : -1,
            "echo_level"            : 0,
            "model_parts_settings"  : {
                "save_model_parts" : true,
                "serializer_trace" : "no_trace",
                "clean_before_save": true
            }
        }""")

    def __init__(self, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem, list_of_model_parts: 'list[Kratos.ModelPart]'):
        Kratos.OutputProcess.__init__(self)
        parameters.ValidateAndAssignDefaults(self.GetDefaultParameters())

        self.optimization_problem = optimization_problem
        self.restart_files_path, self.restart_file_name = SplitRestartFileName(parameters["restart_file_name"].GetString())
        self.restart_save_frequency = parameters["restart_save_frequency"].GetInt()
        if self.restart_save_frequency <= 0:
            raise RuntimeError(f"\"restart_save_frequency\" must be > 0. [ restart_save_frequency = {self.restart_save_frequency} ].")
        self.max_files_to_keep = parameters["max_files_to_keep"].GetInt()
        self.echo_level = parameters["echo_level"].GetInt()

        if "<step>" not in self.restart_file_name:
            raise RuntimeError(f"\"restart_file_name\" should contain the \"<step>\" placeholder [ restart_file_name = \"{self.restart_file_name}\" ].")

        model_parts_settings = parameters["model_parts_settings"]
        self.save_model_parts = model_parts_settings["save_model_parts"].GetBool()
        self.clean_before_save = model_parts_settings["clean_before_save"].GetBool()
        self.serializer_trace = _SERIALIZER_TRACE_TYPES[model_parts_settings["serializer_trace"].GetString()]
        self.list_of_model_parts = list_of_model_parts

        self.restart_files_path.mkdir(parents=True, exist_ok=True)

    def IsOutputStep(self) -> bool:
        return self.optimization_problem.GetStep() % self.restart_save_frequency == 0

    def PrintOutput(self) -> None:
        step = self.optimization_problem.GetStep()

        serializer = Kratos.StreamSerializer(self.serializer_trace)
        serializer.Set(Kratos.Serializer.SHALLOW_GLOBAL_POINTERS_SERIALIZATION)

        # NEIGHBOUR_ELEMENTS/etc. (e.g. from the Helmholtz filter) can't survive a Serializer
        # round-trip (nodes load before the entities they point to -- crash on load). They're
        # recomputable, so clear them for this save (covers both the model-part save below and any
        # tensor, e.g. "control_field", whose container reaches the same nodes) and restore right
        # after so this run's own continuation is unaffected.
        KratosOA.OptAppModelPartUtils.ClearNeighbourEntitiesData(self.list_of_model_parts)
        try:
            model_part_names = []
            if self.save_model_parts and self.list_of_model_parts:
                if self.clean_before_save:
                    # strips the dynamically-created "<OPTIMIZATION_APP_AUTO>"-prefixed sub-model-parts
                    # (pure recomputable scratch/cache, e.g. for overhang-angle sensitivity) so they
                    # don't bloat/complicate the restart file; they are recreated on demand when next
                    # needed.
                    KratosOA.OptAppModelPartUtils.RemoveModelPartsWithCommonReferenceEntitiesBetweenReferenceListAndExaminedList(self.list_of_model_parts)
                for model_part in self.list_of_model_parts:
                    serializer.Save(model_part.Name, model_part)
                    model_part_names.append(model_part.Name)

            scalar_snapshot = {}
            tensor_kinds = {}
            SaveBufferedDictToSerializer(self.optimization_problem.GetProblemDataContainer(), "", serializer, scalar_snapshot, tensor_kinds, self.echo_level)
        finally:
            KratosOA.OptAppModelPartUtils.RestoreNeighbourEntitiesData()

        file_path = (self.restart_files_path / self.restart_file_name.replace("<step>", str(step))).resolve()
        if self.restart_files_path.resolve() not in file_path.parents:
            raise RuntimeError(f"Resolved restart checkpoint path is outside restart_files_path. [ file_path = \"{file_path}\" ].")
        with open(file_path, "wb") as file_output:
            pickle.dump({
                "step": step,
                "serializer": serializer,
                "model_part_names": model_part_names,
                "scalars": scalar_snapshot,
                "tensor_kinds": tensor_kinds,
            }, file_output, protocol=pickle.HIGHEST_PROTOCOL)

        if self.echo_level > 0:
            Kratos.Logger.PrintInfo(self.__class__.__name__, f"Wrote restart checkpoint to \"{file_path}\".")

        self.__PruneOldCheckpoints()

    def __PruneOldCheckpoints(self) -> None:
        if self.max_files_to_keep <= 0:
            return

        checkpoint_glob = self.restart_file_name.replace("<step>", "*")
        checkpoints = sorted(self.restart_files_path.glob(checkpoint_glob), key=lambda p: p.stat().st_mtime)
        for stale_checkpoint in checkpoints[:-self.max_files_to_keep]:
            stale_checkpoint.unlink()
            if self.echo_level > 0:
                Kratos.Logger.PrintInfo(self.__class__.__name__, f"Removed old restart checkpoint \"{stale_checkpoint}\".")
