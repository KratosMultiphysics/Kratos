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

# isinstance against the 3 base types also matches every subtype (combined, Variable, Properties, ...).
_TENSOR_ADAPTOR_BASE_TYPES = (
    Kratos.TensorAdaptors.DoubleTensorAdaptor,
    Kratos.TensorAdaptors.IntTensorAdaptor,
    Kratos.TensorAdaptors.BoolTensorAdaptor,
)

# modules holding every tensor adaptor type registered with Kratos::Serializer.
_TENSOR_ADAPTOR_MODULES = (Kratos.TensorAdaptors, KratosOA.TensorAdaptors)

_SCALAR_LEAF_TYPES = (int, float, bool, str, list, dict, type(None))

def CreateTensorAdaptorLoadTarget(type_name: str):
    """Returns a default-constructed tensor adaptor of type type_name, to Load a saved one into."""
    for module in _TENSOR_ADAPTOR_MODULES:
        if hasattr(module, type_name):
            return getattr(module, type_name)()
    raise RuntimeError(f"Unknown restart tensor adaptor type \"{type_name}\".")

def SaveBufferedDictToSerializer(node: BufferedDict, path: str, serializer: Kratos.StreamSerializer, scalar_snapshot: dict, tensor_kinds: dict, node_metadata: dict, echo_level: int = 0) -> None:
    """Recursively walks every buffer slot and sub item of a BufferedDict.

    Mirrors the traversal BufferedDict.PrintData/__Info already implements (iterate every buffer
    slot, then recurse into sub items). Kratos.TensorAdaptors.* leaves are saved natively into
    `serializer` as themselves (so they stay linked to the Model also saved into the same
    serializer -- see OptimizationProblemRestartOutput.SaveRestartFile); plain picklable leaves go
    into `scalar_snapshot` instead.

    "path" is the BufferedDict.SetValue()-style "/"-separated key of `node` itself (SetValue
    natively resolves such nested keys, see buffered_dict.py). Each leaf's own BufferedDict key is
    `f"{path}/{key}"` (or just `key` at the root); `tensor_kinds`/`scalar_snapshot` are keyed by
    `(bufferdict_key, step_index)` tuples, since the same key can hold different values across
    buffer slots. Tensor leaves additionally get a unique string tag derived from that same pair,
    since Kratos.Serializer.Save/Load need a string tag; `tensor_kinds` maps the pair to
    (tag, exact type name), which is all CreateTensorAdaptorLoadTarget needs on restore.

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
    """
    for step_index in range(node.GetBufferSize()):
        for key, value in node.GetValueItems(step_index).items():
            bufferdict_key = f"{path}/{key}" if path else key
            if isinstance(value, _TENSOR_ADAPTOR_BASE_TYPES):
                serializer_tag = f"{bufferdict_key}:{step_index}"
                serializer.Save(serializer_tag, value)
                tensor_kinds[(bufferdict_key, step_index)] = (serializer_tag, type(value).__name__)
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
        SaveBufferedDictToSerializer(sub_item, sub_path, serializer, scalar_snapshot, tensor_kinds, node_metadata, echo_level)

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
            SaveBufferedDictToSerializer(self.optimization_problem.GetProblemDataContainer(), "", serializer, scalar_snapshot, tensor_kinds, node_metadata, self.echo_level)
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
