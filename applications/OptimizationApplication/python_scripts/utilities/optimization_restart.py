import importlib
import os
import pickle
from pathlib import Path

import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.buffered_dict import BufferedDict
from KratosMultiphysics.OptimizationApplication.utilities.restart_file_naming import SplitRestartFileName, ListCheckpoints

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

_SCALAR_LEAF_TYPES = (int, float, bool, str, list, dict, type(None), numpy.ndarray)

def CheckSerialRestart() -> None:
    """Raises if running with more than one rank: optimization restart is serial-only for now."""
    if Kratos.DataCommunicator.GetDefault().Size() > 1:
        raise RuntimeError("Optimization restart (restart_settings.save_restart/load_restart) supports serial runs only.")

def GetImportedApplicationNames() -> 'list[str]':
    """Returns the names of the installed Kratos applications registered with the Kernel.

    An application's Kernel name equals its KratosMultiphysics.<name> module name, so every package
    next to KratosMultiphysics/__init__.py that the Kernel reports as imported is returned.
    """
    kernel = Kratos.KratosGlobals.Kernel
    kratos_package_path = Path(Kratos.__file__).parent
    return sorted(path.name for path in kratos_package_path.iterdir() if path.is_dir() and kernel.IsImported(path.name))

def ImportApplications(application_names: 'list[str]') -> None:
    """Imports KratosMultiphysics.<name> for every name, so that the elements, conditions and
    constitutive laws of a checkpointed Model are registered before the Serializer loads it."""
    for application_name in application_names:
        importlib.import_module(f"KratosMultiphysics.{application_name}")

def CreateTensorAdaptorLoadTarget(type_name: str):
    """Returns a default-constructed tensor adaptor of type type_name, to Load a saved one into."""
    for module in _TENSOR_ADAPTOR_MODULES:
        if hasattr(module, type_name):
            return getattr(module, type_name)()
    raise RuntimeError(f"Unknown restart tensor adaptor type \"{type_name}\".")

def SaveBufferedDictToSerializer(node: BufferedDict, path: str, serializer: Kratos.StreamSerializer, scalar_snapshot: dict, tensor_kinds: dict, node_metadata: dict, warned_keys: 'set[str]') -> None:
    """Recursively saves a BufferedDict node and all its sub items for a restart checkpoint.

    Tensor adaptors go into serializer (relinked to the saved Model on load) and their
    (tag, type name) into tensor_kinds; plain values (int/float/bool/str/None/list/dict/numpy array) go into
    scalar_snapshot. Both are keyed by (path/key, step_index). node_metadata gets each sub node's
    (buffer_size, clear_buffer_when_advancing). Values of any other type are skipped with a warning,
    once per key (warned_keys holds the keys already warned about).
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
            elif bufferdict_key not in warned_keys:
                warned_keys.add(bufferdict_key)
                Kratos.Logger.PrintWarning(
                    "OptimizationProblemRestartOutput",
                    f"Skipping check-pointing of \"{bufferdict_key}\" [ value = {value} ]: does not know how to "
                    f"check-point a value of type \"{type(value).__name__}\". Only tensor adaptors and plain "
                    "int/float/bool/str/None/list/dict/numpy array values are supported. If this is a new kind of "
                    "cross-iteration state (rather than something re-derived on every run, e.g. from a "
                    "settings/measurement file), either route it through a Kratos.TensorAdaptors.* "
                    "container or store it via ComponentDataView's Buffered/UnBuffered data as one of "
                    "the supported types.")

    for name, sub_item in node.GetSubItems().items():
        sub_path = f"{path}/{name}" if path else name
        node_metadata[sub_path] = (sub_item.GetBufferSize(), sub_item.GetClearBufferWhenAdvancing())
        SaveBufferedDictToSerializer(sub_item, sub_path, serializer, scalar_snapshot, tensor_kinds, node_metadata, warned_keys)

class OptimizationProblemRestartOutput:
    """Writes optimization restart checkpoints.

    Each checkpoint is one pickled file per save step holding a StreamSerializer with the
    "SerializerFlags" record, the whole Kratos.Model and every tensor of the optimization problem's
    BufferedDict, plus the BufferedDict's plain scalars and node metadata (buffer sizes) and the
    imported applications (imported again before the Model is loaded) next to it.
    OptimizationAnalysis restores it verbatim: the Model before any model part controller is
    created, the BufferedDict in Initialize() after every component has sized its data buffers.

    This is deliberately not a Kratos.OutputProcess. It is installed by OptimizationAnalysis via
    OptimizationProblem.SetRestartOutput() and its only caller is Algorithm.SaveRestartFile(), right
    after UpdateControl(), so the saved mesh and BufferedDict agree on the same design.

    Files are written atomically (temp file + rename), and with "max_files_to_keep" > 0 the oldest
    checkpoints (by parsed step) are deleted.
    """
    @classmethod
    def GetDefaultParameters(cls) -> Kratos.Parameters:
        return Kratos.Parameters("""{
            "restart_file_name"     : "Optimization_Restart/restart_<step>.pkl",
            "restart_save_frequency": 1,
            "max_files_to_keep"     : -1,
            "serializer_trace"      : "no_trace",
            "echo_level"            : 0
        }""")

    def __init__(self, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem, model: Kratos.Model):
        CheckSerialRestart()
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
        self.__warned_keys: 'set[str]' = set()

        if "<step>" not in self.restart_file_name:
            raise RuntimeError(f"\"restart_file_name\" should contain the \"<step>\" placeholder [ restart_file_name = \"{self.restart_file_name}\" ].")

        self.restart_files_path.mkdir(parents=True, exist_ok=True)

    def IsSaveStep(self) -> bool:
        """Whether the current step is a multiple of "restart_save_frequency"."""
        return self.optimization_problem.GetStep() % self.restart_save_frequency == 0

    def SaveRestartFile(self) -> None:
        """Writes the checkpoint of the current step and prunes old checkpoints.

        Neighbour entity data and the "<OPTIMIZATION_APP_AUTO>" sub model parts are removed before
        saving (both are recomputable); the neighbour data is restored afterwards, even on failure.
        """
        step = self.optimization_problem.GetStep()

        serializer = Kratos.StreamSerializer(self.serializer_trace)
        flags = Kratos.Serializer.SHALLOW_GLOBAL_POINTERS_SERIALIZATION
        serializer.Set(flags)
        # StreamSerializer's pickle state drops serializer.Set() flags, so save them as the first
        # record and reapply them on load (as kratos/python_scripts/project.py does).
        serializer.Save("SerializerFlags", flags)

        root_model_parts = [self.model[name] for name in self.model.GetModelPartNames() if "." not in name]

        # NEIGHBOUR_* pointers can't survive a Serializer round-trip; they are recomputable, so clear
        # them for the save and restore them right after.
        neighbour_backup = KratosOA.OptAppModelPartUtils.ClearNeighbourEntitiesData(root_model_parts)
        try:
            # "<OPTIMIZATION_APP_AUTO>" sub model parts are recomputable cache, so don't checkpoint them
            KratosOA.OptAppModelPartUtils.RemoveModelPartsWithCommonReferenceEntitiesBetweenReferenceListAndExaminedList(root_model_parts)

            serializer.Save("Model", self.model)

            scalar_snapshot = {}
            tensor_kinds = {}
            node_metadata = {}
            SaveBufferedDictToSerializer(self.optimization_problem.GetProblemDataContainer(), "", serializer, scalar_snapshot, tensor_kinds, node_metadata, self.__warned_keys)
        finally:
            KratosOA.OptAppModelPartUtils.RestoreNeighbourEntitiesData(neighbour_backup)

        file_path = (self.restart_files_path / self.restart_file_name.replace("<step>", str(step))).resolve()
        if self.restart_files_path.resolve() not in file_path.parents:
            raise RuntimeError(f"Resolved restart checkpoint path is outside restart_files_path. [ file_path = \"{file_path}\" ].")

        # write to a temp file in the same directory and rename into place, so a run killed
        # mid-dump never leaves a truncated file that a later "latest" resolution could pick up.
        tmp_file_path = file_path.with_name(file_path.name + ".tmp")
        with open(tmp_file_path, "wb") as file_output:
            pickle.dump({
                "step": step,
                "applications": GetImportedApplicationNames(),
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

        checkpoints = ListCheckpoints(self.restart_files_path, self.restart_file_name)
        for _, stale_checkpoint in checkpoints[:-self.max_files_to_keep]:
            stale_checkpoint.unlink()
            if self.echo_level > 0:
                Kratos.Logger.PrintInfo(self.__class__.__name__, f"Removed old restart checkpoint \"{stale_checkpoint}\".")
