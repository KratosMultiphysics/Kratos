import pickle
from pathlib import Path

import numpy
import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.buffered_dict import BufferedDict
from KratosMultiphysics.OptimizationApplication.utilities.restart_file_naming import GetModelPartRestartFileBaseName, ParseStepFromFileName, SplitRestartFileName

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

_TENSOR_ADAPTOR_TYPES = (
    Kratos.TensorAdaptors.DoubleTensorAdaptor,
    Kratos.TensorAdaptors.IntTensorAdaptor,
    Kratos.TensorAdaptors.BoolTensorAdaptor,
)

_SCALAR_LEAF_TYPES = (int, float, bool, str, list, dict, type(None))

class _UnsupportedLeafType(Exception):
    pass

def ConvertLeafToRestartData(value) -> dict:
    """Converts a single BufferedDict leaf value to a plain, pickle-safe representation.

    Kratos.TensorAdaptors.* are C++-bound and are not picklable directly, so their data is
    copied out to plain numpy arrays here. Everything else stored in the OptimizationProblem's
    BufferedDict tree is expected to already be a plain, picklable Python value.

    Raises _UnsupportedLeafType for anything else (e.g. a raw C++-bound object such as a
    SystemIdentificationApplication Sensor stored via ComponentDataView's UnBuffered data); the
    caller skips those leaves rather than failing the whole checkpoint, since such objects are
    typically re-derived from their own inputs (e.g. read back from a settings/measurement file)
    on every run and are not cross-iteration state that needs restoring.
    """
    if isinstance(value, Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor):
        return {"kind": "combined_tensor", "parts": [numpy.array(ta.data, copy=True) for ta in value.GetTensorAdaptors()]}
    if isinstance(value, _TENSOR_ADAPTOR_TYPES):
        return {"kind": "tensor", "array": numpy.array(value.data, copy=True)}
    if isinstance(value, _SCALAR_LEAF_TYPES):
        return {"kind": "scalar", "value": value}

    raise _UnsupportedLeafType(type(value).__name__)

def SnapshotBufferedDict(node: BufferedDict, echo_level: int = 0) -> dict:
    """Recursively snapshots every buffer slot and sub item of a BufferedDict into a plain dict.

    Mirrors the traversal BufferedDict.PrintData/__Info already implements (iterate every
    buffer slot, then recurse into sub items), just building a nested dict instead of a string.

    Leaves of an unsupported type are omitted from the snapshot (with a warning) instead of
    aborting the checkpoint; see ConvertLeafToRestartData.
    """
    slots = []
    for step_index in range(node.GetBufferSize()):
        slot = {}
        for key, value in node.GetValueItems(step_index).items():
            try:
                slot[key] = ConvertLeafToRestartData(value)
            except _UnsupportedLeafType as exc:
                if echo_level > 0:
                    Kratos.Logger.PrintWarning(
                        "OptimizationProblemRestartOutputProcess",
                        f"Skipping check-pointing of \"{key}\" [ value = {value} ]: does not know how to "
                        f"check-point a value of type \"{exc}\". Only tensor adaptors and plain "
                        "int/float/bool/str/None/list/dict values are supported. If this is a new kind of "
                        "cross-iteration state (rather than something re-derived on every run, e.g. from a "
                        "settings/measurement file), either route it through a Kratos.TensorAdaptors.* "
                        "container or store it via ComponentDataView's Buffered/UnBuffered data as one of "
                        "the supported types.")
        slots.append(slot)

    return {
        "slots": slots,
        "sub_items": {name: SnapshotBufferedDict(sub_item, echo_level) for name, sub_item in node.GetSubItems().items()},
    }

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

        payload = {
            "step": step,
            "data": SnapshotBufferedDict(self.optimization_problem.GetProblemDataContainer(), self.echo_level),
        }

        file_path = (self.restart_files_path / self.restart_file_name.replace("<step>", str(step))).resolve()
        if self.restart_files_path.resolve() not in file_path.parents:
            raise RuntimeError(f"Resolved restart checkpoint path is outside restart_files_path. [ file_path = \"{file_path}\" ].")
        with open(file_path, "wb") as file_output:
            pickle.dump(payload, file_output, protocol=pickle.HIGHEST_PROTOCOL)

        if self.echo_level > 0:
            Kratos.Logger.PrintInfo(self.__class__.__name__, f"Wrote restart checkpoint to \"{file_path}\".")

        if self.save_model_parts and self.list_of_model_parts:
            self.__SaveModelParts(step)

        self.__PruneOldCheckpoints()

    def __SaveModelParts(self, step: int) -> None:
        if self.clean_before_save:
            # strips the dynamically-created "<OPTIMIZATION_APP_AUTO>"-prefixed sub-model-parts
            # (pure recomputable scratch/cache, e.g. for overhang-angle sensitivity) so they don't
            # bloat/complicate the restart file; they are recreated on demand when next needed.
            KratosOA.OptAppModelPartUtils.RemoveModelPartsWithCommonReferenceEntitiesBetweenReferenceListAndExaminedList(self.list_of_model_parts)

        for model_part in self.list_of_model_parts:
            base_name = GetModelPartRestartFileBaseName(self.restart_file_name, step, model_part.Name)
            mp_file_path = (self.restart_files_path / base_name).resolve()
            if self.restart_files_path.resolve() not in mp_file_path.parents:
                raise RuntimeError(f"Resolved restart checkpoint path is outside restart_files_path. [ file_path = \"{mp_file_path}\" ].")
            serializer = Kratos.FileSerializer(str(mp_file_path), self.serializer_trace)
            serializer.Set(Kratos.Serializer.SHALLOW_GLOBAL_POINTERS_SERIALIZATION)
            serializer.Save(model_part.Name, model_part)
            if self.echo_level > 0:
                Kratos.Logger.PrintInfo(self.__class__.__name__, f"Wrote restart checkpoint for model part \"{model_part.Name}\" to \"{mp_file_path}.rest\".")

    def __PruneOldCheckpoints(self) -> None:
        if self.max_files_to_keep <= 0:
            return

        checkpoint_glob = self.restart_file_name.replace("<step>", "*")
        checkpoints = sorted(self.restart_files_path.glob(checkpoint_glob), key=lambda p: p.stat().st_mtime)
        for stale_checkpoint in checkpoints[:-self.max_files_to_keep]:
            stale_step = ParseStepFromFileName(self.restart_file_name, stale_checkpoint.name)
            stale_checkpoint.unlink()
            if self.echo_level > 0:
                Kratos.Logger.PrintInfo(self.__class__.__name__, f"Removed old restart checkpoint \"{stale_checkpoint}\".")

            if self.save_model_parts:
                for model_part in self.list_of_model_parts:
                    base_name = GetModelPartRestartFileBaseName(self.restart_file_name, stale_step, model_part.Name)
                    rest_path = self.restart_files_path / f"{base_name}.rest"
                    if rest_path.is_file():
                        rest_path.unlink()
                        if self.echo_level > 0:
                            Kratos.Logger.PrintInfo(self.__class__.__name__, f"Removed old restart checkpoint for model part \"{model_part.Name}\": \"{rest_path}\".")
