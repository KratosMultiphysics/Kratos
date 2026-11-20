"""Sample data structure and (de)serialization helpers for active learning.

Pure Kratos + numpy: this module never imports torch or physicsnemo.
"""

from dataclasses import dataclass, field
from pathlib import Path

import numpy

import KratosMultiphysics as Kratos


@dataclass
class KratosALSample:
    """One active-learning sample: a design-space point plus its labels.

    This is the queue item type flowing through the physicsnemo active
    learning queues (query_queue -> label_queue -> serialize_queue).

    Attributes:
        sample_id: Unique identifier; also names the per-sample case
            directory of the execution backends.
        parameters: The design-space point, as JSON-path -> value overrides
            applied to the case's ProjectParameters (see
            ApplyParameterOverrides).
        fields: Labeled result fields, filled by the label strategy. Keys are
            "<VARIABLE>__<location>", matching DatasetExportProcess npz keys.
        metadata: Free-form bookkeeping (solver time, retries, ...).
    """
    sample_id: str
    parameters: dict = field(default_factory=dict)
    fields: dict = field(default_factory=dict)
    metadata: dict = field(default_factory=dict)

    @property
    def is_labeled(self) -> bool:
        return bool(self.fields)


def ApplyParameterOverrides(project_parameters: Kratos.Parameters, overrides: dict) -> None:
    """Patches a Parameters object with '/'-separated JSON-path overrides.

    Example path: "processes/loads_process_list/0/Parameters/modulus".
    Integer path segments index into array parameters. The leaf must already
    exist and its current type decides which typed setter is used.
    """
    for path, value in overrides.items():
        current = project_parameters
        segments = [s for s in path.split("/") if s != ""]
        if not segments:
            raise ValueError(f"Empty override path for value {value}.")
        for segment in segments[:-1]:
            current = _Descend(current, segment, path)
        leaf_name = segments[-1]
        leaf = _Descend(current, leaf_name, path)
        _SetTypedValue(current, leaf_name, leaf, value, path)


def _Descend(parameters: Kratos.Parameters, segment: str, full_path: str) -> Kratos.Parameters:
    if segment.lstrip("-").isdigit():
        index = int(segment)
        if not parameters.IsArray():
            raise ValueError(f"Path segment \"{segment}\" of \"{full_path}\" indexes a non-array parameter.")
        if not 0 <= index < parameters.size():
            raise ValueError(f"Index {index} of \"{full_path}\" out of range (size {parameters.size()}).")
        return parameters[index]
    if not parameters.Has(segment):
        raise ValueError(f"Path segment \"{segment}\" of \"{full_path}\" does not exist in the parameters.")
    return parameters[segment]


def _SetTypedValue(parent: Kratos.Parameters, name: str, leaf: Kratos.Parameters, value, full_path: str) -> None:
    # Arrays are addressed by index: parent[index] returns a leaf whose typed
    # setters modify the underlying storage in place, so setting on `leaf`
    # covers both the object-key and the array-index cases.
    if leaf.IsBool() and isinstance(value, bool):
        leaf.SetBool(value)
    elif leaf.IsInt() and isinstance(value, int) and not isinstance(value, bool):
        leaf.SetInt(value)
    elif leaf.IsDouble() and isinstance(value, (int, float)) and not isinstance(value, bool):
        leaf.SetDouble(float(value))
    elif leaf.IsString() and isinstance(value, str):
        leaf.SetString(value)
    else:
        raise ValueError(
            f"Type mismatch for override \"{full_path}\": cannot assign {value!r} "
            f"to the existing parameter (value \"{leaf.WriteJsonString()}\").")


def LoadFieldsFromNpzDirectory(path, last_step_only: bool = True):
    """Loads the fields written by DatasetExportProcess from a directory.

    Args:
        path: Directory containing "<prefix>_<step>.npz" files.
        last_step_only: If True (default) only the highest-step file is read;
            otherwise fields get an extra leading step axis (stacked in step
            order).

    Returns:
        (fields, metadata): fields maps "<VARIABLE>__<location>" -> array;
        metadata holds "TIME" and "STEP" (arrays of all read steps).
    """
    path = Path(path)
    files = sorted(path.glob("*.npz"), key=lambda p: int(p.stem.rsplit("_", 1)[-1]))
    if not files:
        raise FileNotFoundError(f"No .npz result files found in \"{path}\".")
    if last_step_only:
        files = files[-1:]

    per_step = []
    for f in files:
        with numpy.load(f) as data:
            per_step.append({key: numpy.array(data[key]) for key in data.files})

    metadata = {
        "TIME": numpy.array([step.pop("TIME") for step in per_step]),
        "STEP": numpy.array([step.pop("STEP") for step in per_step]),
    }
    if last_step_only:
        return per_step[0], metadata
    keys = per_step[0].keys()
    return {key: numpy.stack([step[key] for step in per_step]) for key in keys}, metadata
