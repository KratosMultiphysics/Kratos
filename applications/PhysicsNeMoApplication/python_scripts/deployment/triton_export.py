"""Exports a trained surrogate as a Triton Inference Server model repository.

Turns a checkpoint into inference-as-a-service: Triton loads a directory tree

    <repository>/<model_name>/config.pbtxt
    <repository>/<model_name>/<version>/model.onnx      (or model.pt)

and serves it over HTTP/gRPC, which TritonInferenceProcess then calls from
inside the solution loop instead of loading a local checkpoint.

Two deliberate choices worth knowing:

- **max_batch_size is 0.** Triton's batching axis assumes the leading
  dimension is a sample index it may freely stack; here the leading dimension
  is the ENTITY (node/element) count of one Kratos case. Declaring 0 disables
  batching and makes the declared `dims` the full tensor shape, with -1 on the
  entity axis so any mesh size is accepted.
- **ONNX is exported through torch.onnx.export, not physicsnemo's
  export_to_onnx_stream.** The upstream helper exposes no dynamic_axes,
  input_names or output_names, so it would freeze the entity count into the
  graph - useless for a mesh-size-agnostic service. (Its double-execution
  bug was fixed in physicsnemo 2.2; the missing dynamic axes were not.)
  training_utils.ExportOnnxModel remains the right tool for the fixed-size
  local artifact OnnxInferenceProcess consumes.

Tensor names and widths come from the model card
(model_registry.LoadModelCard) so the served names match the fields the
Kratos-side process gathers.

torch is imported lazily; nothing here needs Triton itself to be installed.
"""

from pathlib import Path

import numpy

import KratosMultiphysics as Kratos
from KratosMultiphysics.PhysicsNeMoApplication.deployment import model_registry
_PLATFORMS = {"onnx": "onnxruntime_onnx", "torchscript": "pytorch_libtorch"}
_MODEL_FILES = {"onnx": "model.onnx", "torchscript": "model.pt"}
_TRITON_DTYPES = {"float32": "TYPE_FP32", "float64": "TYPE_FP64", "float16": "TYPE_FP16"}


def _TryImportTorch():
    try:
        import torch
        return torch
    except ImportError as e:
        raise ImportError(
            "PhysicsNeMoApplication.triton_export requires torch, which could not be "
            "imported. Install it with e.g. 'pip install torch'.") from e


def _FormatTensorBlock(keyword: str, tensor: dict) -> str:
    dims = ", ".join(str(int(d)) for d in tensor["dims"])
    return (f"{keyword} [\n"
            f"  {{\n"
            f"    name: \"{tensor['name']}\"\n"
            f"    data_type: {tensor['data_type']}\n"
            f"    dims: [ {dims} ]\n"
            f"  }}\n"
            f"]")


def MakeTritonConfig(model_name: str, platform: str, inputs, outputs,
                     max_batch_size: int = 0, instance_group=None,
                     dynamic_batching: bool = False) -> str:
    """Renders a Triton `config.pbtxt`.

    Pure text generation - no torch, no Triton - so a repository's contract
    can be checked anywhere.

    Args:
        model_name: The served model name (the repository sub-directory).
        platform: Triton platform string, e.g. "onnxruntime_onnx".
        inputs / outputs: sequences of {"name", "data_type", "dims"} dicts;
            `dims` entries may be -1 for a free axis.
        max_batch_size: 0 (default) disables Triton batching - correct for
            per-entity tensors, whose leading axis is the mesh size rather
            than a stackable sample index.
        instance_group: Optional list of {"count", "kind"} dicts, e.g.
            [{"count": 2, "kind": "KIND_GPU"}].
        dynamic_batching: Emit an empty `dynamic_batching {}` block. Only
            meaningful with max_batch_size > 0.

    Returns:
        The config.pbtxt text.
    """
    if max_batch_size < 0:
        raise ValueError(f"max_batch_size must be >= 0, got {max_batch_size}.")
    if dynamic_batching and max_batch_size == 0:
        raise ValueError(
            "dynamic_batching requires max_batch_size > 0; with per-entity tensors "
            "(max_batch_size 0) Triton has no batch axis to fill.")
    if not inputs or not outputs:
        raise ValueError("A Triton config needs at least one input and one output.")

    blocks = [f"name: \"{model_name}\"",
              f"platform: \"{platform}\"",
              f"max_batch_size: {int(max_batch_size)}"]
    blocks.extend(_FormatTensorBlock("input", tensor) for tensor in inputs)
    blocks.extend(_FormatTensorBlock("output", tensor) for tensor in outputs)

    if instance_group:
        entries = "\n".join(
            f"  {{\n    count: {int(entry.get('count', 1))}\n"
            f"    kind: {entry.get('kind', 'KIND_AUTO')}\n  }}"
            for entry in instance_group)
        blocks.append(f"instance_group [\n{entries}\n]")
    if dynamic_batching:
        blocks.append("dynamic_batching { }")

    return "\n".join(blocks) + "\n"


def _CardWidths(card, key):
    """Per-field widths from a model card's field list, or None when unknown."""
    if not card or key not in card:
        return None
    names = []
    for entry in card[key]:
        variable_name = entry.get("variable_name")
        if variable_name is None:
            return None
        names.append(variable_name)
    return names


def ExportTritonModelRepository(model, sample_inputs, settings: Kratos.Parameters):
    """Writes a Triton model repository for a trained model.

    Args:
        model: The torch model to serve. Its forward must take one
            (n_entities, total_input_width) tensor and return one
            (n_entities, total_output_width) tensor - the contract the
            deployment processes already use.
        sample_inputs: An example input tensor (or tuple of one) giving the
            widths; its entity axis is exported as free.
        settings: Kratos Parameters:
            {
                "repository_path" : "triton_models",
                "model_name"      : "kratos_surrogate",
                "model_version"   : 1,
                "format"          : "onnx" | "torchscript",
                "card_file"       : "",         // model card to reuse/copy
                "input_name"      : "input",    // when no card names it
                "output_name"     : "output",
                "max_batch_size"  : 0,
                "instance_count"  : 0,          // 0 = omit instance_group
                "instance_kind"   : "KIND_AUTO",
                "opset_version"   : 17
            }

    Returns:
        The path of the written config.pbtxt.
    """
    torch = _TryImportTorch()

    defaults = Kratos.Parameters("""{
        "repository_path" : "triton_models",
        "model_name"      : "kratos_surrogate",
        "model_version"   : 1,
        "format"          : "onnx",
        "card_file"       : "",
        "input_name"      : "input",
        "output_name"     : "output",
        "max_batch_size"  : 0,
        "instance_count"  : 0,
        "instance_kind"   : "KIND_AUTO",
        "opset_version"   : 17
    }""")
    settings = settings.Clone()
    settings.ValidateAndAssignDefaults(defaults)

    model_format = settings["format"].GetString()
    if model_format not in _PLATFORMS:
        raise ValueError(
            f"Unsupported format \"{model_format}\". Use one of {tuple(_PLATFORMS)}.")
    version = settings["model_version"].GetInt()
    if version < 1:
        raise ValueError(f"\"model_version\" must be >= 1 [ model_version = {version} ].")

    model_name = settings["model_name"].GetString()
    repository = Path(settings["repository_path"].GetString())
    version_directory = repository / model_name / str(version)
    version_directory.mkdir(parents=True, exist_ok=True)
    model_file = version_directory / _MODEL_FILES[model_format]

    if isinstance(sample_inputs, (tuple, list)):
        if len(sample_inputs) != 1:
            raise ValueError(
                f"Triton export expects a single input tensor, got {len(sample_inputs)}.")
        sample = sample_inputs[0]
    else:
        sample = sample_inputs
    sample = torch.as_tensor(sample)

    model = model.eval()
    with torch.no_grad():
        reference = model(sample)

    input_name = settings["input_name"].GetString()
    output_name = settings["output_name"].GetString()
    card = None
    if settings["card_file"].GetString():
        card = model_registry.LoadModelCard(settings["card_file"].GetString())
        input_fields = _CardWidths(card, "input_fields")
        output_fields = _CardWidths(card, "output_fields")
        if input_fields:
            input_name = "__".join(input_fields)
        if output_fields:
            output_name = "__".join(output_fields)

    if model_format == "onnx":
        # dynamic entity axis: the served graph must accept any mesh size
        torch.onnx.export(
            model, (sample,), str(model_file),
            input_names=[input_name], output_names=[output_name],
            dynamic_axes={input_name: {0: "n_entities"},
                          output_name: {0: "n_entities"}},
            opset_version=settings["opset_version"].GetInt())
    else:
        torch.jit.script(model).save(str(model_file))

    def triton_dtype(tensor):
        name = str(tensor.dtype).replace("torch.", "")
        if name not in _TRITON_DTYPES:
            raise ValueError(
                f"Unsupported tensor dtype \"{name}\" for Triton; "
                f"use one of {tuple(_TRITON_DTYPES)}.")
        return _TRITON_DTYPES[name]

    inputs = [{"name": input_name, "data_type": triton_dtype(sample),
               "dims": [-1, int(sample.shape[-1])]}]
    outputs = [{"name": output_name, "data_type": triton_dtype(reference),
                "dims": [-1, int(reference.shape[-1])]}]

    instance_group = None
    if settings["instance_count"].GetInt() > 0:
        instance_group = [{"count": settings["instance_count"].GetInt(),
                           "kind": settings["instance_kind"].GetString()}]

    config_text = MakeTritonConfig(
        model_name, _PLATFORMS[model_format], inputs, outputs,
        max_batch_size=settings["max_batch_size"].GetInt(),
        instance_group=instance_group)
    config_file = repository / model_name / "config.pbtxt"
    config_file.write_text(config_text)

    if card is not None:
        model_registry.SaveModelCard(str(model_file), card)
    return str(config_file)
