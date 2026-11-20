"""Model checkpoint loading and device placement.

Single source of truth for how this application loads trained models, so
every process and strategy resolves checkpoints and devices identically.

torch and physicsnemo are optional runtime dependencies, imported lazily
inside LoadModel only.
"""

import numpy

import KratosMultiphysics as Kratos


def _TryImportTorch():
    try:
        import torch
        return torch
    except ImportError as e:
        raise ImportError(
            "PhysicsNeMoApplication.model_registry requires torch, which could not be "
            "imported. Install it with e.g. 'pip install torch'.") from e


def _TryImportPhysicsNemo():
    try:
        import physicsnemo
        return physicsnemo
    except ImportError as e:
        raise ImportError(
            "Loading a \"physicsnemo\" checkpoint requires physicsnemo, which could not be "
            "imported. Install it with e.g. 'pip install nvidia-physicsnemo'.") from e


def ResolveDevice(device_name: str):
    """Resolves a device string ("auto"/"cpu"/"cuda"[:i]) to a torch.device."""
    torch = _TryImportTorch()
    if device_name == "auto":
        device_name = "cuda" if torch.cuda.is_available() else "cpu"
    return torch.device(device_name)


def LoadModel(settings: Kratos.Parameters):
    """Loads a trained model from a checkpoint and places it on a device.

    Settings:
        checkpoint_file: Path to the checkpoint.
        checkpoint_type: "torchscript" (torch.jit.load, default) or
            "physicsnemo" (physicsnemo.Module.from_checkpoint, .mdlus files).
        device: "auto" (default; cuda when available, else cpu), "cpu" or
            "cuda"[:index].
        model_card_policy: "advisory" (default), "strict" or "ignore" - how
            deployment processes treat model-card mismatches (see
            LoadModelWithCardCheck; LoadModel itself never reads the card).
        torch_compile: false (default). When true, the loaded model is
            wrapped with torch.compile(fullgraph=True). Only "physicsnemo"
            checkpoints are compilable - TorchScript modules do not compose
            with dynamo, so "torchscript" + torch_compile raises.
        nvtx_ranges: false (default). When true, enables the NVTX ranges
            around the deployment hot paths (see utilities.nvtx_utils);
            they emit only when CUDA is available.

    Returns:
        (model, device): The model in eval mode on the device, and the
        torch.device it was placed on.
    """
    default_settings = Kratos.Parameters("""{
        "checkpoint_file"   : "PLEASE_SPECIFY_CHECKPOINT_FILE",
        "checkpoint_type"   : "torchscript",
        "device"            : "auto",
        "model_card_policy" : "advisory",
        "torch_compile"     : false,
        "nvtx_ranges"       : false
    }""")
    settings.ValidateAndAssignDefaults(default_settings)

    torch = _TryImportTorch()
    device = ResolveDevice(settings["device"].GetString())

    if settings["nvtx_ranges"].GetBool():
        from KratosMultiphysics.PhysicsNeMoApplication.utilities import nvtx_utils
        nvtx_utils.EnableNvtxRanges()

    checkpoint_file = settings["checkpoint_file"].GetString()
    checkpoint_type = settings["checkpoint_type"].GetString()
    if checkpoint_type == "torchscript":
        model = torch.jit.load(checkpoint_file, map_location=device)
    elif checkpoint_type == "physicsnemo":
        physicsnemo = _TryImportPhysicsNemo()
        model = physicsnemo.Module.from_checkpoint(checkpoint_file)
    else:
        raise ValueError(
            f"Unsupported checkpoint type \"{checkpoint_type}\". "
            "Use \"torchscript\" or \"physicsnemo\".")

    model = model.eval().to(device)
    if settings["torch_compile"].GetBool():
        if checkpoint_type == "torchscript":
            raise ValueError(
                "\"torch_compile\" is not supported for \"torchscript\" checkpoints: "
                "TorchScript modules do not compose with torch.compile. Save the model "
                "as a physicsnemo checkpoint (.mdlus) to compile it.")
        model = torch.compile(model, fullgraph=True)
    return model, device


def _CardPath(checkpoint_file):
    import pathlib
    return pathlib.Path(str(checkpoint_file) + ".card.json")


def SaveModelCard(checkpoint_file, card: dict) -> None:
    """Writes a model card sidecar ("<checkpoint_file>.card.json").

    The card is a free-form JSON dict describing what the checkpoint was
    trained for. Recommended keys, checked by the deployment processes:
    "input_fields" / "output_fields" as [{"variable_name", "data_location"}],
    and "output_normalization" (see LoadOutputNormalization) when the model
    was trained on normalized targets. Anything else ("grid_shape",
    "history_size", training provenance, ...) travels along untouched.

    This function never synthesizes keys - the card written is exactly the
    dict given. Anything auto-derived belongs in whatever builds that dict.
    """
    import json
    with open(_CardPath(checkpoint_file), "w") as f:
        json.dump(card, f, indent=4)


def LoadModelCard(checkpoint_file):
    """Returns the model card dict for a checkpoint, or None if there is none."""
    import json
    path = _CardPath(checkpoint_file)
    if not path.is_file():
        return None
    with open(path) as f:
        return json.load(f)


def _FieldList(card_entry):
    return [(entry.get("variable_name"), entry.get("data_location")) for entry in card_entry]


_MODEL_CARD_POLICIES = ("advisory", "strict", "ignore")


def ValidateFieldsAgainstCard(card, input_specs, output_specs, tag: str,
                              policy: str = "advisory") -> bool:
    """Compares configured field specs with a model card's.

    Policies:
        "advisory" (default): a mismatch produces one detailed
            KRATOS_WARNING (the card may be stale, or the user may
            deliberately remap fields) and never raises.
        "strict": a mismatch raises RuntimeError with the same message.
        "ignore": mismatches are silent.

    Returns True when everything matches (or the card doesn't constrain the
    fields), False when a mismatch was found (and the policy let it pass).

    Args:
        card: The dict from LoadModelCard (None is accepted and passes).
        input_specs/output_specs: [(variable_name, data_location)] as
            configured in the deployment process.
        tag: Logger tag of the calling process.
        policy: "advisory", "strict" or "ignore".
    """
    if policy not in _MODEL_CARD_POLICIES:
        raise ValueError(
            f"Unsupported model card policy \"{policy}\". Use one of {_MODEL_CARD_POLICIES}.")
    if card is None:
        return True
    ok = True
    for key, configured in (("input_fields", list(input_specs)), ("output_fields", list(output_specs))):
        if key not in card:
            continue
        expected = _FieldList(card[key])
        if expected != configured:
            message = (
                f"Configured {key} {configured} do not match the model card's {expected}. "
                "The model may have been trained for different fields")
            if policy == "strict":
                raise RuntimeError(
                    f"{tag}: {message}; the model card policy is \"strict\", so execution stops.")
            if policy == "advisory":
                Kratos.Logger.PrintWarning(
                    tag, f"{message}; the card is advisory, so execution continues.")
            ok = False
    return ok


_NORMALIZATION_TYPES = ("none", "mean_std", "min_max")


def LoadOutputNormalization(model_settings: Kratos.Parameters, checkpoint_file=None):
    """The card's "output_normalization" entry, or None when there is none.

    A model trained on normalized targets emits normalized predictions, and
    writing those onto Kratos variables as if they were physical is wrong
    by whatever the scaling was. The card is where that scaling travels
    with the checkpoint.

    Schema:
        {"type": "mean_std", "mean": [...], "std":  [...]}
        {"type": "min_max",  "min":  [...], "max":  [...],
                             "range": [0.0, 1.0]}     // optional
        {"type": "none"}

    Each array is length 1 (broadcast over all channels) or length
    total_width (per channel, in the concatenated output_fields order).
    "range" is the interval the training normalization mapped onto -
    [0, 1] by default, though DoMINO's convention is [-1, 1].

    Deliberately read regardless of "model_card_policy": "ignore" means
    "do not validate the field lists", and silently dropping the
    de-normalization would reintroduce exactly the bug this exists to
    prevent. Configurations written before this key existed are unaffected,
    since a card without it yields None.

    Returns:
        The entry dict, or None (the identity path).
    """
    return _LoadNormalizationEntry(model_settings, "output_normalization", checkpoint_file)


def LoadInputNormalization(model_settings: Kratos.Parameters, checkpoint_file=None):
    """The card's "input_normalization" entry, or None when there is none.

    The symmetric half of LoadOutputNormalization: a model trained on
    standardized FEATURES expects standardized features, and feeding it
    the raw Kratos fields is wrong by the same silent factor - measured as
    an 18% position drift on the particle path, where
    CreateParticleTrajectoryDataset(normalize=True) standardizes the
    velocity history the deployment process then fed raw. Same schema,
    same broadcast rule (length 1 or the concatenated input_fields width),
    same "read regardless of the policy" rule. The rule every process
    follows: inputs are normalized BEFORE the OOD-guard check, since the
    guard was calibrated on what the model saw in training.
    """
    return _LoadNormalizationEntry(model_settings, "input_normalization", checkpoint_file)


def MakeMeanStdNormalization(mean, std) -> dict:
    """A "mean_std" card entry from per-channel statistics (array-likes)."""
    mean = numpy.asarray(mean, dtype=numpy.float64).reshape(-1)
    std = numpy.asarray(std, dtype=numpy.float64).reshape(-1)
    if mean.size != std.size:
        raise ValueError(
            f"mean has {mean.size} entries but std has {std.size}; they describe the "
            "same channels.")
    return {"type": "mean_std", "mean": mean.tolist(), "std": std.tolist()}


def _LoadNormalizationEntry(model_settings: Kratos.Parameters, key: str, checkpoint_file=None):
    if checkpoint_file is None:
        # Not every process keys its card off "checkpoint_file": ONNX uses
        # "onnx_file", Triton "card_file", and an ensemble names its members
        # in "checkpoint_files". Those pass their own path; a settings block
        # with none of them simply has no card.
        if model_settings.Has("checkpoint_file"):
            checkpoint_file = model_settings["checkpoint_file"].GetString()
        elif model_settings.Has("checkpoint_files"):
            members = model_settings["checkpoint_files"].GetStringArray()
            checkpoint_file = members[0] if members else None
        if not checkpoint_file:
            return None

    card = LoadModelCard(checkpoint_file)
    if not card:
        return None
    normalization = card.get(key)
    if not normalization:
        return None

    kind = normalization.get("type", "none")
    if kind not in _NORMALIZATION_TYPES:
        raise ValueError(
            f"Unsupported \"{key}\" type \"{kind}\" in the model card. "
            f"Supported: {', '.join(_NORMALIZATION_TYPES)}.")
    if kind == "none":
        return None

    required = ("mean", "std") if kind == "mean_std" else ("min", "max")
    for name in required:
        if name not in normalization:
            raise ValueError(
                f"\"{key}\" of type \"{kind}\" needs \"{name}\".")
    return normalization


def _IsTorchTensor(value) -> bool:
    # duck-typed so the numpy path never imports torch
    return hasattr(value, "detach") and hasattr(value, "cpu")


def _NormalizationScaleOffset(normalization, n_channels: int,
                              key: str = "output_normalization"):
    """The validated (scale, offset) of a card entry: physical =
    normalized * scale + offset, each a float64 vector of length 1
    (broadcast) or n_channels."""
    def _Vector(name):
        vector = numpy.asarray(normalization[name], dtype=numpy.float64).reshape(-1)
        if vector.size not in (1, n_channels):
            raise ValueError(
                f"\"{key}\" has {vector.size} entries for \"{name}\" but "
                f"the prediction has {n_channels} channels; the card does not belong to "
                "this model.")
        return vector

    if normalization["type"] == "mean_std":
        scale, offset = _Vector("std"), _Vector("mean")
    else:
        low, high = _Vector("min"), _Vector("max")
        interval = normalization.get("range", [0.0, 1.0])
        span = float(interval[1]) - float(interval[0])
        if span == 0.0:
            raise ValueError(f"\"{key}\" range must not be degenerate.")
        scale = (high - low) / span
        offset = low - float(interval[0]) * scale

    if numpy.any(scale == 0.0):
        raise ValueError(
            f"\"{key}\" has a zero scale, which would make the inverse "
            "undefined (a constant channel needs a scale of 1, not 0).")
    return scale, offset


def _BroadcastAlongAxis(vector, ndim: int, channel_axis: int):
    """Reshapes a length-1-or-C vector to broadcast along channel_axis."""
    if ndim < 2:
        return vector
    shape = [1] * ndim
    shape[channel_axis % ndim] = vector.size
    return vector.reshape(shape)


def ApplyOutputNormalization(prediction, normalization, scale_only: bool = False,
                             channel_axis: int = -1):
    """Inverts a training normalization on a model's prediction.

    Args:
        prediction: An array-like with the channels on channel_axis -
            (n_entities, total_width) for the row-ordered writers, or a
            channels-first (C, *spatial) grid with channel_axis=0. Returned
            unchanged - the same object - when normalization is None, so
            the identity path costs nothing and preserves dtype.
        normalization: A LoadOutputNormalization entry.
        scale_only: Apply the scale but NOT the offset. Required for a
            standard deviation or any other spread: shifting a spread by
            the mean is meaningless, and doing so is the mistake a single
            shared hook most easily makes.
        channel_axis: The axis the card's per-channel vectors run along.

    Returns:
        The de-normalized prediction, with the same type as the input - a
        torch tensor in, a torch tensor out, on the SAME device and dtype
        and with its autograd graph intact: the arithmetic runs in torch
        for a tensor (it used to bounce through numpy, which returned a
        CUDA prediction on the host and cut every gradient - and the
        surrogate response function differentiates through this). numpy
        in, float64 numpy out. Returns the untouched input object when
        there is nothing to do.
    """
    if normalization is None:
        return prediction

    is_tensor = _IsTorchTensor(prediction)
    if not is_tensor:
        prediction = numpy.asarray(prediction, dtype=numpy.float64)
    ndim = int(prediction.ndim)
    n_channels = int(prediction.shape[channel_axis]) if ndim > 1 else 1
    scale, offset = _NormalizationScaleOffset(normalization, n_channels)
    scale = _BroadcastAlongAxis(scale, ndim, channel_axis)
    offset = _BroadcastAlongAxis(offset, ndim, channel_axis)

    if is_tensor:
        torch = _TryImportTorch()
        result = prediction * torch.as_tensor(
            scale, dtype=prediction.dtype, device=prediction.device)
        if not scale_only:
            result = result + torch.as_tensor(
                offset, dtype=prediction.dtype, device=prediction.device)
        return result

    values = prediction * scale
    if not scale_only:
        values = values + offset
    return values


def ApplyInputNormalization(features, normalization, channel_axis: int = -1):
    """Applies a training normalization to the features a model is fed.

    The forward map - (x - mean) / std, or (x - min) / (max - min) onto
    "range" - i.e. the exact inverse of ApplyOutputNormalization for the
    same entry, so the two round-trip to the identity. Same type
    preservation (torch stays torch, on its device and dtype; numpy gives
    float64) and the same identity-object path for None.

    Args:
        features: An array-like with the channels on channel_axis - the
            concatenated (n_entities, total_input_width) features, or a
            channels-first (C, *spatial) grid with channel_axis=0.
        normalization: A LoadInputNormalization entry.
        channel_axis: The axis the card's per-channel vectors run along.
    """
    if normalization is None:
        return features

    is_tensor = _IsTorchTensor(features)
    if not is_tensor:
        features = numpy.asarray(features, dtype=numpy.float64)
    ndim = int(features.ndim)
    n_channels = int(features.shape[channel_axis]) if ndim > 1 else 1
    scale, offset = _NormalizationScaleOffset(normalization, n_channels, key="input_normalization")
    scale = _BroadcastAlongAxis(scale, ndim, channel_axis)
    offset = _BroadcastAlongAxis(offset, ndim, channel_axis)

    if is_tensor:
        torch = _TryImportTorch()
        return (features - torch.as_tensor(offset, dtype=features.dtype, device=features.device)) \
            / torch.as_tensor(scale, dtype=features.dtype, device=features.device)
    return (features - offset) / scale


def LoadModelWithCardCheck(model_settings: Kratos.Parameters, input_specs, output_specs, tag: str):
    """LoadModel plus model-card validation - the single entry point every
    deployment process uses at (lazy) model-load time.

    The policy comes from model_settings["model_card_policy"] ("advisory"
    by default, "strict" to refuse mismatched deployments, "ignore" to
    silence the check).

    Returns:
        (model, device), as LoadModel.
    """
    model, device = LoadModel(model_settings)
    policy = model_settings["model_card_policy"].GetString()
    if policy != "ignore":
        card = LoadModelCard(model_settings["checkpoint_file"].GetString())
        ValidateFieldsAgainstCard(card, input_specs, output_specs, tag, policy)
    return model, device
