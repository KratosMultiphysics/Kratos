"""Fine-tuning a pretrained DoMINO on Kratos data.

Two recipes, both adapting a frozen pretrained checkpoint rather than
training one from scratch.

**Predictor-corrector.** NVIDIA's own recipe, defined as
`Y_finetuned = Y_predictor + Y_corrector`: the pretrained checkpoint is the
frozen predictor, and a trainable network learns its error. Note what that
is upstream and what it is here. Upstream's corrector is *a second full
DoMINO* (~10 M parameters) trained on `ground_truth - base_prediction`; its
lightness is in how fast it converges, not in its size, and at full mesh
resolution it needs far more memory than a single consumer GPU has. What is
shipped here is the same decomposition with a small residual head, and -
exactly as upstream does in its first two stages - the predictor's output is
computed once and **cached**, so it never runs inside the training loop.
That makes the cost independent of the predictor's size.

**LoRA.** Low-rank adapters on the pretrained weights themselves
(`physicsnemo.experimental.peft`). Roughly 1-2 % of the parameters are
trainable, and `MergeAndSave` folds them back into an ordinary `.mdlus`
that `model_registry` loads and `DominoInferenceProcess` deploys with no
change to the **model** settings. Note it still needs the same
de-normalization block as the checkpoint it was adapted from: a merged
model lives in the pretrained model's normalized output space, so
`scaling_factors_file`/`normalization`/`redimensionalize` remain required.
For the same reason `CacheBasePredictions` returns raw normalized output -
form residuals against ground truth in one space consistently, not against
physical Kratos values. This is usually the better option; the
predictor-corrector path exists because it is the recipe the literature and
NVIDIA's documentation describe.

Accuracy claims for either belong to NVIDIA, who describe their own
fine-tuning results as preliminary and report them on 18 training samples.
Nothing here reproduces or endorses a number.

torch and physicsnemo are imported lazily.
"""

import numpy

import KratosMultiphysics as Kratos


def _TryImportTorch():
    try:
        import torch
        return torch
    except ImportError as e:
        raise ImportError(
            "PhysicsNeMoApplication.domino_finetune requires torch, which could not be "
            "imported. Install it with e.g. 'pip install torch'.") from e


def _TryImportPeft():
    try:
        from physicsnemo.experimental.peft import (
            LoRAConfig, apply_lora, merge_lora, save_adapter, load_adapter)
        return LoRAConfig, apply_lora, merge_lora, save_adapter, load_adapter
    except ImportError as e:
        raise ImportError(
            "PhysicsNeMoApplication.domino_finetune requires physicsnemo's experimental "
            "PEFT module for the LoRA recipe, which could not be imported. Install it "
            "with e.g. 'pip install -U nvidia-physicsnemo'.") from e


def CacheBasePredictions(model, batches, device=None):
    """Runs the frozen predictor once per case and returns its outputs.

    This is the step that makes the recipe affordable: the predictor never
    runs inside the training loop, so the corrector's cost does not depend
    on the predictor's size.

    Args:
        model: The frozen pretrained DoMINO.
        batches: Iterable of already-preprocessed sample dicts (what
            DominoInferenceProcess feeds the model).
        device: Optional torch device.

    Returns:
        list of (volume, surface) numpy arrays, either possibly None,
        matching what DoMINO.forward returns.
    """
    torch = _TryImportTorch()
    was_training = model.training
    model.eval()
    cached = []
    try:
        with torch.no_grad():
            for batch in batches:
                if device is not None:
                    batch = {key: value.to(device) if hasattr(value, "to") else value
                             for key, value in batch.items()}
                volume, surface = model(batch)
                cached.append((
                    None if volume is None else volume.detach().cpu().numpy(),
                    None if surface is None else surface.detach().cpu().numpy()))
    finally:
        model.train(was_training)
    return cached


def CreateCorrector(n_features: int, n_outputs: int, hidden: int = 64,
                    n_layers: int = 3, dtype=None):
    """A small residual head predicting the predictor's error.

    Deliberately tiny: it is applied per surface entity, and its whole point
    is to be cheap next to a frozen 10 M-parameter predictor.
    """
    torch = _TryImportTorch()
    if n_layers < 1:
        raise ValueError("n_layers must be at least 1.")

    layers = []
    width = n_features
    for _ in range(n_layers - 1):
        layers += [torch.nn.Linear(width, hidden), torch.nn.GELU()]
        width = hidden
    layers.append(torch.nn.Linear(width, n_outputs))
    corrector = torch.nn.Sequential(*layers)
    # the last layer starts at zero so the untrained corrector is exactly the
    # identity on the predictor - training can only improve on it
    torch.nn.init.zeros_(corrector[-1].weight)
    torch.nn.init.zeros_(corrector[-1].bias)
    return corrector.to(dtype) if dtype is not None else corrector


def _CorrectorDtype(corrector, torch):
    """The dtype the corrector's own weights are in.

    `CreateCorrector` takes a `dtype`, so the inputs have to follow it rather
    than being pinned to float32 - a float64 corrector otherwise fails inside
    the first linear layer with an opaque "mat1 and mat2 must have the same
    dtype".
    """
    for parameter in corrector.parameters():
        return parameter.dtype
    return torch.float32


def TrainCorrector(corrector, features, residuals, epochs: int = 100,
                   learning_rate: float = 1e-3, batch_size: int = 0,
                   device=None, echo_interval: int = 0):
    """Fits the corrector on (features -> residual) pairs.

    `residuals` are `ground_truth - base_prediction`; nothing here recomputes
    the base prediction, which is the point of caching it.

    Returns:
        list of per-epoch mean losses.
    """
    torch = _TryImportTorch()
    dtype = _CorrectorDtype(corrector, torch)
    features = torch.as_tensor(numpy.asarray(features), dtype=dtype)
    residuals = torch.as_tensor(numpy.asarray(residuals), dtype=dtype)
    if features.shape[0] != residuals.shape[0]:
        raise ValueError(
            f"features has {features.shape[0]} rows but residuals has "
            f"{residuals.shape[0]}; they must line up entity for entity.")
    if device is not None:
        corrector, features, residuals = (corrector.to(device), features.to(device),
                                          residuals.to(device))

    optimizer = torch.optim.Adam(corrector.parameters(), lr=learning_rate)
    n_rows = features.shape[0]
    step = batch_size if batch_size > 0 else n_rows
    history = []
    corrector.train()
    for epoch in range(epochs):
        permutation = torch.randperm(n_rows, device=features.device)
        total = 0.0
        for start in range(0, n_rows, step):
            rows = permutation[start:start + step]
            optimizer.zero_grad()
            loss = torch.nn.functional.mse_loss(corrector(features[rows]), residuals[rows])
            loss.backward()
            optimizer.step()
            total += loss.item() * len(rows)
        history.append(total / n_rows)
        if echo_interval and (epoch + 1) % echo_interval == 0:
            Kratos.Logger.PrintInfo(
                "PhysicsNeMoApplication.domino_finetune",
                f"corrector epoch {epoch + 1}/{epochs}: loss {history[-1]:.6e}")
    return history


def ApplyCorrector(corrector, base_prediction, features, device=None):
    """`Y_predictor + Y_corrector`, the combination step."""
    torch = _TryImportTorch()
    features = torch.as_tensor(numpy.asarray(features),
                               dtype=_CorrectorDtype(corrector, torch))
    if device is not None:
        corrector, features = corrector.to(device), features.to(device)
    corrector.eval()
    with torch.no_grad():
        correction = corrector(features).cpu().numpy()
    base = numpy.asarray(base_prediction)
    if correction.shape != base.shape:
        raise ValueError(
            f"corrector returned {list(correction.shape)} but the base prediction is "
            f"{list(base.shape)}.")
    return base + correction


def ApplyLora(model, rank: int = 8, target_pattern: str = r"^agg_model_surf\..*",
              alpha=None, dropout: float = 0.0):
    """Wraps a pretrained DoMINO's linear layers in LoRA adapters.

    Args:
        target_pattern: A regex over **fully-qualified module names**, and
            only `nn.Linear`-like layers are wrapped. This is easy to get
            wrong: `^solution_calculator_surf\\..*` matches **zero** layers,
            because SolutionCalculatorSurface re-uses `nn_basis_surf` and
            `agg_model_surf` by reference rather than owning submodules. The
            addressable groups are `geo_rep_volume`, `geo_rep_surface`,
            `nn_basis_surf`, `fc_p_surf`, `surface_local_geo_encodings`,
            `volume_local_geo_encodings` and `agg_model_surf`.

    Returns:
        (model, n_wrapped, n_trainable). Raises if the pattern matched
        nothing, since silently training zero parameters looks exactly like
        a successful run.
    """
    LoRAConfig, apply_lora, _, _, _ = _TryImportPeft()
    configuration = LoRAConfig(rank=rank, target_pattern=target_pattern,
                               alpha=alpha, lora_dropout=dropout)
    result = apply_lora(model, configuration)

    # upstream already raises when a pattern matches nothing ("apply_lora
    # matched 0 wrappable layers"), which is the behaviour we want; this is a
    # backstop in case that guard ever moves, because silently training zero
    # parameters looks exactly like a successful run.
    n_wrapped = int(getattr(result, "n_wrapped", 0))
    if not n_wrapped:
        raise ValueError(
            f"target_pattern {target_pattern!r} matched no wrappable layer. It is "
            "matched against fully-qualified module names, and only Linear layers are "
            "wrapped; see this function's docstring for the addressable groups.")
    return model, n_wrapped, int(getattr(result, "n_trainable", 0))


def MergeAndSave(model, path: str):
    """Folds the adapters back in and writes a plain `.mdlus`.

    The merged file is an ordinary physicsnemo checkpoint: `model_registry`
    loads it with `"checkpoint_type": "physicsnemo"` and
    `DominoInferenceProcess` deploys it unchanged. That is what makes LoRA
    the least invasive of the two recipes.
    """
    _, _, merge_lora, _, _ = _TryImportPeft()
    merged = merge_lora(model)
    merged.save(str(path))
    return merged
