"""Derivative-informed (Sobolev) training against Kratos's adjoint gradients.

A surrogate fitted on field values alone is graded on values alone, and its
*derivatives* are whatever the fit happened to leave behind. That matters as
soon as the surrogate is used for anything gradient-driven - shape
optimization, sensitivity screening, an inverse problem - because those read
dJ/dX, a quantity the training objective never looked at.

Kratos already computes that quantity exactly. ``AdjointSensitivityProcess``
writes it into a nodal variable, an export process carries it into the
samples like any other field, and ``CreateNpzDataset`` puts it in the
targets. What is left is a loss term that reads it, which is this module:

    dJ_model/dX  vs  dJ_kratos/dX

evaluated by differentiating the *model's own* objective with respect to the
coordinate channels of its input. The idiom is the one
``physics_informed.MakePhysicsLossTerm``'s ``"autodiff"`` path established:
the coordinate channels are detached into a fresh leaf, the input is
rebuilt around them, and the model is re-run - a graph that only *reaches*
the coordinates is not enough, they have to be the leaf. ``create_graph=True``
keeps the second derivative, so the matching term is itself trainable.

Layout: inputs ``(..., N, C_in)`` with the first ``coordinate_channels``
columns the coordinates, targets ``(..., N, C_target)`` with
``gradient_columns`` holding dJ/dX. Both a single ``(N, C)`` mesh and a
batched ``(B, N, C)`` stack work; the objective is summed over everything but
never mixes samples, because the model processes them independently, so one
backward pass yields each sample's own gradient.

torch is imported lazily; module import stays ML-free.
"""

import KratosMultiphysics as Kratos


def _TryImportTorch():
    try:
        import torch
        return torch
    except ImportError as e:
        raise ImportError(
            "PhysicsNeMoApplication.sobolev_training requires torch, which could not be "
            "imported. Install it with e.g. 'pip install torch'.") from e


_REDUCTIONS = ("mse", "relative")


def _ReadSettings(settings: Kratos.Parameters):
    defaults = Kratos.Parameters("""{
        "coordinate_channels" : 3,
        "objective_channels"  : [],
        "objective_weights"   : [],
        "gradient_columns"    : [],
        "weight"              : 1.0,
        "reduction"           : "mse"
    }""")
    settings = settings.Clone()
    settings.ValidateAndAssignDefaults(defaults)

    coordinate_channels = settings["coordinate_channels"].GetInt()
    if coordinate_channels < 1:
        raise ValueError(
            f"\"coordinate_channels\" must be >= 1, got {coordinate_channels}.")

    objective_channels = [int(v) for v in settings["objective_channels"].GetVector()]
    objective_weights = [float(v) for v in settings["objective_weights"].GetVector()]
    if objective_weights and len(objective_weights) != len(objective_channels):
        raise ValueError(
            f"\"objective_weights\" has {len(objective_weights)} entries but "
            f"\"objective_channels\" has {len(objective_channels)}.")

    gradient_columns = [int(v) for v in settings["gradient_columns"].GetVector()]
    if gradient_columns and len(gradient_columns) != coordinate_channels:
        raise ValueError(
            f"\"gradient_columns\" names {len(gradient_columns)} column(s) but the gradient "
            f"has {coordinate_channels} component(s) - one per coordinate channel.")

    reduction = settings["reduction"].GetString()
    if reduction not in _REDUCTIONS:
        raise ValueError(
            f"Unsupported \"reduction\" \"{reduction}\". Use {' or '.join(_REDUCTIONS)}.")

    return (coordinate_channels, objective_channels, objective_weights, gradient_columns,
            settings["weight"].GetDouble(), reduction)


def SensitivityGradient(model, inputs, coordinate_channels: int = 3,
                        objective_channels=(), objective_weights=(),
                        create_graph: bool = False):
    """d(model objective)/d(coordinates), the surrogate's own dJ/dX.

    The quantity the loss term matches, exposed on its own because measuring
    a surrogate's gradient accuracy is the point of training this way: run it
    at a held-out design and compare against the FEM adjoint.

    Args:
        model: A torch module mapping ``(..., N, C_in)`` to ``(..., N, C_out)``.
        inputs: The input tensor; its first ``coordinate_channels`` columns
            are the coordinates.
        objective_channels: Prediction channels forming J; empty means all.
        objective_weights: Per-channel weights; empty means ones.
        create_graph: Keep the graph, so the result is differentiable w.r.t.
            the model's weights (what the loss term needs).

    Returns:
        (..., N, coordinate_channels) dJ/dX, with the same batch shape as
        the input.
    """
    torch = _TryImportTorch()

    # the coordinates must be the LEAF the model input derives from - a graph
    # that merely passes through them does not give torch anything to
    # differentiate against
    coordinates = inputs[..., :coordinate_channels].detach().clone()
    coordinates.requires_grad_(True)
    rebuilt = torch.cat([coordinates, inputs[..., coordinate_channels:]], dim=-1)

    prediction = model(rebuilt)
    if objective_channels:
        selected = prediction[..., list(objective_channels)]
    else:
        selected = prediction
    if objective_weights:
        weights = torch.as_tensor(list(objective_weights), dtype=selected.dtype,
                                  device=selected.device)
        selected = selected * weights
    # Samples never mix (the model maps each independently), so one backward
    # over the summed objective yields every sample's own gradient.
    objective = selected.sum()

    gradient, = torch.autograd.grad(objective, coordinates, create_graph=create_graph)
    return gradient


def MakeSensitivityLossTerm(settings: Kratos.Parameters):
    """Builds a derivative-matching loss term for ``training_utils.TrainModel``.

    The returned callable takes FOUR arguments -
    ``term(model, inputs, prediction, targets)`` - which is how TrainModel
    knows to hand it the batch targets; three-argument terms are unaffected.

    Settings:
        coordinate_channels (3)
            How many leading input columns are the coordinates.
        objective_channels ([])
            Prediction channels forming the model's objective J; empty means
            every channel. This must be the same objective the stored
            gradient was computed for - use the channels matching the
            ``adjoint_bridge.MakeObjectiveWeights`` block that produced it.
        objective_weights ([])
            Per-channel weights for those channels; empty means ones.
        gradient_columns ([])
            Target columns holding dJ/dX, one per coordinate channel. Empty
            means the LAST ``coordinate_channels`` columns of the target,
            which is where appending the sensitivity field to
            ``CreateNpzDataset``'s ``output_keys`` puts it.
        weight (1.0)
            Multiplies the term before it is added to the data loss.
        reduction ("mse")
            ``"mse"`` - mean squared error on the gradient. ``"relative"`` -
            that error divided by the reference gradient's own energy, which
            is what makes the term scale-free when dJ/dX is orders of
            magnitude smaller or larger than the field itself (a shape
            gradient usually is).

    Returns:
        term(model, inputs, prediction, targets) -> scalar tensor.
    """
    torch = _TryImportTorch()
    (coordinate_channels, objective_channels, objective_weights, gradient_columns,
     weight, reduction) = _ReadSettings(settings)

    def SensitivityLossTerm(model, inputs, prediction, targets):
        if inputs.shape[-1] < coordinate_channels:
            raise ValueError(
                f"\"coordinate_channels\" is {coordinate_channels} but the inputs have "
                f"only {inputs.shape[-1]} channel(s).")
        columns = gradient_columns or list(
            range(targets.shape[-1] - coordinate_channels, targets.shape[-1]))
        if min(columns) < 0 or max(columns) >= targets.shape[-1]:
            raise ValueError(
                f"\"gradient_columns\" {columns} do not fit targets with "
                f"{targets.shape[-1]} column(s).")

        gradient = SensitivityGradient(
            model, inputs, coordinate_channels, objective_channels, objective_weights,
            create_graph=True)
        reference = targets[..., columns].to(gradient.dtype)

        difference = gradient - reference
        if reduction == "mse":
            value = difference.square().mean()
        else:
            # a zero reference gradient would divide by zero; the fallback is
            # the plain mean square, which is the right limit
            energy = reference.square().sum()
            value = difference.square().sum() / torch.clamp(energy, min=1e-30)
        return weight * value

    return SensitivityLossTerm
