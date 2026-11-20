"""Virtual Foundry GraphNet: sintering and deformation surrogates.

VFGN (physicsnemo.models.vfgn) is a Learning-to-Simulate model for metal
binder-jetting sintering: particles shrink and deform under a temperature
schedule. It shares the shipped Lagrangian particle machinery - the same
radius graph, the same history window, the same semi-implicit integration -
but differs in three ways that this module absorbs:

- it consumes a POSITION SEQUENCE (N, T, 3), not a velocity window (N, K*3);
- it takes `senders`/`receivers` as separate 1-D tensors, not an edge_index
  or a PyG graph (its own `graph_mode` argument is dead code - the model
  never builds a graph, the caller always supplies one);
- it REQUIRES normalization statistics, and silently produces NaNs without
  usable ones.

**`VFGNLearnedSimulator.forward()` is unusable in physicsnemo 2.2.** Its
shape guard demands a 2-D `next_positions` while the body's arithmetic needs
`(N, predict_length, 3)`; every input either raises or returns a
shape-mismatched pair. Verified across the matrix, and pinned by a test so we
learn when it is fixed. Rollout therefore goes through the public
`inference()`, which works, and training through the same encode/decode
composition `forward()` would have performed.

Two more upstream contracts encoded here: `num_dimensions` must equal
`3 * predict_length` (the decoder's output width and the reshape must agree),
and the model gains lazily-created parameters on its first forward - so an
optimizer must be built AFTER a warm-up pass or it will miss half the model.

torch/physicsnemo/torch_scatter are optional runtime dependencies.
"""

import numpy

import KratosMultiphysics as Kratos

_STATS_KEYS = ("velocity", "acceleration", "context")


def _TryImportTorch():
    try:
        import torch
        return torch
    except ImportError as e:
        raise ImportError(
            "PhysicsNeMoApplication.vfgn_bridge requires torch, which could not be "
            "imported. Install it with e.g. 'pip install torch'.") from e


def _TryImportVfgn():
    try:
        import torch_scatter  # noqa: F401  (VFGN's message passing needs it)
    except ImportError as e:
        raise ImportError(
            "PhysicsNeMoApplication.vfgn_bridge requires torch_scatter, which could not "
            "be imported. Install it with e.g. 'pip install torch_scatter'.") from e
    try:
        from physicsnemo.models.vfgn import VFGNLearnedSimulator
        return VFGNLearnedSimulator
    except ImportError as e:
        raise ImportError(
            "PhysicsNeMoApplication.vfgn_bridge requires physicsnemo, which could not be "
            "imported. Install it with e.g. 'pip install -U nvidia-physicsnemo'.") from e


def MakeNormalizationStats(velocity_mean, velocity_std,
                           acceleration_mean, acceleration_std,
                           context_mean=None, context_std=None):
    """The duck-typed stats dict VFGN expects.

    Upstream reads `.mean`/`.std` attributes, so a plain nested dict does NOT
    work. A zero std is rejected here rather than silently producing inf/nan:
    upstream applies an epsilon to the context statistics only.
    """
    import types

    torch = _TryImportTorch()

    def entry(mean, std, name, dtype):
        mean = torch.as_tensor(numpy.asarray(mean), dtype=dtype)
        std = torch.as_tensor(numpy.asarray(std), dtype=dtype)
        if float(std.abs().min()) <= 0.0:
            raise ValueError(
                f"\"{name}\" standard deviation contains a zero entry; VFGN divides by it "
                "without an epsilon, which would silently produce NaNs. Use 1.0 for a "
                "channel that does not vary.")
        return types.SimpleNamespace(mean=mean, std=std)

    stats = {
        # positions are float64 upstream, and the stats broadcast against them
        "velocity": entry(velocity_mean, velocity_std, "velocity", torch.float64),
        "acceleration": entry(acceleration_mean, acceleration_std,
                              "acceleration", torch.float64),
    }
    if context_mean is not None:
        stats["context"] = entry(context_mean, context_std, "context", torch.float32)
    else:
        import types as _types
        stats["context"] = _types.SimpleNamespace(
            mean=torch.zeros(1), std=torch.ones(1))
    return stats


def StatsToCard(stats) -> dict:
    """Normalization statistics as JSON for a model card.

    The shipped particle path computes these (CreateParticleTrajectoryDataset)
    but never persisted them, so a model trained on normalized data was
    deployed against raw data. Round-tripping them through the card closes
    that.
    """
    return {key: {"mean": numpy.asarray(stats[key].mean).tolist(),
                  "std": numpy.asarray(stats[key].std).tolist()}
            for key in _STATS_KEYS if key in stats}


def StatsFromCard(card: dict):
    """Rebuilds the stats dict written by StatsToCard."""
    if not card or not all(key in card for key in ("velocity", "acceleration")):
        raise ValueError(
            "The model card carries no \"velocity\"/\"acceleration\" normalization "
            "statistics; VFGN cannot run without them (see StatsToCard).")
    context = card.get("context")
    return MakeNormalizationStats(
        card["velocity"]["mean"], card["velocity"]["std"],
        card["acceleration"]["mean"], card["acceleration"]["std"],
        None if context is None else context["mean"],
        None if context is None else context["std"])


def StatsFromTrajectoryDataset(dataset):
    """Stats from a CreateParticleTrajectoryDataset, which already computes them."""
    for attribute in ("feature_mean", "feature_std", "target_mean", "target_std"):
        if not hasattr(dataset, attribute):
            raise ValueError(
                f"The dataset carries no \"{attribute}\"; pass a dataset built by "
                "torch_dataset.CreateParticleTrajectoryDataset.")
    # the feature block is the velocity window, oldest first: its last three
    # channels are the most recent velocity, which is what VFGN normalizes
    velocity_mean = numpy.asarray(dataset.feature_mean)[-3:]
    velocity_std = numpy.asarray(dataset.feature_std)[-3:]
    return MakeNormalizationStats(
        velocity_mean, velocity_std,
        numpy.asarray(dataset.target_mean), numpy.asarray(dataset.target_std))


def CreateVfgnSimulator(settings: Kratos.Parameters, normalization_stats):
    """Builds a VFGNLearnedSimulator with the contracts upstream leaves implicit.

    Settings:
        {
            "predict_length"   : 1,     // future steps per call
            "num_seq"          : 5,     // position-history length (>= 3)
            "num_particle_types" : 3,
            "particle_type_embedding_size" : 16,
            "connectivity_param" : 0.015,
            "boundaries"       : [[0.0, 1.0], [0.0, 1.0], [0.0, 1.0]]
        }
    """
    VFGNLearnedSimulator = _TryImportVfgn()

    defaults = Kratos.Parameters("""{
        "predict_length"               : 1,
        "num_seq"                      : 5,
        "num_particle_types"           : 3,
        "particle_type_embedding_size" : 16,
        "connectivity_param"           : 0.015,
        "boundaries"                   : [[0.0, 1.0], [0.0, 1.0], [0.0, 1.0]]
    }""")
    settings = Kratos.Parameters("{}") if settings is None else settings.Clone()
    settings.ValidateAndAssignDefaults(defaults)

    predict_length = settings["predict_length"].GetInt()
    if predict_length < 1:
        raise ValueError(f"\"predict_length\" must be >= 1, got {predict_length}.")
    num_seq = settings["num_seq"].GetInt()
    if num_seq < 3:
        raise ValueError(
            f"\"num_seq\" must be >= 3, got {num_seq}: the model differentiates the "
            "position history twice.")

    for key in ("velocity", "acceleration"):
        if normalization_stats is None or key not in normalization_stats:
            raise ValueError(
                f"VFGN needs \"{key}\" normalization statistics; without them upstream "
                "raises deep inside the encoder. Build them with MakeNormalizationStats "
                "or StatsFromTrajectoryDataset.")

    boundaries = [[float(v) for v in settings["boundaries"][i].GetVector()]
                  for i in range(settings["boundaries"].size())]

    return VFGNLearnedSimulator(
        # the decoder emits num_dimensions and reshapes to (predict_length, 3),
        # so these must agree or the two disagree silently
        num_dimensions=3 * predict_length,
        num_seq=num_seq,
        boundaries=boundaries,
        num_particle_types=settings["num_particle_types"].GetInt(),
        particle_type_embedding_size=settings["particle_type_embedding_size"].GetInt(),
        normalization_stats=normalization_stats,
        connectivity_param=settings["connectivity_param"].GetDouble())


def _AsGraphTensors(torch, position_sequence, edge_index, particle_types, global_context):
    n_nodes = int(position_sequence.shape[0])
    edge_index = torch.as_tensor(numpy.asarray(edge_index, dtype=numpy.int64))
    senders, receivers = edge_index[0], edge_index[1]
    n_particles = torch.tensor([n_nodes])
    n_edges = torch.tensor([int(edge_index.shape[1])])
    if particle_types is None:
        particle_types = torch.zeros(n_nodes, dtype=torch.int64)
    else:
        particle_types = torch.as_tensor(
            numpy.asarray(particle_types, dtype=numpy.int64))
    if global_context is not None:
        global_context = torch.as_tensor(
            numpy.asarray(global_context, dtype=numpy.float32)).reshape(1, -1)
    return senders, receivers, n_particles, n_edges, particle_types, global_context


def RunVfgnRollout(model, position_sequence, edge_index, predict_length: int,
                   particle_types=None, global_context=None):
    """Predicts the next positions from a position history.

    Uses the public `inference()`; `forward()` is unusable in 2.2 (see the
    module docstring).

    Args:
        model: A VFGNLearnedSimulator.
        position_sequence: (N, T, 3) history, oldest first.
        edge_index: (2, E) graph, e.g. from particle_bridge - the model builds
            no graph of its own despite its `graph_mode` argument.
        predict_length: Future steps; must match the model's construction.

    Returns:
        (N, predict_length, 3) float64 numpy array of predicted positions.
    """
    torch = _TryImportTorch()

    position_sequence = torch.as_tensor(
        numpy.asarray(position_sequence, dtype=numpy.float64))
    if position_sequence.dim() != 3 or position_sequence.shape[2] != 3:
        raise ValueError(
            f"position_sequence must be (N, T, 3), got {tuple(position_sequence.shape)}.")

    senders, receivers, n_particles, n_edges, particle_types, global_context = \
        _AsGraphTensors(torch, position_sequence, edge_index, particle_types, global_context)

    with torch.no_grad():
        predicted = model.inference(
            position_sequence, n_particles, n_edges, senders, receivers,
            int(predict_length), global_context, particle_types)
    return predicted.detach().cpu().numpy()


def ComputeVfgnLoss(model, position_sequence, next_positions, edge_index,
                    predict_length: int, particle_types=None, global_context=None,
                    noise_std: float = 6.7e-4):
    """One training step's loss, bypassing the broken `forward()`.

    Runs exactly the encode -> process -> decode composition `forward()` would
    have performed, with the same random-walk input noise and the same
    noise-corrected target, and returns the MSE between predicted and target
    normalized accelerations.

    Note:
        The model creates parameters lazily on its first forward, so build the
        optimizer AFTER one call of this function (or after a warm-up
        rollout) - otherwise it optimizes roughly half the model.
    """
    torch = _TryImportTorch()

    position_sequence = torch.as_tensor(
        numpy.asarray(position_sequence, dtype=numpy.float64))
    next_positions = torch.as_tensor(
        numpy.asarray(next_positions, dtype=numpy.float64))
    if next_positions.shape[1] != predict_length:
        raise ValueError(
            f"next_positions must be (N, predict_length, 3) with predict_length="
            f"{predict_length}, got {tuple(next_positions.shape)}.")

    senders, receivers, n_particles, n_edges, particle_types, global_context = \
        _AsGraphTensors(torch, position_sequence, edge_index, particle_types, global_context)

    noise = model.get_random_walk_noise_for_position_sequence(
        position_sequence, float(noise_std))
    noisy_sequence = position_sequence + noise

    graph = model.EncodingFeature(
        noisy_sequence, n_particles, n_edges, senders, receivers,
        global_context, particle_types)
    predicted = model._graph_network(*graph)

    # the target is corrected by the most recent noise so that input noise
    # cancels in the finite-difference acceleration, as upstream does
    most_recent_noise = noise[:, -1].unsqueeze(1).tile([1, predict_length, 1])
    target = model._inverse_decoder_postprocessor(
        next_positions + most_recent_noise, noisy_sequence)
    return torch.nn.functional.mse_loss(predicted, target)
