"""Mesh-reduced temporal attention over ROM trajectories.

Pairs RomApplication's POD reduction (rom_bridge) with physicsnemo's
decoder-only temporal transformer
(``physicsnemo.models.mesh_reduced.temporal_model.Sequence_Model``): the POD
basis is the (linear, exact-ordering) encoder, the attention model learns
the DYNAMICS of the reduced coordinates q(t), and predicted trajectories
reconstruct to full-order fields via
``rom_bridge.ReconstructFromReducedSpace``.

Alignment contract (important): a context token is ALWAYS prepended - a
zeros (B, 1, 1) tensor when the user has no case parameters. With it, output
slot i of ``forward(z, context)`` predicts ``z_{i+1}``; without one the
model's internal ``[:, 1:]`` slice silently drops the first-step prediction.
Training is teacher-forced (``mse(model(z[:, :-1], ctx), z[:, 1:])``) and
autoregressive rollout goes through the model's own ``sample``.

Checkpointing note: Sequence_Model is a plain torch.nn.Module (not a
physicsnemo Module) and TorchScript cannot script it, so
``training_utils.SaveTrainedModel`` / ``model_registry.LoadModel`` do not
apply - use Save/LoadRomTemporalModel (state_dict + settings) below.

torch and physicsnemo are optional runtime dependencies, imported lazily.
"""

import types

import numpy

import KratosMultiphysics as Kratos
from KratosMultiphysics.PhysicsNeMoApplication.deployment import model_registry
def _TryImportTorch():
    try:
        import torch
        return torch
    except ImportError as e:
        raise ImportError(
            "PhysicsNeMoApplication.rom_temporal requires torch, which could not be "
            "imported. Install it with e.g. 'pip install torch'.") from e


def _TryImportSequenceModel():
    try:
        from physicsnemo.models.mesh_reduced.temporal_model import Sequence_Model
        return Sequence_Model
    except ImportError as e:
        raise ImportError(
            "PhysicsNeMoApplication.rom_temporal requires physicsnemo, which could not be "
            "imported. Install it with e.g. 'pip install nvidia-physicsnemo'.") from e


_DEFAULT_SETTINGS = """{
    "input_dim"                  : 0,
    "context_dim"                : 1,
    "num_layers_decoder"         : 3,
    "num_heads"                  : 8,
    "dim_feedforward_scale"      : 4,
    "num_layers_context_encoder" : 2,
    "num_layers_input_encoder"   : 2,
    "num_layers_output_encoder"  : 2,
    "dropout_rate"               : 0.0,
    "device"                     : "auto"
}"""


def CreateSequenceModel(settings: Kratos.Parameters, device=None):
    """Creates a Sequence_Model over reduced coordinates.

    Args:
        settings: Kratos Parameters; defaults:
            input_dim (REQUIRED > 0; the basis's n_modes), context_dim (1),
            num_layers_decoder (3), num_heads (8), dim_feedforward_scale (4),
            num_layers_context_encoder/input_encoder/output_encoder (2),
            dropout_rate (0.0), device ("auto").
        device: Overrides the settings device when given.

    Returns:
        The model on the resolved device. Its ``dist`` handle is the
        initialized physicsnemo DistributedManager when one exists, else a
        minimal device-only shim (the model only reads ``dist.device``).
    """
    Sequence_Model = _TryImportSequenceModel()

    settings.ValidateAndAssignDefaults(Kratos.Parameters(_DEFAULT_SETTINGS))
    input_dim = settings["input_dim"].GetInt()
    if input_dim < 1:
        raise ValueError(f"\"input_dim\" must be >= 1 (the number of ROM modes), got {input_dim}.")
    num_heads = settings["num_heads"].GetInt()
    if input_dim % num_heads != 0:
        raise ValueError(
            f"\"input_dim\" ({input_dim}, the number of ROM modes) must be divisible by "
            f"\"num_heads\" ({num_heads}) - the attention embedding dimension is input_dim.")

    if device is None:
        device = model_registry.ResolveDevice(settings["device"].GetString())

    from physicsnemo.distributed.manager import DistributedManager
    if DistributedManager.is_initialized():
        dist = DistributedManager()
    else:
        dist = types.SimpleNamespace(device=device)

    model = Sequence_Model(
        input_dim=input_dim,
        input_context_dim=settings["context_dim"].GetInt(),
        dist=dist,
        dropout_rate=settings["dropout_rate"].GetDouble(),
        num_layers_decoder=settings["num_layers_decoder"].GetInt(),
        num_heads=settings["num_heads"].GetInt(),
        dim_feedforward_scale=settings["dim_feedforward_scale"].GetInt(),
        num_layers_context_encoder=settings["num_layers_context_encoder"].GetInt(),
        num_layers_input_encoder=settings["num_layers_input_encoder"].GetInt(),
        num_layers_output_encoder=settings["num_layers_output_encoder"].GetInt())
    return model.to(device)


def CreateRomTrajectoryDataset(q_trajectories, contexts=None):
    """Creates a torch Dataset of ROM-coordinate trajectories.

    Args:
        q_trajectories: (S, T, M) array, or a sequence of (T, M) arrays with
            EQUAL T (project snapshot series with
            rom_bridge.ProjectToReducedSpace and transpose to (T, M)).
        contexts: Optional per-trajectory case parameters, (S, C) or
            (S, 1, C); None uses zeros (S, 1, 1) - a context token is always
            fed to the model (see the module docstring).

    Returns:
        A torch.utils.data.Dataset yielding (z (T, M) float32, ctx (1, C)
        float32).
    """
    torch = _TryImportTorch()

    trajectories = [numpy.asarray(q, dtype=numpy.float32) for q in q_trajectories]
    lengths = {q.shape for q in trajectories}
    if len(lengths) != 1 or trajectories[0].ndim != 2:
        raise ValueError(
            f"All trajectories must share one (T, M) shape; got {sorted(lengths)}.")
    stacked = numpy.stack(trajectories)  # (S, T, M)

    if contexts is None:
        context_array = numpy.zeros((len(trajectories), 1, 1), dtype=numpy.float32)
    else:
        context_array = numpy.asarray(contexts, dtype=numpy.float32)
        if context_array.ndim == 2:
            context_array = context_array[:, None, :]
        if context_array.shape[0] != len(trajectories) or context_array.ndim != 3:
            raise ValueError(
                f"contexts must be (S, C) or (S, 1, C) with S = {len(trajectories)}; got "
                f"{list(numpy.asarray(contexts).shape)}.")

    class RomTrajectoryDataset(torch.utils.data.Dataset):
        def __len__(self):
            return stacked.shape[0]

        def __getitem__(self, index):
            return torch.from_numpy(stacked[index]), torch.from_numpy(context_array[index])

    return RomTrajectoryDataset()


def TrainRomTemporalModel(model, dataset, settings: Kratos.Parameters):
    """Teacher-forced training of a Sequence_Model on ROM trajectories.

    Loss per batch: mse(model(z[:, :-1], ctx), z[:, 1:]) - output slot i
    predicts z_{i+1} because the context token is prepended.

    Args:
        model: A Sequence_Model (CreateSequenceModel).
        dataset: A CreateRomTrajectoryDataset dataset.
        settings: Kratos Parameters; defaults:
            epochs (100), batch_size (8), learning_rate (1e-3),
            device ("auto"), shuffle (true), echo_interval (0), seed (-1).

    Returns:
        list[float]: mean training loss per epoch; model ends in eval mode.
    """
    torch = _TryImportTorch()

    default_settings = Kratos.Parameters("""{
        "epochs"        : 100,
        "batch_size"    : 8,
        "learning_rate" : 1e-3,
        "device"        : "auto",
        "shuffle"       : true,
        "echo_interval" : 0,
        "seed"          : -1
    }""")
    settings.ValidateAndAssignDefaults(default_settings)

    seed = settings["seed"].GetInt()
    if seed >= 0:
        torch.manual_seed(seed)

    device = model_registry.ResolveDevice(settings["device"].GetString())
    model = model.to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=settings["learning_rate"].GetDouble())
    loader = torch.utils.data.DataLoader(
        dataset,
        batch_size=settings["batch_size"].GetInt(),
        shuffle=settings["shuffle"].GetBool())
    echo_interval = settings["echo_interval"].GetInt()

    history = []
    model.train()
    for epoch in range(settings["epochs"].GetInt()):
        epoch_loss = 0.0
        batches = 0
        for z, context in loader:
            z = z.to(device)
            context = context.to(device)
            optimizer.zero_grad()
            prediction = model(z[:, :-1], context)
            loss = torch.nn.functional.mse_loss(prediction, z[:, 1:])
            loss.backward()
            optimizer.step()
            epoch_loss += loss.item()
            batches += 1
        history.append(epoch_loss / max(batches, 1))
        if echo_interval > 0 and (epoch + 1) % echo_interval == 0:
            Kratos.Logger.PrintInfo(
                "TrainRomTemporalModel",
                f"epoch {epoch + 1}/{settings['epochs'].GetInt()}: loss = {history[-1]:.6e}")
    model.eval()
    return history


def PredictRomTrajectory(model, initial_q, steps: int, context=None) -> numpy.ndarray:
    """Autoregressive rollout of the reduced coordinates.

    Args:
        model: A trained Sequence_Model.
        initial_q: The prompt - (M,) one state or (T0, M) several.
        steps: Number of future states to generate.
        context: Optional case parameters, (C,) or (1, C); None uses the
            zeros token.

    Returns:
        (T0 + steps, M) float64 numpy - prompt plus generated states. Feed
        rows (transposed to (M, T)) to rom_bridge.ReconstructFromReducedSpace.
    """
    torch = _TryImportTorch()

    prompt = numpy.asarray(initial_q, dtype=numpy.float32)
    if prompt.ndim == 1:
        prompt = prompt[None, :]
    if prompt.ndim != 2:
        raise ValueError(f"initial_q must be (M,) or (T0, M); got {list(prompt.shape)}.")

    if context is None:
        context_array = numpy.zeros((1, 1, 1), dtype=numpy.float32)
    else:
        context_array = numpy.asarray(context, dtype=numpy.float32).reshape(1, 1, -1)

    parameter = next(model.parameters())
    device = parameter.device
    z0 = torch.from_numpy(prompt)[None].to(device)             # (1, T0, M)
    ctx = torch.from_numpy(context_array).to(device)           # (1, 1, C)
    trajectory = model.sample(z0, int(steps), ctx)             # (1, T0+steps, M), no_grad inside
    return numpy.asarray(trajectory[0].cpu().to(torch.float64).numpy())


def SaveRomTemporalModel(model, settings: Kratos.Parameters, checkpoint_file) -> None:
    """Saves a Sequence_Model as {settings json, state_dict}.

    Sequence_Model is neither TorchScript-scriptable nor a physicsnemo
    Module, so the generic SaveTrainedModel/LoadModel paths do not apply.
    """
    torch = _TryImportTorch()
    torch.save({
        "settings": settings.WriteJsonString(),
        "state_dict": model.state_dict(),
    }, str(checkpoint_file))


def LoadRomTemporalModel(checkpoint_file, device="cpu"):
    """Rebuilds a Sequence_Model saved by SaveRomTemporalModel.

    Returns:
        (model, settings): the model in eval mode on the device, and the
        Kratos Parameters it was built with.
    """
    torch = _TryImportTorch()
    checkpoint = torch.load(str(checkpoint_file), map_location="cpu", weights_only=False)
    settings = Kratos.Parameters(checkpoint["settings"])
    model = CreateSequenceModel(settings, device=torch.device(device))
    model.load_state_dict(checkpoint["state_dict"])
    return model.eval(), settings
