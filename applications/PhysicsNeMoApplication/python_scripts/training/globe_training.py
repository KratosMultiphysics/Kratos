"""Training loop for GLOBE, whose inputs are dicts rather than a batch.

TrainModel's loop batches a Dataset of (inputs, targets) TENSORS through a
DataLoader. GLOBE takes a dict of boundary meshes, a dict of reference
lengths and a point cloud, and its "batch" is one case; a collate function
cannot stack meshes. So it gets its own small loop, which is otherwise the
same shape as TrainModel's: optimizer, epochs, seeded, mean loss per epoch
returned, model left in eval mode.

A case is a tuple

    (boundary_meshes, reference_lengths, coordinates, targets)

as globe_bridge.BuildGlobeBoundaryMeshes and the usual gathers produce it;
build the list once and reuse it, since tessellating the boundaries every
epoch costs more than the forward pass.

torch is imported lazily.
"""

import numpy

import KratosMultiphysics as Kratos
from KratosMultiphysics.PhysicsNeMoApplication.bridges import globe_bridge
from KratosMultiphysics.PhysicsNeMoApplication.deployment import model_registry


def _TryImportTorch():
    try:
        import torch
        return torch
    except ImportError as e:
        raise ImportError(
            "PhysicsNeMoApplication.globe_training requires torch, which could not be "
            "imported. Install it with e.g. 'pip install torch'.") from e


def TrainGlobe(model, cases, settings: Kratos.Parameters, output_names=None):
    """Trains a GLOBE model on a list of boundary-value cases.

    Args:
        model: A GLOBE instance.
        cases: [(boundary_meshes, reference_lengths, coordinates, targets)],
            targets being (N, C_out) array-like.
        settings: Kratos Parameters; defaults:
            epochs (100), learning_rate (1e-3), optimizer ("adam"),
            loss ("mse" | "l1"), device ("cpu" - GLOBE's cluster tree is
            built per forward and the CPU path is the verified one),
            shuffle (true), echo_interval (0), seed (-1).
        output_names: The model's output field names in target-column
            order; defaults to the model's own output_field_ranks order.

    Returns:
        list[float]: mean training loss per epoch.
    """
    torch = _TryImportTorch()

    default_settings = Kratos.Parameters("""{
        "epochs"        : 100,
        "learning_rate" : 1e-3,
        "optimizer"     : "adam",
        "loss"          : "mse",
        "device"        : "cpu",
        "shuffle"       : true,
        "echo_interval" : 0,
        "seed"          : -1
    }""")
    settings.ValidateAndAssignDefaults(default_settings)
    if not cases:
        raise ValueError("TrainGlobe needs at least one case.")

    seed = settings["seed"].GetInt()
    if seed >= 0:
        torch.manual_seed(seed)

    device = model_registry.ResolveDevice(settings["device"].GetString())
    model = model.to(device)
    if output_names is None:
        output_names = list(getattr(model, "output_field_ranks", {}).keys())
        if not output_names:
            raise ValueError(
                "output_names is None and the model exposes no output_field_ranks; "
                "name the output fields explicitly.")

    optimizer_name = settings["optimizer"].GetString()
    learning_rate = settings["learning_rate"].GetDouble()
    if optimizer_name == "adam":
        optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate)
    elif optimizer_name == "sgd":
        optimizer = torch.optim.SGD(model.parameters(), lr=learning_rate)
    else:
        raise ValueError(
            f"Unsupported optimizer \"{optimizer_name}\". Use \"adam\" or \"sgd\".")
    loss_name = settings["loss"].GetString()
    if loss_name == "mse":
        loss_fn = torch.nn.functional.mse_loss
    elif loss_name == "l1":
        loss_fn = torch.nn.functional.l1_loss
    else:
        raise ValueError(f"Unsupported loss \"{loss_name}\". Use \"mse\" or \"l1\".")

    echo_interval = settings["echo_interval"].GetInt()
    shuffle = settings["shuffle"].GetBool()
    generator = numpy.random.default_rng(None if seed < 0 else seed)

    history = []
    model.train()
    for epoch in range(settings["epochs"].GetInt()):
        order = generator.permutation(len(cases)) if shuffle else range(len(cases))
        epoch_loss = 0.0
        for index in order:
            boundary_meshes, reference_lengths, coordinates, targets = cases[index]
            optimizer.zero_grad()
            prediction = globe_bridge.RunGlobeForward(
                model, coordinates, boundary_meshes, reference_lengths,
                output_names, enable_grad=True)
            target = torch.as_tensor(
                numpy.asarray(targets), dtype=prediction.dtype, device=prediction.device)
            if target.ndim == 1:
                target = target[:, None]
            loss = loss_fn(prediction, target)
            loss.backward()
            optimizer.step()
            epoch_loss += loss.item()
        history.append(epoch_loss / len(cases))
        if echo_interval > 0 and (epoch + 1) % echo_interval == 0:
            Kratos.Logger.PrintInfo(
                "TrainGlobe",
                f"epoch {epoch + 1}/{settings['epochs'].GetInt()}: "
                f"loss = {history[-1]:.6e}")
    model.eval()
    return history
