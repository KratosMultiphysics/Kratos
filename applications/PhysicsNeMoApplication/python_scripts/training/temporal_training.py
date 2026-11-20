"""Temporal training schemes over Kratos transient trajectories.

The four ways to learn a time-dependent field from `(T, N, W)` state
trajectories (crash/deformation surrogates, transient thermal, any
autoregressive rollout), sharing ONE window convention with
rollout_utils.EvaluateRollout and TimeSeriesInferenceProcess: a history of
K states concatenated per node, OLDEST FIRST, into (N, K*W). A model
trained here therefore deploys through TimeSeriesInferenceProcess with no
adapter.

- "single_step"     : window -> next state. The cheap workhorse; feed the
                      dataset to training_utils.TrainModel.
- "time_conditional": initial window + normalized time -> state at t.
                      No error accumulation, but no dynamics either.
- "one_shot"        : initial window -> final state. The crash-outcome
                      framing (deformation at the end of the event).
- autoregressive    : TrainAutoregressive rolls the model forward R steps
                      through its own predictions and backpropagates
                      through the whole rollout (BPTT), optionally with
                      per-step gradient checkpointing. This is what
                      actually stabilizes long rollouts, and needs its own
                      loop: TrainModel's epoch callbacks run under
                      eval()+no_grad and cannot train.

torch is imported lazily; physicsnemo is not needed at all.
"""

import numpy

import KratosMultiphysics as Kratos

from KratosMultiphysics.PhysicsNeMoApplication.deployment import model_registry
def _TryImportTorch():
    try:
        import torch
        return torch
    except ImportError as e:
        raise ImportError(
            "PhysicsNeMoApplication.temporal_training requires torch, which could not "
            "be imported. Install it with e.g. 'pip install torch'.") from e


def _AsTrajectories(states):
    """Normalizes one (T, N, W) trajectory or a sequence of them to a list."""
    if isinstance(states, numpy.ndarray) and states.ndim == 3:
        trajectories = [states]
    else:
        trajectories = [numpy.asarray(trajectory, dtype=float) for trajectory in states]
    for trajectory in trajectories:
        if trajectory.ndim != 3:
            raise ValueError(
                f"Each trajectory must be (T, N, W), got shape {trajectory.shape}.")
    return trajectories


def MakeWindow(history):
    """(N, K*W) window from K states, OLDEST FIRST - the shared convention.

    Args:
        history: sequence of K (N, W) arrays/tensors, oldest first.
    """
    if hasattr(history[0], "detach"):
        torch = _TryImportTorch()
        return torch.cat(list(history), dim=1)
    return numpy.concatenate(list(history), axis=1)


CreateTrajectoryWindowDataset_SCHEMES = ("single_step", "time_conditional", "one_shot")


def CreateTrajectoryWindowDataset(states, settings: Kratos.Parameters):
    """Windows trajectories into (inputs, targets) pairs for TrainModel.

    Args:
        states: one (T, N, W) trajectory or a sequence of them (they may
            differ in T but must share N and W).
        settings: Kratos Parameters:
            {
                "scheme"       : "single_step" | "time_conditional" | "one_shot",
                "history_size" : 2
            }

    Returns:
        A torch Dataset whose items are row-batched (N, ...) tensors:
            single_step     : ((N, K*W) window, (N, W) next state)
            time_conditional: ((N, K*W + 1) window + normalized time in
                              [0, 1] as a broadcast extra channel,
                              (N, W) state at that time)
            one_shot        : ((N, K*W) initial window, (N, W) final state)
    """
    torch = _TryImportTorch()

    defaults = Kratos.Parameters("""{
        "scheme"       : "single_step",
        "history_size" : 2
    }""")
    settings = settings.Clone()
    settings.ValidateAndAssignDefaults(defaults)
    scheme = settings["scheme"].GetString()
    history_size = settings["history_size"].GetInt()
    if scheme not in CreateTrajectoryWindowDataset_SCHEMES:
        raise ValueError(
            f"Unknown scheme \"{scheme}\". Use one of {CreateTrajectoryWindowDataset_SCHEMES}.")
    if history_size < 1:
        raise ValueError(f"history_size must be >= 1, got {history_size}.")

    trajectories = _AsTrajectories(states)
    samples = []
    for trajectory in trajectories:
        n_steps = trajectory.shape[0]
        if n_steps <= history_size:
            raise ValueError(
                f"A trajectory has {n_steps} states but history_size is {history_size}: "
                "at least history_size + 1 states are needed.")
        if scheme == "single_step":
            for index in range(history_size, n_steps):
                window = MakeWindow(trajectory[index - history_size:index])
                samples.append((window, trajectory[index]))
        elif scheme == "one_shot":
            window = MakeWindow(trajectory[:history_size])
            samples.append((window, trajectory[-1]))
        else:  # time_conditional
            window = MakeWindow(trajectory[:history_size])
            last = n_steps - 1
            for index in range(history_size, n_steps):
                time_channel = numpy.full((window.shape[0], 1), index / last)
                samples.append((numpy.concatenate([window, time_channel], axis=1),
                                trajectory[index]))

    class TrajectoryWindowDataset(torch.utils.data.Dataset):
        def __len__(self):
            return len(samples)

        def __getitem__(self, index):
            inputs, targets = samples[index]
            return (torch.tensor(inputs, dtype=torch.float32),
                    torch.tensor(targets, dtype=torch.float32))

    dataset = TrajectoryWindowDataset()
    dataset.scheme = scheme
    dataset.history_size = history_size
    return dataset


def RolloutPredictions(model, initial_history, steps, checkpoint=False):
    """Differentiable autoregressive rollout - EvaluateRollout's twin.

    Same window convention (deque of K states, oldest-first concatenation),
    but inside the autograd graph so the multi-step error can be
    backpropagated (BPTT).

    Args:
        model: torch Module mapping (N, K*W) -> (N, W).
        initial_history: sequence of K (N, W) tensors, oldest first.
        steps: number of steps to roll forward.
        checkpoint: recompute each step's activations in backward
            (torch.utils.checkpoint, use_reentrant=False) - trades compute
            for memory on long rollouts.

    Returns:
        (steps, N, W) tensor of predictions, in the graph.
    """
    torch = _TryImportTorch()
    from torch.utils.checkpoint import checkpoint as torch_checkpoint

    history = list(initial_history)
    predictions = []
    for _ in range(steps):
        window = MakeWindow(history)
        if checkpoint and model.training and torch.is_grad_enabled():
            prediction = torch_checkpoint(model, window, use_reentrant=False)
        else:
            prediction = model(window)
        predictions.append(prediction)
        history = history[1:] + [prediction]
    return torch.stack(predictions)


def TrainAutoregressive(model, states, settings: Kratos.Parameters):
    """Backpropagation-through-time training on trajectory rollouts.

    Each sample is a rollout window: the model is seeded with K true states
    and then consumes its OWN predictions for `rollout_steps` steps; the
    loss is the mean squared error over the whole rollout, so gradients
    flow through every step. Long rollouts trade memory for compute via
    per-step gradient checkpointing.

    Args:
        model: torch Module mapping (N, K*W) -> (N, W).
        states: one (T, N, W) trajectory or a sequence of them.
        settings: Kratos Parameters; defaults:
            epochs (50), rollout_steps (4), history_size (2),
            learning_rate (1e-3), optimizer ("adam"|"sgd"),
            gradient_checkpointing (false), device ("auto"|"cpu"|"cuda"),
            shuffle (true), echo_interval (0 = silent), seed (-1).

    Returns:
        list[float]: mean loss per epoch. The model ends on the resolved
        device in eval() mode - the same contract as TrainModel.
    """
    torch = _TryImportTorch()

    defaults = Kratos.Parameters("""{
        "epochs"                 : 50,
        "rollout_steps"          : 4,
        "history_size"           : 2,
        "learning_rate"          : 1e-3,
        "optimizer"              : "adam",
        "gradient_checkpointing" : false,
        "device"                 : "auto",
        "shuffle"                : true,
        "echo_interval"          : 0,
        "seed"                   : -1
    }""")
    settings = settings.Clone()
    settings.ValidateAndAssignDefaults(defaults)

    epochs = settings["epochs"].GetInt()
    rollout_steps = settings["rollout_steps"].GetInt()
    history_size = settings["history_size"].GetInt()
    learning_rate = settings["learning_rate"].GetDouble()
    optimizer_name = settings["optimizer"].GetString()
    checkpoint = settings["gradient_checkpointing"].GetBool()
    shuffle = settings["shuffle"].GetBool()
    echo_interval = settings["echo_interval"].GetInt()
    seed = settings["seed"].GetInt()
    if rollout_steps < 1 or history_size < 1:
        raise ValueError("rollout_steps and history_size must be >= 1.")

    if seed >= 0:
        torch.manual_seed(seed)
    device = model_registry.ResolveDevice(settings["device"].GetString())
    model = model.to(device)

    trajectories = _AsTrajectories(states)
    windows = []
    for trajectory in trajectories:
        tensor = torch.tensor(trajectory, dtype=torch.float32, device=device)
        span = history_size + rollout_steps
        if tensor.shape[0] < span:
            raise ValueError(
                f"A trajectory has {tensor.shape[0]} states but history_size + "
                f"rollout_steps = {span} are needed.")
        for start in range(tensor.shape[0] - span + 1):
            windows.append((tensor[start:start + history_size],
                            tensor[start + history_size:start + span]))
    if not windows:
        raise ValueError("No rollout windows could be formed from the given trajectories.")

    if optimizer_name == "adam":
        optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate)
    elif optimizer_name == "sgd":
        optimizer = torch.optim.SGD(model.parameters(), lr=learning_rate)
    else:
        raise ValueError(f"Unknown optimizer \"{optimizer_name}\". Use \"adam\" or \"sgd\".")

    model.train()
    history = []
    for epoch in range(epochs):
        order = torch.randperm(len(windows)).tolist() if shuffle else range(len(windows))
        total = 0.0
        for index in order:
            seed_states, targets = windows[index]
            optimizer.zero_grad()
            predictions = RolloutPredictions(
                model, list(seed_states), rollout_steps, checkpoint=checkpoint)
            loss = torch.nn.functional.mse_loss(predictions, targets)
            loss.backward()
            optimizer.step()
            total += loss.item()
        mean_loss = total / len(windows)
        history.append(mean_loss)
        if echo_interval and (epoch + 1) % echo_interval == 0:
            Kratos.Logger.PrintInfo(
                "TrainAutoregressive", f"epoch {epoch + 1}/{epochs}: loss = {mean_loss:.6e}")

    model.eval()
    return history
