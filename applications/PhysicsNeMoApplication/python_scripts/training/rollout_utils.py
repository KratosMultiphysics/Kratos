"""Multi-step rollout evaluation for autoregressive time-series surrogates.

A next-state surrogate that looks accurate one step ahead can still drift or
blow up when fed its own predictions — the standard failure mode of
autoregressive models. EvaluateRollout runs exactly the rollout
TimeSeriesInferenceProcess would perform (same window contract: the last K
states concatenated along channels, oldest first) against a ground-truth
trajectory and returns the per-step error-growth curve.

torch is imported lazily; module import stays ML-free.
"""

import collections

import numpy


def _TryImportTorch():
    try:
        import torch
        return torch
    except ImportError as e:
        raise ImportError(
            "PhysicsNeMoApplication.rollout_utils requires torch, which could not be "
            "imported. Install it with e.g. 'pip install torch'.") from e


def EvaluateRollout(model, states, history_size: int, device: str = "cpu",
                    metric_names=None, mc_samples: int = 0, seed: int = -1):
    """Rolls a next-state model forward against a ground-truth trajectory.

    The history is seeded with the first history_size TRUE states; from then
    on the model is fed its own predictions (a genuine autoregressive
    rollout, not teacher forcing).

    Args:
        model: Model mapping (N, K*W) windows to (N, W) next states — the
            TimeSeriesInferenceProcess contract.
        states: Ground-truth trajectory, array-like of shape (T, N, W).
        history_size: The window length K (>= 2, matching the process).
        device: torch device for the forward passes.
        metric_names: Optional list of validation_metrics_process
            SUPPORTED_METRICS names evaluated per rollout step against the
            true state (e.g. ["rmse", "relative_l2"]).
        mc_samples: 0 (default) rolls deterministically. >= 2 runs MC
            dropout at every step (uncertainty_utils.MonteCarloPredict —
            the model needs dropout-like layers): the MEAN prediction is
            fed back autoregressively and the per-step std is returned,
            giving multi-step uncertainty growth alongside the error curve.
        seed: >= 0 seeds torch's RNG before an MC rollout.

    Returns:
        (predictions, errors) — predictions (T-K, N, W), the rollout from
        step K onwards, and errors (T-K,), the per-step RMS error against
        the true states. When metric_names or mc_samples are used the
        return is (predictions, errors, extras) with
        extras["per_step_metrics"] = {name: (T-K,) array} and/or
        extras["std"] = (T-K, N, W) MC standard deviations.
    """
    torch = _TryImportTorch()

    states = numpy.asarray(states, dtype=float)
    if states.ndim != 3:
        raise ValueError(f"states must have shape (T, N, W); got {states.shape}.")
    if history_size < 2:
        raise ValueError(f"history_size must be >= 2 [ history_size = {history_size} ].")
    if states.shape[0] <= history_size:
        raise ValueError(
            f"Need more than history_size={history_size} states to roll out; got {states.shape[0]}.")

    if metric_names:
        from KratosMultiphysics.PhysicsNeMoApplication.processes.validation_metrics_process import ComputeMetricValues
    if mc_samples:
        from KratosMultiphysics.PhysicsNeMoApplication.deployment import uncertainty_utils
    parameter = next(model.parameters(), None)
    dtype = parameter.dtype if parameter is not None else torch.float64

    history = collections.deque(states[:history_size], maxlen=history_size)
    predictions = []
    errors = []
    stds = []
    per_step_metrics = collections.defaultdict(list)
    with torch.no_grad():
        for step in range(history_size, states.shape[0]):
            window = numpy.concatenate(list(history), axis=1)  # (N, K*W), oldest first
            window_tensor = torch.from_numpy(window).to(device, dtype)
            if mc_samples:
                mean, std = uncertainty_utils.MonteCarloPredict(
                    model, lambda m: m(window_tensor).cpu(), mc_samples, seed)
                seed = -1  # seed only the first step; keep the samples varied
                prediction = mean.double().numpy()
                stds.append(std.double().numpy())
            else:
                prediction = model(window_tensor).cpu().double().numpy()
            predictions.append(prediction)
            errors.append(float(numpy.sqrt(numpy.mean((prediction - states[step]) ** 2))))
            if metric_names:
                values = ComputeMetricValues(
                    torch.from_numpy(prediction), torch.from_numpy(states[step]), metric_names)
                for name, value in values.items():
                    per_step_metrics[name].append(value)
            history.append(prediction)  # autoregressive: feed the prediction back

    result = (numpy.stack(predictions), numpy.array(errors))
    if metric_names or mc_samples:
        extras = {}
        if metric_names:
            extras["per_step_metrics"] = {
                name: numpy.array(values) for name, values in per_step_metrics.items()}
        if mc_samples:
            extras["std"] = numpy.stack(stds)
        result = result + (extras,)
    return result
