import contextlib
"""Parameters-driven training loop and checkpoint saving.

Removes the boilerplate every surrogate needs: TrainModel runs a standard
supervised loop over any (inputs, targets) Dataset (CreateNpzDataset output,
a TensorDataset, ...) configured through Kratos Parameters, and
SaveTrainedModel writes the checkpoint in one of the two formats
model_registry.LoadModel reads (physicsnemo .mdlus or TorchScript),
optionally with a model card sidecar.

torch is imported lazily; module import stays ML-free.
"""

import inspect

import KratosMultiphysics as Kratos
from KratosMultiphysics.PhysicsNeMoApplication.deployment import model_registry
def _TryImportTorch():
    try:
        import torch
        return torch
    except ImportError as e:
        raise ImportError(
            "PhysicsNeMoApplication.training_utils requires torch, which could not be "
            "imported. Install it with e.g. 'pip install torch'.") from e


def _GuardCalibrationRows(inputs):
    """A training batch as (rows, channels) guard features.

    Must mirror what the deployment processes feed GuardCheck.Check, or the
    calibrated guard scores a different quantity than it later checks:

    * nodal batches (B, N, C) - channels LAST - are checked as (N, C), so
      calibration rows are the points: (B*N, C);
    * grid batches (B, C, *spatial) - channels FIRST, ndim >= 4 - are
      checked as (prod(spatial), C) by the grid processes, so the channel
      axis is moved last before flattening;
    * plain (B, C) batches are their own rows.

    Grids used to be pooled over everything but the LAST axis, which for a
    (B, C, D, H, W) batch is the W axis - a guard calibrated on that flags
    every deployment input, in-distribution ones included.
    """
    if inputs.ndim >= 4:
        return inputs.movedim(1, -1).reshape(-1, inputs.shape[1])
    return inputs.reshape(-1, inputs.shape[-1])


def _WantsTargets(loss_term) -> bool:
    """Whether an extra loss term takes the batch targets as a 4th argument.

    Most terms grade the prediction against physics and need nothing else
    (physics_informed, differentiable_residual); a derivative-matching term
    needs the stored adjoint gradient, which travels in the targets. Rather
    than break the three-argument contract those terms were written against,
    the arity is resolved here - once, before the epoch loop, never per batch.

    A callable whose signature cannot be read (a C callable, an exotic
    wrapper) is treated as three-argument: that is the older contract, so
    guessing it is the safe direction.
    """
    try:
        parameters = list(inspect.signature(loss_term).parameters.values())
    except (TypeError, ValueError):
        return False
    if any(p.kind is inspect.Parameter.VAR_POSITIONAL for p in parameters):
        return True
    positional = [p for p in parameters
                  if p.kind in (inspect.Parameter.POSITIONAL_ONLY,
                                inspect.Parameter.POSITIONAL_OR_KEYWORD)]
    return len(positional) >= 4


_OPTIMIZERS = ("adam", "sgd", "muon")
_SCHEDULERS = ("none", "step", "cosine")
_AMP_DTYPES = ("bfloat16", "float16")
_LOGGER_BACKENDS = ("console", "mlflow", "wandb")


def _TryImportStaticCapture():
    try:
        from physicsnemo.utils import StaticCaptureTraining
        return StaticCaptureTraining
    except ImportError as e:
        raise ImportError(
            "\"static_capture\" requires physicsnemo, which could not be imported. "
            "Install it with e.g. 'pip install -U nvidia-physicsnemo'.") from e


def _TryImportCheckpointing():
    try:
        from physicsnemo.utils import load_checkpoint, save_checkpoint
        return save_checkpoint, load_checkpoint
    except ImportError as e:
        raise ImportError(
            "Resumable checkpoints require physicsnemo, which could not be imported. "
            "Install it with e.g. 'pip install -U nvidia-physicsnemo'.") from e


def _TryImportMuon():
    try:
        from physicsnemo.optim import CombinedOptimizer, Muon
        return Muon, CombinedOptimizer
    except ImportError as e:
        raise ImportError(
            "The \"muon\" optimizer requires physicsnemo >= 2.2, which could not be "
            "imported. Install it with e.g. 'pip install -U nvidia-physicsnemo'.") from e


def _TryImportLaunchLogger():
    try:
        from physicsnemo.utils import LaunchLogger
        return LaunchLogger
    except ImportError as e:
        raise ImportError(
            "\"launch_logger\" requires physicsnemo, which could not be imported. "
            "Install it with e.g. 'pip install -U nvidia-physicsnemo'.") from e


def _ReadPerformanceSettings(settings: Kratos.Parameters) -> dict:
    """The "performance" block of TrainModel, validated, as a plain dict."""
    scheduler_settings = Kratos.Parameters("{}")
    if settings.Has("scheduler"):
        scheduler_settings = settings["scheduler"].Clone()
        settings.RemoveValue("scheduler")
    settings.ValidateAndAssignDefaults(Kratos.Parameters("""{
        "amp"                  : false,
        "amp_dtype"            : "bfloat16",
        "static_capture"       : false,
        "cuda_graphs"          : false,
        "gradient_clip_norm"   : 0.0,
        "profile"              : false,
        "profile_output"       : "train_profile.json",
        "launch_logger"        : false,
        "logger_backend"       : "console",
        "checkpoint_directory" : "",
        "checkpoint_interval"  : 0,
        "resume"               : false
    }"""))
    scheduler_settings.ValidateAndAssignDefaults(Kratos.Parameters("""{
        "type"      : "none",
        "step_size" : 1,
        "gamma"     : 0.1,
        "t_max"     : 0
    }"""))

    performance = {
        "amp": settings["amp"].GetBool(),
        "amp_dtype": settings["amp_dtype"].GetString(),
        "static_capture": settings["static_capture"].GetBool(),
        "cuda_graphs": settings["cuda_graphs"].GetBool(),
        "gradient_clip_norm": settings["gradient_clip_norm"].GetDouble(),
        "profile": settings["profile"].GetBool(),
        "profile_output": settings["profile_output"].GetString(),
        "launch_logger": settings["launch_logger"].GetBool(),
        "logger_backend": settings["logger_backend"].GetString(),
        "checkpoint_directory": settings["checkpoint_directory"].GetString(),
        "checkpoint_interval": settings["checkpoint_interval"].GetInt(),
        "resume": settings["resume"].GetBool(),
        "scheduler": {
            "type": scheduler_settings["type"].GetString(),
            "step_size": scheduler_settings["step_size"].GetInt(),
            "gamma": scheduler_settings["gamma"].GetDouble(),
            "t_max": scheduler_settings["t_max"].GetInt(),
        },
    }
    if performance["amp_dtype"] not in _AMP_DTYPES:
        raise ValueError(
            f"Unsupported amp_dtype \"{performance['amp_dtype']}\". Use one of {_AMP_DTYPES}.")
    if performance["logger_backend"] not in _LOGGER_BACKENDS:
        raise ValueError(
            f"Unsupported logger_backend \"{performance['logger_backend']}\". Use one of "
            f"{_LOGGER_BACKENDS}.")
    if performance["scheduler"]["type"] not in _SCHEDULERS:
        raise ValueError(
            f"Unsupported scheduler \"{performance['scheduler']['type']}\". Use one of "
            f"{_SCHEDULERS}.")
    if performance["resume"] and not performance["checkpoint_directory"]:
        raise ValueError(
            "\"resume\" needs a \"checkpoint_directory\" to resume from.")
    return performance


def _BuildOptimizer(model, optimizer_name: str, learning_rate: float):
    """adam, sgd, or muon - the last split across parameter ranks.

    physicsnemo's Muon orthogonalizes each update matrix, so it REJECTS
    1-D parameters outright (biases, norm gains): "Muon only supports 2D
    parameters". The matrices go to Muon and everything else to Adam,
    through physicsnemo's CombinedOptimizer, which is a genuine
    torch Optimizer - schedulers and checkpoints accept it unchanged.
    """
    torch = _TryImportTorch()
    parameters = [p for p in model.parameters() if p.requires_grad]
    if optimizer_name == "adam":
        return torch.optim.Adam(parameters, lr=learning_rate)
    if optimizer_name == "sgd":
        return torch.optim.SGD(parameters, lr=learning_rate)
    if optimizer_name == "muon":
        Muon, CombinedOptimizer = _TryImportMuon()
        matrices = [p for p in parameters if p.ndim >= 2]
        others = [p for p in parameters if p.ndim < 2]
        if not matrices:
            raise ValueError(
                "\"muon\" found no parameter of rank >= 2 to optimize; it updates "
                "weight MATRICES only. Use \"adam\" for this model.")
        optimizers = [Muon(matrices, lr=learning_rate)]
        if others:
            optimizers.append(torch.optim.Adam(others, lr=learning_rate))
        return CombinedOptimizer(optimizers)
    raise ValueError(
        f"Unsupported optimizer \"{optimizer_name}\". Use one of {_OPTIMIZERS}.")


def _BuildScheduler(optimizer, scheduler_settings: dict, epochs: int):
    """A learning-rate schedule stepped once per epoch, or None."""
    torch = _TryImportTorch()
    kind = scheduler_settings["type"]
    if kind == "none":
        return None
    if kind == "step":
        return torch.optim.lr_scheduler.StepLR(
            optimizer, step_size=scheduler_settings["step_size"],
            gamma=scheduler_settings["gamma"])
    return torch.optim.lr_scheduler.CosineAnnealingLR(
        optimizer, T_max=scheduler_settings["t_max"] or epochs)


class _TrainingStep:
    """One optimizer step, eager or statically captured.

    The two paths are kept apart deliberately rather than merged: physicsnemo's
    StaticCaptureTraining owns zero_grad, backward, the scaler and the step
    itself, and accepts ONLY a physicsnemo Module. The eager path does those
    things by hand, with torch.autocast for AMP and a GradScaler only where
    it does anything (float16 on CUDA).
    """

    def __init__(self, model, optimizer, objective, performance: dict, device) -> None:
        torch = _TryImportTorch()
        self._torch = torch
        self._model = model
        self._optimizer = optimizer
        self._objective = objective
        self._clip = performance["gradient_clip_norm"]
        self.scaler = None
        self._static = None

        amp = performance["amp"]
        amp_dtype = torch.bfloat16 if performance["amp_dtype"] == "bfloat16" else torch.float16
        if performance["static_capture"]:
            StaticCaptureTraining = _TryImportStaticCapture()
            try:
                self._static = StaticCaptureTraining(
                    model=model, optim=optimizer,
                    use_graphs=performance["cuda_graphs"] and device.type == "cuda",
                    use_amp=amp, amp_type=amp_dtype,
                    gradient_clip_norm=self._clip or None)(objective)
            except ValueError as error:
                raise ValueError(
                    f"\"static_capture\" failed: {error}. physicsnemo's "
                    "StaticCaptureTraining accepts only a physicsnemo Module; wrap a "
                    "plain torch model with physicsnemo.Module.from_torch, or leave "
                    "static_capture off and use \"amp\" alone.") from error
            return

        if amp:
            self._autocast = lambda: torch.autocast(device_type=device.type, dtype=amp_dtype)
            if amp_dtype == torch.float16 and device.type == "cuda":
                self.scaler = torch.amp.GradScaler("cuda")
        else:
            self._autocast = contextlib.nullcontext

    def __call__(self, inputs, targets) -> float:
        if self._static is not None:
            return float(self._static(inputs, targets))

        torch = self._torch
        self._optimizer.zero_grad()
        with self._autocast():
            loss = self._objective(inputs, targets)
        if self.scaler is not None:
            self.scaler.scale(loss).backward()
            if self._clip > 0.0:
                self.scaler.unscale_(self._optimizer)
                torch.nn.utils.clip_grad_norm_(self._model.parameters(), self._clip)
            self.scaler.step(self._optimizer)
            self.scaler.update()
        else:
            loss.backward()
            if self._clip > 0.0:
                torch.nn.utils.clip_grad_norm_(self._model.parameters(), self._clip)
            self._optimizer.step()
        return loss.item()


def _HasCheckpoint(directory: str) -> bool:
    from pathlib import Path

    return bool(directory) and any(Path(directory).glob("checkpoint.*.pt"))


def _ProfilerContext(performance: dict):
    """torch.profiler around the whole run, exported as a chrome trace.

    torch's own profiler rather than physicsnemo's Profiler registry: the
    registry is a manager in front of this same profiler plus others, and
    configuring it well is its own subject. The trace opens in
    chrome://tracing or Perfetto.
    """
    if not performance["profile"]:
        return contextlib.nullcontext()
    torch = _TryImportTorch()
    activities = [torch.profiler.ProfilerActivity.CPU]
    if torch.cuda.is_available():
        activities.append(torch.profiler.ProfilerActivity.CUDA)
    output = performance["profile_output"]

    @contextlib.contextmanager
    def profiled():
        with torch.profiler.profile(activities=activities) as profiler:
            yield profiler
        profiler.export_chrome_trace(output)
        Kratos.Logger.PrintInfo("TrainModel", f"Wrote the training profile to \"{output}\".")

    return profiled()


def TrainModel(model, dataset, settings: Kratos.Parameters, epoch_callbacks=None,
               extra_loss_terms=None):
    """Trains a model on an (inputs, targets) dataset.

    Args:
        model: A torch.nn.Module (physicsnemo Modules included).
        dataset: A torch.utils.data.Dataset yielding (inputs, targets).
        settings: Kratos Parameters; defaults:
            epochs (100), batch_size (32), learning_rate (1e-3),
            optimizer ("adam"|"sgd"), loss ("mse"|"l1"),
            device ("auto"|"cpu"|"cuda"), shuffle (true),
            echo_interval (0 = silent, N = log every N epochs),
            seed (-1 = leave the RNG alone),
            concrete_dropout_reg_weight (0.0 = off; > 0 adds that multiple
                of the summed physicsnemo.nn.ConcreteDropout regularization
                losses to the objective, making the dropout rates learnable),
            target_channels ([] = the whole target): restrict the DATA loss
                to these target columns. Needed when the targets carry more
                than the model predicts - a sample exported with an adjoint
                sensitivity field alongside the solution has the gradient in
                its trailing columns, there for a loss term rather than for
                the model to reproduce (see sobolev_training).
            optimizer also accepts "muon": physicsnemo's Muon on the weight
                matrices and Adam on everything of rank < 2, which Muon
                rejects outright.
            performance ({} = off): the training performance layer -
                amp/amp_dtype (torch.autocast, a GradScaler only for float16
                on CUDA), static_capture/cuda_graphs (physicsnemo's
                StaticCaptureTraining; physicsnemo Modules only),
                gradient_clip_norm, profile/profile_output (a chrome trace),
                launch_logger/logger_backend (physicsnemo's LaunchLogger:
                console, mlflow or wandb), checkpoint_directory/
                checkpoint_interval/resume (physicsnemo's resumable
                checkpoints - model, optimizer, scheduler and scaler), and a
                scheduler block {type: none|step|cosine, step_size, gamma,
                t_max}. Resuming returns the history of the epochs run in
                THIS call, so a resumed run's history continues where the
                interrupted one stopped.
            ood_guard ({} = off; {"guard_file": "...", "buffer_size": 0 =
                len(dataset), "knn_k", "sensitivity"} calibrates an OOD
                guard on the training inputs during the first epoch and
                saves it to guard_file - the sidecar the deployment
                processes' "ood_guard" blocks load; see ood_guard_utils).
        epoch_callbacks: Optional iterable of callables
            cb(epoch, model, history) invoked after every epoch (plain
            Python argument - callables do not serialize into Parameters).
            The canonical physics-informed monitor: run the model on a
            held-out case, write the prediction into the case's model part,
            and log solver_residuals.ResidualEvaluator.ComputeResidualNorm()
            - the real PDE residual of the current surrogate. The plain
            ResidualEvaluator is assembled outside the autodiff graph
            (logging/ranking/early stopping only); for a gradient-carrying
            EXACT-residual loss term, use
            differentiable_residual.MakeExactResidualLossTerm via
            extra_loss_terms instead.
        extra_loss_terms: Optional iterable of callables
            term(model, inputs, prediction) -> scalar tensor, ADDED to the
            data loss every batch (gradient-carrying, unlike the plain
            evaluator above). Canonical sources:
            physics_informed.MakePhysicsLossTerm (analytic strong-form
            residuals via physicsnemo.sym) and
            differentiable_residual.MakeExactResidualLossTerm (the exact
            discrete residual through the real FEM assembly).
            A term declaring a FOURTH positional argument is called
            term(model, inputs, prediction, targets) instead - the batch
            targets, for terms that grade against stored data rather than
            against physics alone. sobolev_training.MakeSensitivityLossTerm
            is the one that needs it (the adjoint gradient rides in the
            targets). The arity is resolved once, before training starts.

    Returns:
        list[float]: mean training loss per epoch. The model ends up on the
        resolved device, in eval mode.
    """
    torch = _TryImportTorch()

    default_settings = Kratos.Parameters("""{
        "epochs"                      : 100,
        "batch_size"                  : 32,
        "learning_rate"               : 1e-3,
        "optimizer"                   : "adam",
        "loss"                        : "mse",
        "device"                      : "auto",
        "shuffle"                     : true,
        "echo_interval"               : 0,
        "seed"                        : -1,
        "concrete_dropout_reg_weight" : 0.0,
        "target_channels"             : [],
        "streaming"                   : false,
        "warm_restart"                : {},
        "ood_guard"                   : {},
        "performance"                 : {}
    }""")
    settings.ValidateAndAssignDefaults(default_settings)

    # optional concrete-dropout regularization (learnable dropout rates)
    concrete_reg_weight = settings["concrete_dropout_reg_weight"].GetDouble()

    # optional OOD-guard calibration: collect the training inputs and save
    # the guard sidecar for the deployment processes' "ood_guard" blocks
    guard_settings = settings["ood_guard"]
    guard_settings.ValidateAndAssignDefaults(Kratos.Parameters("""{
        "guard_file"  : "",
        "buffer_size" : 0,
        "knn_k"       : 10,
        "sensitivity" : 1.5
    }"""))
    guard_file = guard_settings["guard_file"].GetString()
    guard = None

    seed = settings["seed"].GetInt()
    if seed >= 0:
        torch.manual_seed(seed)

    device = model_registry.ResolveDevice(settings["device"].GetString())
    model = model.to(device)

    # Re-initialization goes here deliberately: after the device move so the
    # noise is drawn on the right device, and BEFORE the optimizer exists so
    # its moment estimates are never stale with respect to perturbed weights.
    warm_restart = settings["warm_restart"]
    if warm_restart.Has("shrink") or warm_restart.Has("perturb"):
        _ApplyWarmRestart(model, warm_restart, seed)

    optimizer = _BuildOptimizer(
        model, settings["optimizer"].GetString(), settings["learning_rate"].GetDouble())
    performance = _ReadPerformanceSettings(settings["performance"])
    epochs = settings["epochs"].GetInt()
    scheduler = _BuildScheduler(optimizer, performance["scheduler"], epochs)

    loss_name = settings["loss"].GetString()
    if loss_name == "mse":
        loss_fn = torch.nn.functional.mse_loss
    elif loss_name == "l1":
        loss_fn = torch.nn.functional.l1_loss
    else:
        raise ValueError(f"Unsupported loss \"{loss_name}\". Use \"mse\" or \"l1\".")

    streaming = settings["streaming"].GetBool()
    if streaming:
        # A live stream is single-pass and has no length: shuffling needs a
        # sampler (impossible), and a second epoch would drain an exhausted
        # queue and record a spurious zero-loss epoch.
        if settings["shuffle"].GetBool():
            raise ValueError(
                "\"shuffle\" is not possible with \"streaming\": an iterable dataset "
                "has no sampler. Set \"shuffle\": false (samples arrive in solver order).")
        if settings["epochs"].GetInt() != 1:
            raise ValueError(
                f"\"streaming\" requires \"epochs\": 1, got "
                f"{settings['epochs'].GetInt()}: the stream is consumed once, so later "
                "epochs would see an empty queue and log a false zero loss.")
        loader = torch.utils.data.DataLoader(
            dataset, batch_size=None)   # the dataset already emits whole batches
    else:
        loader = torch.utils.data.DataLoader(
            dataset,
            batch_size=settings["batch_size"].GetInt(),
            shuffle=settings["shuffle"].GetBool())
    echo_interval = settings["echo_interval"].GetInt()

    resolved_loss_terms = [(term, _WantsTargets(term)) for term in (extra_loss_terms or ())]
    target_channels = [int(v) for v in settings["target_channels"].GetVector()]

    def ComputeObjective(inputs, targets):
        prediction = model(inputs)
        loss = loss_fn(prediction,
                       targets[..., target_channels] if target_channels else targets)
        for loss_term, wants_targets in resolved_loss_terms:
            loss = loss + (loss_term(model, inputs, prediction, targets) if wants_targets
                           else loss_term(model, inputs, prediction))
        if concrete_reg_weight > 0.0:
            from KratosMultiphysics.PhysicsNeMoApplication.deployment import uncertainty_utils
            loss = loss + concrete_reg_weight * uncertainty_utils.CollectConcreteDropoutLosses(model)
        return loss

    step = _TrainingStep(model, optimizer, ComputeObjective, performance, device)

    checkpoint_directory = performance["checkpoint_directory"]
    start_epoch = 0
    if performance["resume"]:
        if _HasCheckpoint(checkpoint_directory):
            _, load_checkpoint = _TryImportCheckpointing()
            start_epoch = int(load_checkpoint(
                checkpoint_directory, models=model, optimizer=optimizer,
                scheduler=scheduler, scaler=step.scaler, device=device))
            Kratos.Logger.PrintInfo(
                "TrainModel", f"Resumed from \"{checkpoint_directory}\" at epoch {start_epoch}.")
        else:
            Kratos.Logger.PrintInfo(
                "TrainModel", f"No checkpoint in \"{checkpoint_directory}\"; starting fresh.")

    launch_logger = None
    if performance["launch_logger"]:
        launch_logger = _TryImportLaunchLogger()
        backend = performance["logger_backend"]
        launch_logger.initialize(use_mlflow=backend == "mlflow", use_wandb=backend == "wandb")

    history = []
    model.train()
    with _ProfilerContext(performance):
        for epoch in range(start_epoch, epochs):
            epoch_loss = 0.0
            batches = 0
            epoch_log = (launch_logger("train", epoch=epoch + 1) if launch_logger is not None
                         else contextlib.nullcontext())
            with epoch_log as log:
                for inputs, targets in loader:
                    inputs = inputs.to(device)
                    targets = targets.to(device)
                    loss_value = step(inputs, targets)
                    epoch_loss += loss_value
                    batches += 1
                    if log is not None:
                        log.log_minibatch({"loss": loss_value})
                    if guard_file and epoch == 0:  # calibrate on the first pass over the data
                        from KratosMultiphysics.PhysicsNeMoApplication.deployment import ood_guard_utils
                        with torch.no_grad():
                            # the guard is a CPU-resident deployment artifact (the
                            # inference processes score CPU features against it), so
                            # feed it CPU tensors regardless of the training device -
                            # collect() otherwise fails on CUDA with a device mismatch
                            rows = _GuardCalibrationRows(
                                inputs.detach().to("cpu", torch.float32))
                        if guard is None:
                            feature_width = int(rows.shape[-1])
                            rows_per_sample = max(1, rows.shape[0] // max(1, inputs.shape[0]))
                            buffer_size = guard_settings["buffer_size"].GetInt()
                            if buffer_size <= 0:
                                if streaming:
                                    raise ValueError(
                                        "The OOD guard needs an explicit \"buffer_size\" when "
                                        "\"streaming\" is set: a live stream has no len().")
                                buffer_size = min(4096, len(dataset) * rows_per_sample)
                            guard = ood_guard_utils.CreateOODGuard(
                                buffer_size, feature_width,
                                guard_settings["knn_k"].GetInt(),
                                guard_settings["sensitivity"].GetDouble())
                            guard_quota = max(1, -(-buffer_size * int(inputs.shape[0]) //
                                                   max(1, len(dataset))))
                        with torch.no_grad():
                            if rows.shape[0] > guard_quota:  # spread the buffer over the epoch
                                rows = rows[torch.randperm(rows.shape[0])[:guard_quota]]
                            sample_latents = rows.mean(dim=0, keepdim=True)
                            guard.collect(rows, sample_latents)
            if scheduler is not None:
                scheduler.step()
            history.append(epoch_loss / max(batches, 1))
            if echo_interval > 0 and (epoch + 1) % echo_interval == 0:
                Kratos.Logger.PrintInfo(
                    "TrainModel", f"epoch {epoch + 1}/{settings['epochs'].GetInt()}: loss = {history[-1]:.6e}")
            if epoch_callbacks:
                model.eval()
                with torch.no_grad():
                    for callback in epoch_callbacks:
                        callback(epoch, model, history)
                model.train()
            if (checkpoint_directory and performance["checkpoint_interval"] > 0
                    and ((epoch + 1) % performance["checkpoint_interval"] == 0 or epoch + 1 == epochs)):
                from pathlib import Path
                save_checkpoint, _ = _TryImportCheckpointing()
                Path(checkpoint_directory).mkdir(parents=True, exist_ok=True)
                save_checkpoint(checkpoint_directory, models=model, optimizer=optimizer,
                                scheduler=scheduler, scaler=step.scaler, epoch=epoch + 1)
    model.eval()
    if guard is not None:
        from KratosMultiphysics.PhysicsNeMoApplication.deployment import ood_guard_utils
        guard.compute_threshold()
        ood_guard_utils.SaveGuard(guard, guard_file)
        Kratos.Logger.PrintInfo("TrainModel", f"Saved OOD guard to \"{guard_file}\".")
    return history


def _GatherShardedStateDict(model):
    """Full, unsharded state dict for an FSDP2 model - or None if not sharded.

    FSDP2 replaces parameters with DTensors. Serializing those directly
    produces a checkpoint that reports success and then cannot be loaded
    ("aten.copy_.default got mixed torch.Tensor and DTensor"), and across
    ranks each writes only its own shard. Gathering first is the fix.

    Two details that matter:
    - a forward whose backward never ran leaves the parameters unsharded as
      plain Parameters, so reshard() must run before inspecting them, or a
      rank can silently contribute full tensors while another contributes
      shards;
    - full_state_dict=True does not by itself materialize every DTensor, so
      any survivor gets an explicit full_tensor().
    """
    torch = _TryImportTorch()
    try:
        from torch.distributed.tensor import DTensor
        from torch.distributed.checkpoint.state_dict import (
            get_model_state_dict, StateDictOptions)
    except ImportError:
        return None

    if hasattr(model, "reshard"):
        model.reshard()
    if not any(isinstance(parameter, DTensor) for parameter in model.parameters()):
        return None

    state_dict = get_model_state_dict(
        model, options=StateDictOptions(full_state_dict=True, cpu_offload=True))
    return {name: (value.full_tensor() if isinstance(value, DTensor) else value)
            for name, value in state_dict.items()}


def SaveTrainedModel(model, checkpoint_file, card=None) -> str:
    """Saves a model in a format model_registry.LoadModel can read.

    physicsnemo Modules save natively (.mdlus, checkpoint_type
    "physicsnemo"); anything else is scripted to TorchScript
    (checkpoint_type "torchscript"). When a card dict is given, the model
    card sidecar is written alongside (see model_registry.SaveModelCard).

    FSDP2-sharded models are gathered to full tensors first and written by
    rank 0 only; saving their DTensors directly yields a checkpoint that
    reports success and cannot be loaded. The gather is collective, so
    every rank must call this.

    Returns:
        The checkpoint_type string to use when loading.
    """
    torch = _TryImportTorch()
    checkpoint_file = str(checkpoint_file)

    checkpoint_type = None
    try:
        import physicsnemo
        if isinstance(model, physicsnemo.Module):
            if not checkpoint_file.endswith(".mdlus"):
                raise ValueError(
                    f"physicsnemo modules must be saved with a \".mdlus\" extension "
                    f"(got \"{checkpoint_file}\").")
            full_state_dict = _GatherShardedStateDict(model)
            if full_state_dict is None:
                model.save(checkpoint_file)
            else:
                # sharded: every rank must reach the gather above (it is
                # collective), but only rank 0 writes, and it writes the
                # gathered full tensors rather than its own DTensors
                import torch.distributed as distributed
                is_writer = (not distributed.is_initialized()
                             or distributed.get_rank() == 0)
                if is_writer:
                    model.save(checkpoint_file, _state_dict=full_state_dict)
                if distributed.is_initialized():
                    distributed.barrier()
            checkpoint_type = "physicsnemo"
    except ImportError:
        pass

    if checkpoint_type is None:
        if _GatherShardedStateDict(model) is not None:
            raise RuntimeError(
                f"Cannot save \"{checkpoint_file}\": the model has sharded (DTensor) "
                "parameters but is not a physicsnemo Module, and TorchScript cannot "
                "represent them. Unshard it first, or make it a physicsnemo Module.")
        try:
            torch.jit.script(model).save(checkpoint_file)
        except Exception as e:
            raise RuntimeError(
                f"Could not save the model to \"{checkpoint_file}\": it is not a physicsnemo "
                "Module and TorchScript scripting failed. Make the model scriptable or save "
                f"it manually. Original error: {e}") from e
        checkpoint_type = "torchscript"

    if card is not None:
        model_registry.SaveModelCard(checkpoint_file, card)
    return checkpoint_type


def ExportOnnxModel(model, sample_inputs, onnx_file, card=None) -> str:
    """Exports a trained model to an .onnx file for OnnxInferenceProcess.

    The export runs physicsnemo.deploy.onnx.export_to_onnx_stream on the
    model with the given sample inputs (a torch tensor or tuple of tensors
    with the deployment-time shapes, e.g. one gathered input batch), and
    writes the resulting byte stream. When a card dict is given, the model
    card sidecar ("<onnx_file>.card.json") is written alongside - the same
    sidecar format every deployment process validates.

    Note: some operators (e.g. the FFTs inside FNO-style models) are not
    supported by the CPU ONNX Runtime; MLP/conv models export and run
    everywhere.

    Returns:
        The onnx_file path as a string.
    """
    onnx_file = str(onnx_file)
    if not onnx_file.endswith(".onnx"):
        raise ValueError(
            f"ONNX models must be exported with a \".onnx\" extension (got \"{onnx_file}\").")

    from KratosMultiphysics.PhysicsNeMoApplication.deployment import onnx_utils
    export_to_onnx_stream = onnx_utils._TryImportOnnxExport()

    try:
        stream = export_to_onnx_stream(model, sample_inputs)
    except ModuleNotFoundError as e:
        if "onnxscript" in str(e):
            raise ImportError(
                "torch's ONNX exporter additionally requires onnxscript, which could not "
                "be imported. Install it with e.g. 'pip install onnxscript'.") from e
        raise
    with open(onnx_file, "wb") as f:
        f.write(stream)

    if card is not None:
        model_registry.SaveModelCard(onnx_file, card)
    return onnx_file


def _ApplyWarmRestart(model, settings: Kratos.Parameters, seed: int):
    """Shrink-and-perturb re-initialization before a warm restart.

    Re-seeding a trained model part-way toward its initialization keeps what
    it learned while restoring plasticity - the fix for a surrogate that has
    to absorb Kratos data from a new geometry family without forgetting the
    old one, and without the pathologies of training a converged network on
    a shifted distribution.

    Settings:
        {
            "shrink"  : 0.5,   // theta <- shrink * theta + perturb * noise
            "perturb" : 0.1,
            "noise"   : "scaled_normal",   // or "normal"
            "include_all_parameters" : false
        }

    Note:
        `perturb` is RELATIVE under the default "scaled_normal" (the noise is
        scaled by each tensor's own standard deviation) and ABSOLUTE under
        "normal", where 0.1 is enormous for typical weights.

        By default only float parameters with more than one dimension are
        touched. Upstream applies itself to everything, which crashes on
        integer parameters and halves LayerNorm/BatchNorm gains toward zero;
        "include_all_parameters" restores that behaviour if you want it.
    """
    torch = _TryImportTorch()
    try:
        from physicsnemo.nn import shrink_and_perturb_
    except ImportError as e:
        raise ImportError(
            "TrainModel's \"warm_restart\" requires physicsnemo >= 2.2 "
            "(shrink_and_perturb_ landed in 2.2), which could not be imported. Install "
            "it with e.g. 'pip install -U nvidia-physicsnemo'.") from e

    defaults = Kratos.Parameters("""{
        "shrink"                 : 0.5,
        "perturb"                : 0.1,
        "noise"                  : "scaled_normal",
        "include_all_parameters" : false
    }""")
    settings = settings.Clone()
    settings.ValidateAndAssignDefaults(defaults)

    shrink = settings["shrink"].GetDouble()
    perturb = settings["perturb"].GetDouble()
    noise = settings["noise"].GetString()
    if noise not in ("scaled_normal", "normal"):
        raise ValueError(
            f"Unsupported warm-restart noise \"{noise}\". Use \"scaled_normal\" or \"normal\".")
    if shrink < 0.0 or perturb < 0.0:
        raise ValueError(
            f"\"shrink\" and \"perturb\" must be >= 0 [ shrink = {shrink}, perturb = {perturb} ].")
    if shrink > 1.0:
        # upstream accepts this and silently amplifies every weight
        raise ValueError(
            f"\"shrink\" must be <= 1 [ shrink = {shrink} ]: values above 1 amplify the "
            "weights instead of shrinking them toward initialization.")

    generator = None
    if seed >= 0:
        generator = torch.Generator(device=next(model.parameters()).device)
        generator.manual_seed(seed)

    include = None
    if not settings["include_all_parameters"].GetBool():
        include = lambda name, parameter: (parameter.is_floating_point()
                                           and parameter.dim() > 1)
    return shrink_and_perturb_(model, shrink, perturb, noise=noise,
                               include=include, generator=generator)
