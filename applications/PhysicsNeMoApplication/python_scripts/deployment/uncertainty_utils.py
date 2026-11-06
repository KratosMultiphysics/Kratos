"""Uncertainty quantification helpers for deployed surrogates.

MonteCarloPredict turns any model containing dropout-like layers (standard
torch Dropout or physicsnemo.nn.ConcreteDropout) into an uncertainty
estimator: N stochastic forward passes with only the dropout layers in train
mode, returning the prediction mean and standard deviation. Combined with
the deployment processes' "uncertainty" settings block this generalizes the
diffusion bridge's ensemble-mean/uncertainty outputs to any deployed model.

CollectConcreteDropoutLosses/GetConcreteDropoutRates delegate to
physicsnemo.nn's concrete-dropout helpers so TrainModel can add the
regularization term that makes the dropout rates learnable.

torch is imported lazily; module import stays ML-free.
"""


def _TryImportTorch():
    try:
        import torch
        return torch
    except ImportError as e:
        raise ImportError(
            "PhysicsNeMoApplication.uncertainty_utils requires torch, which could not be "
            "imported. Install it with e.g. 'pip install torch'.") from e


def _TryImportConcreteDropout():
    try:
        from physicsnemo.nn import collect_concrete_dropout_losses, get_concrete_dropout_rates
        return collect_concrete_dropout_losses, get_concrete_dropout_rates
    except ImportError as e:
        raise ImportError(
            "Concrete-dropout support requires physicsnemo, which could not be "
            "imported. Install it with e.g. 'pip install nvidia-physicsnemo'.") from e


def _IsDropoutLike(module) -> bool:
    # covers torch.nn.Dropout/Dropout2d/..., physicsnemo.nn.ConcreteDropout,
    # and their TorchScript recompilations (original_name preserves the class)
    name = getattr(module, "original_name", type(module).__name__)
    return "Dropout" in name


def FindDropoutModules(model):
    """Returns the dropout-like submodules of a model (scripted included)."""
    return [module for module in model.modules() if _IsDropoutLike(module)]


def MonteCarloPredict(model, forward_fn, num_samples: int, seed: int = -1):
    """MC-dropout prediction: mean and std over stochastic forward passes.

    The model stays in eval mode except for its dropout-like layers, which
    are flipped to train mode for the sampling (and restored afterwards) -
    the standard MC-dropout recipe, matching ConcreteDropout's documented
    inference mode.

    Args:
        model: The model (must contain at least one dropout-like layer -
            otherwise every sample is identical and the std meaningless).
        forward_fn: Callable model -> prediction tensor; the deployment
            process's own forward (interface dispatch included).
        num_samples: Number of stochastic samples (>= 2).
        seed: >= 0 seeds torch's RNG for reproducible sampling.

    Returns:
        (mean, std): tensors with the prediction's shape.
    """
    torch = _TryImportTorch()
    if num_samples < 2:
        raise ValueError(f"num_samples must be >= 2 [ num_samples = {num_samples} ].")

    dropout_modules = FindDropoutModules(model)
    if not dropout_modules:
        raise ValueError(
            "MC-dropout uncertainty needs a model with dropout-like layers "
            "(torch.nn.Dropout*, physicsnemo.nn.ConcreteDropout); this model has none, "
            "so every sample would be identical.")

    if seed >= 0:
        torch.manual_seed(seed)
    for module in dropout_modules:
        module.train()
    try:
        samples = torch.stack([forward_fn(model) for _ in range(num_samples)])
    finally:
        for module in dropout_modules:
            module.eval()
    return samples.mean(dim=0), samples.std(dim=0)


def CollectConcreteDropoutLosses(model):
    """Sum of the regularization losses of every ConcreteDropout layer
    (physicsnemo.nn.collect_concrete_dropout_losses)."""
    collect, _ = _TryImportConcreteDropout()
    return collect(model)


def GetConcreteDropoutRates(model) -> dict:
    """{module_name: learned dropout probability} for every ConcreteDropout
    layer (physicsnemo.nn.get_concrete_dropout_rates)."""
    _, rates = _TryImportConcreteDropout()
    return rates(model)


def _TryImportGpytorch():
    try:
        import gpytorch
        return gpytorch
    except ImportError as e:
        raise ImportError(
            "PhysicsNeMoApplication.uncertainty_utils' GP head requires gpytorch, which "
            "could not be imported. Install it with 'pip install gpytorch' (or "
            "'pip install nvidia-physicsnemo[uq-extras]').") from e


def _TryImportFieldGpHead():
    """physicsnemo's field-valued variational GP head.

    Probes the always-present _GPYTORCH_AVAILABLE flag rather than importing
    the symbol: without gpytorch the package-level import raises a bare
    "cannot import name ..." with no remedy in it, and upstream's own
    actionable message only fires at construction time.
    """
    try:
        import physicsnemo.experimental.uq as experimental_uq
    except ImportError as e:
        raise ImportError(
            "PhysicsNeMoApplication.uncertainty_utils' GP head requires physicsnemo >= "
            "2.2, which could not be imported. Install it with e.g. "
            "'pip install -U nvidia-physicsnemo'.") from e

    if not getattr(experimental_uq, "_GPYTORCH_AVAILABLE", False):
        raise ImportError(
            "PhysicsNeMoApplication.uncertainty_utils' GP head requires gpytorch, which "
            "could not be imported. Install it with 'pip install gpytorch' (or "
            "'pip install nvidia-physicsnemo[uq-extras]').")
    # physicsnemo.experimental carries no API-stability guarantee (as with the
    # OOD guard); pin a version if a deployment depends on this.
    from physicsnemo.experimental.uq import FieldVariationalGPHead
    return FieldVariationalGPHead


def CreateGpHead(input_dim: int, n_train: int, settings=None):
    """Builds a FieldVariationalGPHead for per-node field predictions.

    Args:
        input_dim: Width of the backbone features the head consumes.
        n_train: Number of training POINTS (not geometries) - it normalizes
            the ELBO, so a wrong value silently mis-balances fit against KL.
        settings: Kratos Parameters:
            { "num_tasks": 1, "n_inducing": 64, "mlp_hidden": [32],
              "confidence_z": 1.96 }

    Returns:
        The (untrained) head.
    """
    import KratosMultiphysics as Kratos

    FieldVariationalGPHead = _TryImportFieldGpHead()

    defaults = Kratos.Parameters("""{
        "num_tasks"    : 1,
        "n_inducing"   : 64,
        "mlp_hidden"   : [32],
        "confidence_z" : 1.96
    }""")
    settings = Kratos.Parameters("{}") if settings is None else settings.Clone()
    settings.ValidateAndAssignDefaults(defaults)

    hidden = [int(width) for width in settings["mlp_hidden"].GetVector()]
    return FieldVariationalGPHead(
        int(input_dim), n_train=int(n_train),
        num_tasks=settings["num_tasks"].GetInt(),
        n_inducing=settings["n_inducing"].GetInt(),
        mlp_hidden=hidden or None,
        confidence_z=settings["confidence_z"].GetDouble())


def FitGpHead(features, targets, settings=None, head=None):
    """Fits a GP head on a trained backbone's features.

    physicsnemo ships the head but no training loop, and is explicit that
    skipping the recipe collapses the variance rather than merely costing
    accuracy. This implements it: inducing points seeded from REAL features,
    an auxiliary MSE on the posterior mean (the ELBO alone can buy
    likelihood by inflating variance), and a KL ramp.

    Args:
        features: (N, D) backbone features - the latent embedding, not the
            raw inputs and not the model's outputs.
        targets: (N, T) values to predict.
        settings: Kratos Parameters, the CreateGpHead keys plus
            { "epochs": 200, "learning_rate": 0.05, "mse_weight": 1.0,
              "kl_ramp_epochs": 50, "seed": -1 }
        head: Optional pre-built head (otherwise one is created).

    Returns:
        (head, history) - history being the per-epoch negative ELBO.

    Note:
        Features are cast to float32: the head's feature MLP is float32 even
        when use_double keeps the GP internals in float64, so float64
        features raise a dtype mismatch.
    """
    import KratosMultiphysics as Kratos

    torch = _TryImportTorch()

    defaults = Kratos.Parameters("""{
        "num_tasks"       : 1,
        "n_inducing"      : 64,
        "mlp_hidden"      : [32],
        "confidence_z"    : 1.96,
        "epochs"          : 200,
        "learning_rate"   : 0.05,
        "mse_weight"      : 1.0,
        "kl_ramp_epochs"  : 50,
        "seed"            : -1
    }""")
    settings = Kratos.Parameters("{}") if settings is None else settings.Clone()
    settings.ValidateAndAssignDefaults(defaults)

    seed = settings["seed"].GetInt()
    if seed >= 0:
        torch.manual_seed(seed)

    features = torch.as_tensor(features, dtype=torch.float32)
    targets = torch.as_tensor(targets, dtype=torch.float32)
    if features.dim() != 2 or targets.dim() != 2:
        raise ValueError(
            f"features must be (N, D) and targets (N, T), got {tuple(features.shape)} "
            f"and {tuple(targets.shape)}.")
    if features.shape[0] != targets.shape[0]:
        raise ValueError(
            f"features and targets disagree on the point count: {features.shape[0]} "
            f"vs {targets.shape[0]}.")

    if head is None:
        head_settings = settings.Clone()
        for key in ("epochs", "learning_rate", "mse_weight", "kl_ramp_epochs", "seed"):
            head_settings.RemoveValue(key)
        if head_settings["num_tasks"].GetInt() != targets.shape[1]:
            head_settings["num_tasks"].SetInt(int(targets.shape[1]))
        head = CreateGpHead(int(features.shape[1]), int(features.shape[0]), head_settings)

    # Seed the inducing points from real features - the random default "sits
    # nowhere near the backbone's feature distribution". These are RAW
    # features: the head applies its own DKL transform internally.
    n_inducing = settings["n_inducing"].GetInt()
    selection = torch.randperm(features.shape[0])[:min(n_inducing, features.shape[0])]
    head.set_inducing_points(features[selection].detach())

    optimizer = torch.optim.Adam(head.parameters(),
                                 lr=settings["learning_rate"].GetDouble())
    mse_weight = settings["mse_weight"].GetDouble()
    ramp = max(1, settings["kl_ramp_epochs"].GetInt())

    history = []
    head.train()
    for epoch in range(settings["epochs"].GetInt()):
        optimizer.zero_grad()
        mean, negative_elbo = head.forward_and_loss(
            features, targets, beta=min(1.0, (epoch + 1) / ramp))
        loss = negative_elbo
        if mse_weight > 0.0:
            loss = loss + mse_weight * torch.nn.functional.mse_loss(mean, targets)
        loss.backward()
        optimizer.step()
        history.append(float(negative_elbo))
    head.eval()
    return head, history


def PredictWithGpHead(head, features, epistemic_only: bool = False):
    """(mean, std) from a fitted GP head, in this application's float64.

    Args:
        head: A fitted head.
        features: (N, D) backbone features.
        epistemic_only: Return the epistemic standard deviation (the latent
            GP term, excluding the likelihood noise floor) instead of the
            total - the "where is the model uncertain" signal.

    Returns:
        (mean (N, T), std (N, T)) float64 tensors.
    """
    torch = _TryImportTorch()

    features = torch.as_tensor(features, dtype=torch.float32)
    with torch.no_grad():
        prediction = head.predict(features)
    variance = (prediction.epistemic_variance if epistemic_only else prediction.variance)
    return (prediction.mean.detach().to(torch.float64),
            variance.detach().clamp_min(0.0).sqrt().to(torch.float64))


def SaveGpHead(head, head_file, config=None) -> None:
    """Writes a fitted head to a sidecar next to its checkpoint.

    Mirrors ood_guard_utils.SaveGuard. A GP head cannot go through
    SaveTrainedModel: gpytorch modules are not TorchScript-scriptable.
    """
    torch = _TryImportTorch()
    payload = {"config": dict(config or {}), "state_dict": head.state_dict()}
    torch.save(payload, str(head_file))


def LoadGpHead(head_file, input_dim: int = 0, n_train: int = 0, settings=None):
    """Rebuilds a head saved by SaveGpHead.

    The stored config supplies input_dim/n_train when they were saved with
    it; the explicit arguments override.
    """
    torch = _TryImportTorch()

    payload = torch.load(str(head_file), weights_only=False)
    config = payload.get("config", {})
    resolved_input_dim = int(input_dim or config.get("input_dim", 0))
    resolved_n_train = int(n_train or config.get("n_train", 1))
    if not resolved_input_dim:
        raise ValueError(
            f"\"{head_file}\" carries no input_dim in its config; pass input_dim "
            "explicitly.")

    if settings is None and config:
        import KratosMultiphysics as Kratos
        settings = Kratos.Parameters("{}")
        for key in ("num_tasks", "n_inducing", "confidence_z"):
            if key in config:
                entry = settings.AddEmptyValue(key)
                if isinstance(config[key], float):
                    entry.SetDouble(float(config[key]))
                else:
                    entry.SetInt(int(config[key]))
        if "mlp_hidden" in config:
            settings.AddEmptyValue("mlp_hidden").SetVector(
                Kratos.Vector([float(w) for w in config["mlp_hidden"]]))

    head = CreateGpHead(resolved_input_dim, resolved_n_train, settings)
    head.load_state_dict(payload["state_dict"])
    head.eval()
    return head
