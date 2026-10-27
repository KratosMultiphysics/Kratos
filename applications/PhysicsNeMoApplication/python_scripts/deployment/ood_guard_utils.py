"""Out-of-distribution guardrails for deployed surrogates.

Bridges physicsnemo.experimental.guardrails.embedded.OODGuard (an
embedding-statistics module: per-channel global bounds plus a geometry-latent
kNN distance check, state held in torch buffers) to the application's
deployment processes:

- Calibration: TrainModel's optional "ood_guard" block collects the training
  inputs into a guard and saves it as a sidecar file next to the checkpoint
  (CalibrateGuardFromTensors/SaveGuard for manual pipelines).
- Deployment: every process accepting an "ood_guard" settings block loads the
  sidecar lazily (GuardCheck) and checks the gathered input features each
  inference. Upstream's check() only LOGS warnings, so the check is captured
  from the upstream logger and translated into a boolean plus a Kratos
  warning ("advisory"), a RuntimeError ("strict"), or nothing ("ignore").

Bridge semantics for a (N, C) feature tensor: the guard's global bounds see
the features as one (1, N, C) embedding (per-channel min/max over all
entities), and the geometry-kNN check sees the entity-mean (1, C) vector.

torch/physicsnemo are imported lazily; module import stays ML-free.
The guardrails live in physicsnemo.experimental - no API-stability
guarantee.
"""

import logging

import KratosMultiphysics as Kratos

_OOD_LOGGER_NAME = "physicsnemo.experimental.guardrails.embedded.ood_guard"
_GUARD_POLICIES = ("advisory", "strict", "ignore")


def _TryImportOODGuard():
    try:
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")  # experimental-namespace import warning
            from physicsnemo.experimental.guardrails.embedded import OODGuard
        return OODGuard
    except ImportError as e:
        raise ImportError(
            "PhysicsNeMoApplication.ood_guard_utils requires physicsnemo, which could not "
            "be imported. Install it with e.g. 'pip install nvidia-physicsnemo'.") from e


def _TryImportTorch():
    try:
        import torch
        return torch
    except ImportError as e:
        raise ImportError(
            "PhysicsNeMoApplication.ood_guard_utils requires torch, which could not be "
            "imported. Install it with e.g. 'pip install torch'.") from e


def CreateOODGuard(buffer_size: int, feature_width: int, knn_k: int = 10,
                   sensitivity: float = 1.5):
    """Creates an OODGuard sized for (N, feature_width) feature tensors
    (both the global-bounds and the geometry-kNN check enabled)."""
    OODGuard = _TryImportOODGuard()
    return OODGuard(
        buffer_size=buffer_size,
        global_dim=feature_width,
        geometry_embed_dim=feature_width,
        knn_k=knn_k,
        sensitivity=sensitivity)


def _SplitEmbeddings(features):
    """(N, C) features -> ((1, N, C) global embedding, (1, C) pooled latent)."""
    if features.ndim != 2:
        features = features.reshape(features.shape[0], -1)
    return features[None], features.mean(dim=0, keepdim=True)


def CalibrateGuardFromTensors(guard, samples) -> None:
    """Collects an iterable of (N, C) feature tensors (one per training
    sample) into the guard and computes the kNN threshold."""
    torch = _TryImportTorch()
    with torch.no_grad():
        for sample in samples:
            global_embedding, latent = _SplitEmbeddings(sample.detach().to(torch.float32))
            guard.collect(global_embedding, latent)
    guard.compute_threshold()


def SaveGuard(guard, path) -> None:
    """Writes a guard to a sidecar file (state_dict plus reconstruction dims)."""
    torch = _TryImportTorch()
    torch.save({
        "config": {
            "buffer_size": guard.buffer_size,
            "global_dim": None if guard.global_min is None else int(guard.global_min.shape[0]),
            "geometry_embed_dim": None if guard.geo_embeddings is None else int(guard.geo_embeddings.shape[1]),
            "knn_k": guard.knn_k,
            "sensitivity": guard.sensitivity,
        },
        "state_dict": guard.state_dict(),
    }, str(path))


def LoadGuard(path):
    """Reads a guard sidecar written by SaveGuard."""
    torch = _TryImportTorch()
    OODGuard = _TryImportOODGuard()
    payload = torch.load(str(path), map_location="cpu", weights_only=True)
    guard = OODGuard(**payload["config"])
    guard.load_state_dict(payload["state_dict"])
    return guard


def CheckFeatures(guard, features) -> list:
    """Runs guard.check on (N, C) features; returns the upstream warning
    messages (empty when the features look in-distribution).

    Upstream's check() only logs - a handler on the upstream logger captures
    the warnings so callers get an actual signal.
    """
    torch = _TryImportTorch()

    class _Capture(logging.Handler):
        def __init__(self):
            super().__init__(level=logging.WARNING)
            self.messages = []

        def emit(self, record):
            self.messages.append(record.getMessage())

    capture = _Capture()
    upstream_logger = logging.getLogger(_OOD_LOGGER_NAME)
    upstream_logger.addHandler(capture)
    try:
        with torch.no_grad():
            global_embedding, latent = _SplitEmbeddings(features.detach().to(torch.float32))
            guard.check(global_embedding, latent)
    finally:
        upstream_logger.removeHandler(capture)
    return capture.messages


class GuardCheck:
    """Per-process OOD guard: built from an "ood_guard" settings block,
    lazily loads the guard sidecar and checks gathered features.

    Settings:
        guard_file: The SaveGuard sidecar ("" disables the guard entirely).
        policy: "advisory" (default; one Kratos warning per flagged
            inference), "strict" (raise RuntimeError) or "ignore" (skip the
            check).

    The last outcome is exposed as ``last_flagged`` (None until the first
    check).
    """

    def __init__(self, settings: Kratos.Parameters) -> None:
        default_settings = Kratos.Parameters("""{
            "guard_file" : "",
            "policy"     : "advisory"
        }""")
        settings.ValidateAndAssignDefaults(default_settings)
        self.guard_file = settings["guard_file"].GetString()
        self.policy = settings["policy"].GetString()
        if self.policy not in _GUARD_POLICIES:
            raise ValueError(
                f"Unsupported OOD guard policy \"{self.policy}\". "
                f"Use one of {_GUARD_POLICIES}.")
        self.last_flagged = None
        self._guard = None

    @property
    def enabled(self) -> bool:
        return bool(self.guard_file) and self.policy != "ignore"

    def Check(self, features, tag: str) -> bool:
        """Checks (N, C) features; returns True when flagged as OOD."""
        if not self.enabled:
            return False
        if self._guard is None:
            self._guard = LoadGuard(self.guard_file)
        messages = CheckFeatures(self._guard, features)
        self.last_flagged = bool(messages)
        if messages:
            summary = " | ".join(messages)
            if self.policy == "strict":
                raise RuntimeError(
                    f"{tag}: input flagged as out-of-distribution by the OOD guard "
                    f"({summary}); the policy is \"strict\", so execution stops.")
            Kratos.Logger.PrintWarning(
                tag, f"Input flagged as out-of-distribution ({summary}); "
                "the guard is advisory, so execution continues.")
        return self.last_flagged
