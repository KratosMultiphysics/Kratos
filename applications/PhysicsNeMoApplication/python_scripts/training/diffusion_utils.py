"""Training and sampling helpers for conditional diffusion field models.

Bridges physicsnemo's diffusion stack (CorrDiff-style downscaling) to the
grid data this application exports:

- TrainDiffusionModel: an EDM loss loop over (condition, target) grid pairs
  (CreateGridPairDataset output) for a preconditioned denoiser such as
  physicsnemo.diffusion.preconditioners.EDMPrecondSuperResolution wrapping a
  SongUNet - a physicsnemo Module, so SaveTrainedModel writes a regular
  .mdlus checkpoint.
- GenerateEnsemble: repeated reverse-diffusion sampling conditioned on one
  grid, returning an (S, C_out, *spatial) ensemble whose mean is the
  prediction and whose spread is a calibrated uncertainty field.
- The CorrDiff TWO-STAGE recipe: TrainDiffusionModel's "regression" loss
  trains a deterministic mean stage (CorrDiffRegressionUNet), "residual"
  trains the denoiser on target - regression(condition) with the frozen
  stage-1 model (physicsnemo's RegressionLoss/ResidualLoss);
  TrainCorrDiffPair runs both, and RunRegressionMean +
  DiffusionInferenceProcess's "regression_settings" add the mean back at
  inference (mean shifts, ensemble spread untouched).

The same machinery covers the documented variations: TopoDiff-style
generative design (condition = constraint masks) and flow reconstruction
from sparse data (condition = masked observations) differ only in what the
condition channels contain. Volumetric denoisers (the "unet3d" interface
around physicsnemo.experimental.models.diffusion_unets.DiffusionUNet3D) run
the same train/sample path on full 5-D grids - no thin-axis squeeze.

WHICH UPSTREAM API. physicsnemo 2.2 deprecated the stack this bridge was
first written against (metrics.legacy_losses, samplers.
legacy_deterministic_sampler, preconditioners.legacy all warn about it), so
the default path is now the PROTOCOL API: WrapDiffusionModel adapts a
denoiser to physicsnemo.diffusion.DiffusionModel, TrainDiffusionModel's
"edm_sr" runs MSEDSMLoss over an EDMNoiseScheduler, and GenerateEnsemble
runs samplers.sample. The "api" setting keeps the legacy path reachable,
and "auto" (the default) uses it for the two CorrDiff losses, which have no
protocol counterpart in 2.2. WrapDenoiser (the legacy net(x, img_lr, sigma)
adapter) stays for that path.

Two things the protocol API buys, both in GenerateEnsemble:

- "guidance": diffusion posterior sampling. An already trained denoiser is
  steered at SAMPLING time toward measurements ("data_consistency") or
  toward any differentiable forward operator ("model_consistency") - the
  exact discrete FEM residual of physics.diffusion_residual_operator being
  the interesting one, since it makes the solver's own physics grade the
  generator with no retraining.
- "patching": multi-diffusion. A model trained at patch resolution samples
  a grid larger than it ever saw. 2-D only upstream.

torch and physicsnemo are optional runtime dependencies, imported lazily.
"""

import functools
import inspect
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
            "PhysicsNeMoApplication.diffusion_utils requires torch, which could not be "
            "imported. Install it with e.g. 'pip install torch'.") from e


def _TryImportPhysicsNemoDiffusion():
    try:
        from physicsnemo.diffusion.metrics.legacy_losses import (
            EDMLossSR, RegressionLoss, ResidualLoss)
        from physicsnemo.diffusion.samplers import deterministic_sampler
        return EDMLossSR, RegressionLoss, ResidualLoss, deterministic_sampler
    except ImportError as e:
        raise ImportError(
            "PhysicsNeMoApplication.diffusion_utils requires physicsnemo, which could not "
            "be imported. Install it with e.g. 'pip install nvidia-physicsnemo'.") from e


def _TryImportPhysicsNemoDiffusionProtocol():
    """physicsnemo 2.2's protocol diffusion API - the successor of the
    legacy samplers/losses/preconditioners this bridge used to run on.

    Every legacy module warns that it "will be deprecated in a future
    release"; the protocol API (noise schedulers, preconditioners, the
    DiffusionModel/Predictor/Denoiser protocols, samplers.sample and
    MSEDSMLoss) is what DPS guidance and multi-diffusion patching are
    written against, so it is the default path here.
    """
    try:
        from physicsnemo.diffusion import DiffusionModel
        from physicsnemo.diffusion.metrics import MSEDSMLoss
        from physicsnemo.diffusion.noise_schedulers import EDMNoiseScheduler
        from physicsnemo.diffusion.preconditioners import EDMPreconditioner
        from physicsnemo.diffusion.samplers import sample
        from physicsnemo.diffusion.guidance import (
            DPSScorePredictor, DataConsistencyDPSGuidance, ModelConsistencyDPSGuidance)
    except ImportError as e:
        raise ImportError(
            "The protocol diffusion API requires physicsnemo >= 2.2, which could not be "
            "imported. Install it with e.g. 'pip install -U nvidia-physicsnemo', or set "
            "\"api\" : \"legacy\" to stay on the deprecated modules.") from e
    return types.SimpleNamespace(
        DiffusionModel=DiffusionModel,
        MSEDSMLoss=MSEDSMLoss,
        EDMNoiseScheduler=EDMNoiseScheduler,
        EDMPreconditioner=EDMPreconditioner,
        sample=sample,
        DPSScorePredictor=DPSScorePredictor,
        DataConsistencyDPSGuidance=DataConsistencyDPSGuidance,
        ModelConsistencyDPSGuidance=ModelConsistencyDPSGuidance)


def _TryImportPhysicsNemoMultiDiffusion():
    """physicsnemo 2.2's multi-diffusion: a model trained on patches tiles a
    latent larger than its own training resolution. 2-D only upstream."""
    try:
        from physicsnemo.diffusion.multi_diffusion import (
            MultiDiffusionModel2D, MultiDiffusionMSEDSMLoss, MultiDiffusionPredictor)
    except ImportError as e:
        raise ImportError(
            "Multi-diffusion patching requires physicsnemo >= 2.2, which could not be "
            "imported. Install it with e.g. 'pip install -U nvidia-physicsnemo'.") from e
    return types.SimpleNamespace(
        MultiDiffusionModel2D=MultiDiffusionModel2D,
        MultiDiffusionMSEDSMLoss=MultiDiffusionMSEDSMLoss,
        MultiDiffusionPredictor=MultiDiffusionPredictor)


def _TryImportTensorDict():
    try:
        from tensordict import TensorDict
        return TensorDict
    except ImportError as e:
        raise ImportError(
            "The \"unet3d\" denoiser interface requires tensordict (a physicsnemo "
            "dependency), which could not be imported. Install it with e.g. "
            "'pip install nvidia-physicsnemo'.") from e


def WrapDenoiser(model, interface: str = "dit", out_channels: int = 0):
    """Adapts a denoiser to the EDM sampler/loss contract net(x, img_lr, sigma).

    The shipped samplers and EDMLossSR drive super-resolution-style
    denoisers as ``net(x, img_lr, sigma)``. physicsnemo.models.dit.DiT
    speaks ``dit(x, t, condition=None)`` instead; the "dit" wrapper maps
    one onto the other by concatenating the conditioning grid into the
    input channels (construct the DiT with in_channels = C_out + C_cond and
    out_channels = C_out) and broadcasting sigma to the per-sample timestep
    tensor. The wrapped module exposes img_out_channels, so GenerateEnsemble
    and DiffusionInferenceProcess work unchanged; note the raw DiT acts as
    the denoiser D(x, sigma) directly (no EDM pre/post-scaling), so train it
    through the same wrapper (TrainDiffusionModel accepts it as-is).

    The "unet3d" interface does the same for the volumetric
    physicsnemo.experimental.models.diffusion_unets.DiffusionUNet3D, which
    speaks ``unet(x, t, condition=TensorDict)``: the conditioning grid is
    NOT concatenated but passed natively as ``condition["volume"]``
    (construct the model with x_channels = C_out and vol_cond_channels =
    C_cond), sigma broadcasts to the timestep tensor exactly as for DiT.
    Latents are 5-D (B, C, D, H, W) and each spatial extent must be a power
    of 2 or a multiple of 2**(num_levels - 1) - the model validates this
    itself.

    RoPE, invalid-region masking and alternative attention backends are DiT
    construction choices (block_kwargs/attn_kwargs/attention_backend) - the
    wrapper only standardizes the forward call.

    Args:
        model: The denoiser to wrap (a DiT for interface "dit", a
            DiffusionUNet3D for "unet3d").
        interface: "dit" or "unet3d".
        out_channels: The number of predicted channels; 0 reads the DiT's
            out_channels / the DiffusionUNet3D's x_channels attribute.

    Returns:
        A torch.nn.Module with the net(x, img_lr, sigma) interface.
    """
    torch = _TryImportTorch()
    if interface not in ("dit", "unet3d"):
        raise ValueError(
            f"Unsupported denoiser interface \"{interface}\". Use \"dit\" "
            "(physicsnemo.models.dit.DiT) or \"unet3d\" (physicsnemo.experimental."
            "models.diffusion_unets.DiffusionUNet3D) - EDM-preconditioned denoisers "
            "already speak net(x, img_lr, sigma).")

    if interface == "unet3d":
        # DiffusionUNet3D predicts as many channels as its latent input
        resolved_out_channels = out_channels or getattr(model, "x_channels", 0)
        if not resolved_out_channels:
            raise ValueError(
                "\"out_channels\" is 0 and the model exposes no x_channels; "
                "set it explicitly.")
        TensorDict = _TryImportTensorDict()

        class Unet3dDenoiser(torch.nn.Module):
            # EDM interface attributes the samplers read from the net
            sigma_min = 0.0
            sigma_max = float("inf")

            def __init__(self, unet):
                super().__init__()
                self.unet = unet
                self.img_out_channels = resolved_out_channels

            @staticmethod
            def round_sigma(sigma):
                return torch.as_tensor(sigma)

            def forward(self, x, img_lr, sigma, class_labels=None, **kwargs):
                if x.ndim != 5:
                    raise ValueError(
                        "The \"unet3d\" denoiser is volumetric: expected a 5-D "
                        f"(B, C, D, H, W) latent, got shape {tuple(x.shape)}. "
                        "Planar (squeezed) grids need a 2D interface such as \"dit\".")
                parameter = next(self.unet.parameters(), None)
                dtype = parameter.dtype if parameter is not None else x.dtype
                t = torch.atleast_1d(torch.as_tensor(sigma, device=x.device, dtype=dtype))
                if t.numel() == 1:
                    t = t.expand(x.shape[0])
                condition = TensorDict(
                    {"volume": img_lr.to(dtype)}, batch_size=[x.shape[0]])
                denoised = self.unet(
                    x.to(dtype), t.reshape(x.shape[0]), condition=condition)
                return denoised.to(x.dtype)

        return Unet3dDenoiser(model)

    resolved_out_channels = out_channels or getattr(model, "out_channels", 0)
    if not resolved_out_channels:
        raise ValueError(
            "\"out_channels\" is 0 and the model exposes no out_channels; set it explicitly.")

    class DitDenoiser(torch.nn.Module):
        # EDM interface attributes the samplers read from the net
        sigma_min = 0.0
        sigma_max = float("inf")

        def __init__(self, dit):
            super().__init__()
            self.dit = dit
            self.img_out_channels = resolved_out_channels

        @staticmethod
        def round_sigma(sigma):
            return torch.as_tensor(sigma)

        def forward(self, x, img_lr, sigma, class_labels=None, **kwargs):
            # the samplers run their loop in float64; cast to the DiT's
            # parameter dtype at the boundary and back
            parameter = next(self.dit.parameters(), None)
            dtype = parameter.dtype if parameter is not None else x.dtype
            t = torch.atleast_1d(torch.as_tensor(sigma, device=x.device, dtype=dtype))
            if t.numel() == 1:
                t = t.expand(x.shape[0])
            denoised = self.dit(
                torch.cat([x, img_lr], dim=1).to(dtype), t.reshape(x.shape[0]))
            return denoised.to(x.dtype)

    return DitDenoiser(model)


_PROTOCOL_ADAPTERS = {}


def _AdapterClasses():
    """The protocol adapter classes, built once on first use.

    They are torch.nn.Module subclasses rather than physicsnemo Modules:
    physicsnemo's DiffusionModel is a runtime-checkable Protocol, so duck
    typing is enough, and a plain module keeps the wrapped checkpoint's own
    class untouched (wrap AFTER loading, exactly as the legacy WrapDenoiser
    is used).
    """
    if _PROTOCOL_ADAPTERS:
        return _PROTOCOL_ADAPTERS
    torch = _TryImportTorch()

    class LegacyDenoiserAdapter(torch.nn.Module):
        """A net(x, img_lr, sigma) denoiser as a protocol DiffusionModel.

        The whole super-resolution family - EDMPrecondSuperResolution and
        anything WrapDenoiser produces - takes the condition POSITIONALLY
        and second. The protocol passes it as the "condition" keyword, so
        the two contracts differ only in argument order and naming.
        """

        def __init__(self, net, out_channels: int = 0):
            super().__init__()
            self.net = net
            self.img_out_channels = out_channels or getattr(net, "img_out_channels", 0)

        def forward(self, x, t, condition=None, **kwargs):
            if condition is None:
                raise ValueError(
                    "The super-resolution denoiser interface is conditional: no "
                    "\"condition\" reached it. The sampler does not forward a "
                    "condition on its own - bind it into the predictor.")
            return self.net(x, condition, t, **kwargs)

    class DitDiffusionModel(torch.nn.Module):
        """physicsnemo.models.dit.DiT as a protocol DiffusionModel.

        The condition is concatenated into the input channels (construct
        the DiT with in_channels = C_out + C_cond).
        """

        def __init__(self, dit, out_channels: int = 0):
            super().__init__()
            self.dit = dit
            self.img_out_channels = out_channels or getattr(dit, "out_channels", 0)

        def forward(self, x, t, condition=None, **kwargs):
            parameter = next(self.dit.parameters(), None)
            dtype = parameter.dtype if parameter is not None else x.dtype
            t = torch.atleast_1d(torch.as_tensor(t, device=x.device, dtype=dtype))
            if t.numel() == 1:
                t = t.expand(x.shape[0])
            denoised = self.dit(
                torch.cat([x, condition], dim=1).to(dtype), t.reshape(x.shape[0]))
            return denoised.to(x.dtype)

    class TopoDiffDiffusionModel(torch.nn.Module):
        """physicsnemo.models.topodiff.TopoDiff as a protocol DiffusionModel.

        TopoDiff is built for generative design and takes its constraint
        channels NATIVELY as a second argument - forward(x, cons, timesteps)
        - rather than concatenated into the input like DiT. Construct it with
        in_channels = C_out + C_cond all the same, because its own forward
        concatenates the two before the encoder.
        """

        def __init__(self, topodiff, out_channels: int = 0):
            super().__init__()
            self.topodiff = topodiff
            self.img_out_channels = out_channels or getattr(topodiff, "out_channels", 0)

        def forward(self, x, t, condition=None, **kwargs):
            if condition is None:
                raise ValueError(
                    "The \"topodiff\" interface is conditional: it generates a design "
                    "FROM constraint channels, and no \"condition\" reached it.")
            parameter = next(self.topodiff.parameters(), None)
            dtype = parameter.dtype if parameter is not None else x.dtype
            t = torch.atleast_1d(torch.as_tensor(t, device=x.device, dtype=dtype))
            if t.numel() == 1:
                t = t.expand(x.shape[0])
            generated = self.topodiff(
                x.to(dtype), condition.to(dtype), t.reshape(x.shape[0]))
            return generated.to(x.dtype)

    class Unet3dDiffusionModel(torch.nn.Module):
        """The volumetric DiffusionUNet3D as a protocol DiffusionModel.

        The condition is NOT concatenated: it goes in natively as the
        model's TensorDict "volume" key.
        """

        def __init__(self, unet, out_channels: int = 0):
            super().__init__()
            self.unet = unet
            self.img_out_channels = out_channels or getattr(unet, "x_channels", 0)

        def forward(self, x, t, condition=None, **kwargs):
            if x.ndim != 5:
                raise ValueError(
                    "The \"unet3d\" denoiser is volumetric: expected a 5-D "
                    f"(B, C, D, H, W) latent, got shape {tuple(x.shape)}. "
                    "Planar (squeezed) grids need a 2D interface such as \"dit\".")
            TensorDict = _TryImportTensorDict()
            parameter = next(self.unet.parameters(), None)
            dtype = parameter.dtype if parameter is not None else x.dtype
            t = torch.atleast_1d(torch.as_tensor(t, device=x.device, dtype=dtype))
            if t.numel() == 1:
                t = t.expand(x.shape[0])
            volume = TensorDict({"volume": condition.to(dtype)}, batch_size=[x.shape[0]])
            denoised = self.unet(x.to(dtype), t.reshape(x.shape[0]), condition=volume)
            return denoised.to(x.dtype)

    _PROTOCOL_ADAPTERS.update({
        "edm": LegacyDenoiserAdapter,
        "dit": DitDiffusionModel,
        "unet3d": Unet3dDiffusionModel,
        "topodiff": TopoDiffDiffusionModel,
    })
    return _PROTOCOL_ADAPTERS


_DENOISER_INTERFACES = ("edm", "dit", "unet3d", "topodiff", "protocol")


def WrapDiffusionModel(model, interface: str = "edm", out_channels: int = 0,
                       preconditioner: str = "none", sigma_data: float = 0.5):
    """Adapts a denoiser to physicsnemo 2.2's DiffusionModel protocol.

    The protocol contract is ``model(x, t, condition=None)`` with t of
    shape (B,) - what MSEDSMLoss, samplers.sample and the DPS guidance
    machinery all drive. This is the protocol-API counterpart of
    WrapDenoiser (which targets the legacy net(x, img_lr, sigma) contract
    and is kept for "api" : "legacy").

    Args:
        model: The denoiser. For "edm" anything speaking
            net(x, img_lr, sigma) - EDMPrecondSuperResolution, or a
            WrapDenoiser result. For "dit" a physicsnemo.models.dit.DiT,
            for "unet3d" a DiffusionUNet3D.
        interface: "edm", "dit", "unet3d" or "topodiff".
        out_channels: 0 reads the model's own img_out_channels /
            out_channels / x_channels attribute.
        preconditioner: "none" (default) leaves the model as the raw
            x0-predictor it was trained as; "edm" wraps it in
            EDMPreconditioner (EDM skip/output scaling). It is REFUSED for
            "edm", whose models are preconditioned already - wrapping one
            again would precondition twice and silently produce garbage.
            Train and deploy must agree on this choice.
        sigma_data: The EDM sigma_data of the preconditioner.

    Returns:
        A torch.nn.Module satisfying physicsnemo.diffusion.DiffusionModel
        and exposing img_out_channels.
    """
    adapters = _AdapterClasses()
    if interface not in adapters:
        raise ValueError(
            f"Unsupported denoiser interface \"{interface}\". Use \"edm\" "
            "(net(x, img_lr, sigma) denoisers), \"dit\" (physicsnemo.models.dit.DiT) "
            "\"unet3d\" (physicsnemo.experimental.models.diffusion_unets."
            "DiffusionUNet3D) or \"topodiff\" (physicsnemo.models.topodiff.TopoDiff).")
    if preconditioner not in ("none", "edm"):
        raise ValueError(
            f"Unsupported preconditioner \"{preconditioner}\". Use \"none\" (the model "
            "predicts x0 directly) or \"edm\" (EDMPreconditioner scaling).")
    if preconditioner == "edm" and interface == "edm":
        raise ValueError(
            "The \"edm\" interface is already an EDM-preconditioned denoiser; wrapping "
            "it in EDMPreconditioner would precondition it twice. Use "
            "\"preconditioner\" : \"none\" here.")

    adapted = adapters[interface](model, out_channels)
    if not adapted.img_out_channels:
        raise ValueError(
            "\"out_channels\" is 0 and the model exposes no img_out_channels / "
            "out_channels / x_channels; set it explicitly.")
    if preconditioner == "none":
        return adapted

    protocol = _TryImportPhysicsNemoDiffusionProtocol()
    preconditioned = protocol.EDMPreconditioner(adapted, sigma_data=sigma_data)
    # the preconditioner does not forward the wrapped model's attributes
    preconditioned.img_out_channels = adapted.img_out_channels
    return preconditioned


def _AsProtocolModel(model, interface: str = "edm", out_channels: int = 0,
                     preconditioner: str = "none", sigma_data: float = 0.5):
    """Adapts a model to the protocol unless it already speaks it.

    The interface is ALWAYS explicit, never sniffed from the signature.
    physicsnemo.models.dit.DiT and the volumetric DiffusionUNet3D both
    declare a "condition" parameter, so a signature probe would call them
    "protocol models" - but DiT's condition is a (B, condition_dim) label
    VECTOR (a conditioning grid belongs in its input channels) and
    DiffusionUNet3D's is a TensorDict. Guessing would feed each the wrong
    thing, in one case without any shape error.
    """
    if interface not in _DENOISER_INTERFACES:
        raise ValueError(
            f"Unsupported denoiser interface \"{interface}\". Supported: "
            f"{', '.join(_DENOISER_INTERFACES)} (\"protocol\" = the model already "
            "speaks model(x, t, condition=...)).")
    if interface == "protocol":
        return model
    return WrapDiffusionModel(model, interface, out_channels, preconditioner, sigma_data)


_DIFFUSION_APIS = ("auto", "protocol", "legacy")
# the two CorrDiff losses have no protocol-API counterpart in physicsnemo 2.2
_LEGACY_ONLY_LOSSES = ("regression", "residual")


def _ResolveDiffusionApi(api: str, loss_name: str) -> str:
    """Which physicsnemo diffusion API a call runs on.

    "auto" keeps every shipped recipe working unchanged: the conditional
    EDM loss moves to the protocol API, the two CorrDiff stages stay on the
    legacy losses because physicsnemo 2.2 ships no protocol equivalent
    (metrics/losses.py has only MSEDSMLoss and WeightedMSEDSMLoss).
    """
    if api not in _DIFFUSION_APIS:
        raise ValueError(
            f"Unsupported diffusion api \"{api}\". Use \"auto\" (protocol where it "
            "exists), \"protocol\" or \"legacy\".")
    if api == "auto":
        return "legacy" if loss_name in _LEGACY_ONLY_LOSSES else "protocol"
    if api == "protocol" and loss_name in _LEGACY_ONLY_LOSSES:
        raise ValueError(
            f"Loss \"{loss_name}\" is a CorrDiff stage and physicsnemo 2.2 ships no "
            "protocol-API equivalent of RegressionLoss/ResidualLoss. Use \"api\" : "
            "\"auto\" (or \"legacy\") for the two-stage recipe.")
    return api


def TrainDiffusionModel(model, dataset, settings: Kratos.Parameters, regression_model=None):
    """Trains a conditional diffusion (or CorrDiff-stage) model on
    (condition, target) pairs.

    Args:
        model: The trainable model. For "edm_sr" and "residual": a denoiser
            named by "denoiser_interface" - by default the
            super-resolution contract net(x, img_lr, sigma), e.g.
            EDMPrecondSuperResolution wrapping a SongUNet. For
            "regression": a deterministic mean predictor with the
            (x_zeros, img_lr) interface -
            physicsnemo.models.diffusion_unets.CorrDiffRegressionUNet.
            2D spatial layout (C, H, W): planar Kratos cases use the
            thin-axis idiom (CreateGridPairDataset squeeze_axis).
            Volumetric (C, D, H, W) samples train with "edm_sr" on the
            protocol API, whose loss is rank-agnostic; the 2D-only CorrDiff
            losses reject them with a clear error.
        dataset: A torch Dataset yielding (condition_grid, target_grid)
            float pairs, e.g. CreateGridPairDataset output.
        settings: Kratos Parameters; defaults:
            epochs (100), batch_size (8), learning_rate (1e-4),
            device ("auto"), shuffle (true), echo_interval (0), seed (-1),
            loss ("edm_sr" | "regression" | "residual"),
            api ("auto" | "protocol" | "legacy"), denoiser_interface
            ("edm" | "dit" | "unet3d" | "protocol"), preconditioner
            ("none" | "edm"), P_mean (-1.2 for edm_sr, 0.0 for residual -
            upstream's defaults; explicit values always win), P_std (1.2),
            sigma_data (0.5), sigma_min (0.002), sigma_max (80.0),
            patching ({}).

            "api" picks the physicsnemo API. "auto" (the default) runs
            "edm_sr" on the 2.2 PROTOCOL API (EDMNoiseScheduler +
            MSEDSMLoss) and the two CorrDiff losses on the legacy ones,
            which have no protocol equivalent in 2.2. "protocol" with a
            CorrDiff loss is refused rather than silently downgraded.

            "patching" = {"patch_shape": [H, W], "patch_num": N} trains
            through physicsnemo's multi-diffusion: each step draws N random
            patches, so the model learns at patch resolution and can later
            sample a larger grid (GenerateEnsemble's "patching"). 2-D and
            "edm_sr" only.
        regression_model: REQUIRED for loss "residual": the trained
            CorrDiff regression stage. It is moved to the device, frozen
            (.eval() + requires_grad_(False) - upstream's ResidualLoss
            does NOT no_grad it) and used to compute the residual targets
            target - regression(condition) the denoiser learns. NOTE:
            ResidualLoss passes embedding_selector/global_index kwargs to
            the denoiser, so the "residual" model must wrap
            SongUNetPosEmbd (model_type="SongUNetPosEmbd"; its
            N_grid_channels, default 4, count toward img_in_channels -
            upstream CorrDiff's own sizing).

    Returns:
        list[float]: mean training loss per epoch. The model ends up on the
        resolved device, in eval mode.
    """
    torch = _TryImportTorch()

    user_set_p_mean = settings.Has("P_mean")
    patching_settings = Kratos.Parameters("{}")
    if settings.Has("patching"):
        patching_settings = settings["patching"].Clone()
        settings.RemoveValue("patching")
    default_settings = Kratos.Parameters("""{
        "epochs"             : 100,
        "batch_size"         : 8,
        "learning_rate"      : 1e-4,
        "device"             : "auto",
        "shuffle"            : true,
        "echo_interval"      : 0,
        "seed"               : -1,
        "loss"               : "edm_sr",
        "api"                : "auto",
        "denoiser_interface" : "edm",
        "preconditioner"     : "none",
        "P_mean"             : -1.2,
        "P_std"              : 1.2,
        "sigma_data"         : 0.5,
        "sigma_min"          : 0.002,
        "sigma_max"          : 80.0
    }""")
    settings.ValidateAndAssignDefaults(default_settings)
    loss_name = settings["loss"].GetString()
    if loss_name not in ("edm_sr", "regression", "residual"):
        raise ValueError(
            f"Unsupported diffusion loss \"{loss_name}\". Use \"edm_sr\" (conditional "
            "EDM), \"regression\" (CorrDiff mean stage) or \"residual\" (CorrDiff "
            "denoiser stage, needs regression_model).")
    api = _ResolveDiffusionApi(settings["api"].GetString(), loss_name)

    seed = settings["seed"].GetInt()
    if seed >= 0:
        torch.manual_seed(seed)

    device = model_registry.ResolveDevice(settings["device"].GetString())
    model = model.to(device)
    parameter = next(model.parameters(), None)
    dtype = parameter.dtype if parameter is not None else torch.float32

    # upstream's per-loss P_mean defaults differ (edm_sr -1.2, residual 0.0)
    p_mean = settings["P_mean"].GetDouble()
    if not user_set_p_mean and loss_name == "residual":
        p_mean = 0.0
    sigma_data = settings["sigma_data"].GetDouble()

    # unbatched (C, D, H, W) samples mean volumetric 5-D batches
    volumetric = False
    sample_shape = None
    try:
        sample_shape = tuple(getattr(dataset[0][1], "shape", ()))
        volumetric = len(sample_shape) == 4
    except TypeError:
        pass  # non-indexable dataset: assume the 2D image layout
    if volumetric and loss_name != "edm_sr":
        raise ValueError(
            f"Loss \"{loss_name}\" is 2D-only (upstream CorrDiff is an image recipe "
            "whose losses hard-code the 4-D batch rank); volumetric (C, D, H, W) "
            "samples train with \"edm_sr\".")
    if volumetric and api == "legacy":
        raise ValueError(
            "The legacy EDM loss hard-codes the 4-D image rank (its noise draw is "
            "randn([B, 1, 1, 1])), so volumetric (C, D, H, W) samples cannot train on "
            "it. Use \"api\" : \"protocol\" (or \"auto\"), whose MSEDSMLoss is "
            "rank-agnostic.")

    patched = patching_settings.Has("patch_shape") and patching_settings["patch_shape"].size() > 0
    if patched:
        if api != "protocol" or loss_name != "edm_sr":
            raise ValueError(
                "\"patching\" runs through physicsnemo's multi-diffusion, which exists "
                "only on the protocol API: use the default \"api\" with loss "
                "\"edm_sr\".")
        if volumetric:
            raise ValueError(
                "Multi-diffusion is 2-D upstream, so volumetric (C, D, H, W) samples "
                "cannot be patched; train them whole-grid.")

    if api == "protocol":
        protocol = _TryImportPhysicsNemoDiffusionProtocol()
        scheduler = protocol.EDMNoiseScheduler(
            sigma_min=settings["sigma_min"].GetDouble(),
            sigma_max=settings["sigma_max"].GetDouble(),
            sigma_data=sigma_data, P_mean=p_mean, P_std=settings["P_std"].GetDouble())
        protocol_model = _AsProtocolModel(
            model, settings["denoiser_interface"].GetString(),
            preconditioner=settings["preconditioner"].GetString(), sigma_data=sigma_data)
        if patched:
            multi = _TryImportPhysicsNemoMultiDiffusion()
            if sample_shape is None or len(sample_shape) != 3:
                raise ValueError(
                    "\"patching\" needs indexable (C, H, W) samples to read the global "
                    f"spatial shape from; got {sample_shape}.")
            wrapper = multi.MultiDiffusionModel2D(
                protocol_model, global_spatial_shape=tuple(sample_shape[1:]),
                condition_patch=True)
            wrapper.set_random_patching(
                tuple(int(n) for n in patching_settings["patch_shape"].GetVector()),
                patching_settings["patch_num"].GetInt()
                if patching_settings.Has("patch_num") else 1)
            loss_object = multi.MultiDiffusionMSEDSMLoss(wrapper, scheduler)
        else:
            loss_object = protocol.MSEDSMLoss(protocol_model, scheduler)

        def ComputeLoss(conditions, targets):
            # MSEDSMLoss reduces itself; the condition is a keyword, never
            # positional - the samplers never forward one on their own
            return loss_object(targets.to(device, dtype), condition=conditions.to(device, dtype))
    else:
        EDMLossSR, RegressionLoss, ResidualLoss, _ = _TryImportPhysicsNemoDiffusion()
        if loss_name == "edm_sr":
            legacy_loss = EDMLossSR(
                P_mean=p_mean, P_std=settings["P_std"].GetDouble(), sigma_data=sigma_data)
        elif loss_name == "regression":
            legacy_loss = RegressionLoss()
        else:  # residual
            if regression_model is None:
                raise ValueError(
                    "Loss \"residual\" needs the trained CorrDiff regression stage via "
                    "the regression_model argument.")
            regression_model = regression_model.to(device).eval()
            regression_model.requires_grad_(False)  # upstream does not no_grad it
            legacy_loss = ResidualLoss(
                regression_net=regression_model, P_mean=p_mean,
                P_std=settings["P_std"].GetDouble(), sigma_data=sigma_data)

        def ComputeLoss(conditions, targets):
            return legacy_loss(
                net=model, img_clean=targets.to(device, dtype),
                img_lr=conditions.to(device, dtype)).mean()

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
        for conditions, targets in loader:
            optimizer.zero_grad()
            loss = ComputeLoss(conditions, targets)
            loss.backward()
            optimizer.step()
            epoch_loss += loss.item()
            batches += 1
        history.append(epoch_loss / max(batches, 1))
        if echo_interval > 0 and (epoch + 1) % echo_interval == 0:
            Kratos.Logger.PrintInfo(
                "TrainDiffusionModel",
                f"epoch {epoch + 1}/{settings['epochs'].GetInt()}: loss = {history[-1]:.6e}")
    model.eval()
    return history


def TrainCorrDiffPair(regression_model, diffusion_model, dataset, settings: Kratos.Parameters):
    """The CorrDiff two-stage recipe in one call.

    Stage 1 trains the deterministic regression mean (loss "regression",
    "regression_epochs"/"regression_learning_rate" overriding the shared
    values when given); stage 2 trains the denoiser on the RESIDUAL
    target - regression(condition) with the frozen stage-1 model (loss
    "residual"). Inference then adds the regression mean back - see
    RunRegressionMean and DiffusionInferenceProcess's "regression_settings".

    Args:
        regression_model: e.g. CorrDiffRegressionUNet (trained in place).
        diffusion_model: the EDM-preconditioned denoiser (trained in place).
        dataset: (condition, target) pairs shared by both stages.
        settings: TrainDiffusionModel settings for the DIFFUSION stage,
            plus optional "regression_epochs" and
            "regression_learning_rate".

    Returns:
        (regression_history, diffusion_history)
    """
    settings = settings.Clone()
    regression_epochs = None
    if settings.Has("regression_epochs"):
        regression_epochs = settings["regression_epochs"].GetInt()
        settings.RemoveValue("regression_epochs")
    regression_learning_rate = None
    if settings.Has("regression_learning_rate"):
        regression_learning_rate = settings["regression_learning_rate"].GetDouble()
        settings.RemoveValue("regression_learning_rate")

    def with_overrides(loss, epochs=None, learning_rate=None):
        stage = settings.Clone()
        for key, value in (("loss", loss), ("epochs", epochs), ("learning_rate", learning_rate)):
            if value is None:
                continue
            if stage.Has(key):
                stage.RemoveValue(key)
            entry = stage.AddEmptyValue(key)
            if isinstance(value, str):
                entry.SetString(value)
            elif isinstance(value, int):
                entry.SetInt(value)
            else:
                entry.SetDouble(value)
        return stage

    regression_history = TrainDiffusionModel(
        regression_model, dataset,
        with_overrides("regression", regression_epochs, regression_learning_rate))
    diffusion_history = TrainDiffusionModel(
        diffusion_model, dataset, with_overrides("residual"),
        regression_model=regression_model)
    return regression_history, diffusion_history


def RunRegressionMean(model, condition_grid, output_channels: int = 0):
    """The CorrDiff regression stage's mean prediction for one condition.

    Args:
        model: The trained regression model ((x_zeros, img_lr) interface,
            e.g. CorrDiffRegressionUNet).
        condition_grid: (C_in, *spatial) array-like condition.
        output_channels: 0 reads the model's img_out_channels attribute.

    Returns:
        (C_out, *spatial) float64 numpy array.
    """
    torch = _TryImportTorch()

    if output_channels == 0:
        output_channels = getattr(model, "img_out_channels", 0)
        if not output_channels:
            raise ValueError(
                "output_channels is 0 and the model exposes no img_out_channels; "
                "set it explicitly.")

    parameter = next(model.parameters(), None)
    device = parameter.device if parameter is not None else torch.device("cpu")
    dtype = parameter.dtype if parameter is not None else torch.float32

    condition = torch.as_tensor(numpy.asarray(condition_grid)).to(device, dtype)[None]
    zeros = torch.zeros((1, output_channels) + tuple(condition.shape[2:]),
                        device=device, dtype=dtype)
    with torch.no_grad():
        mean = model(zeros, condition)
    return mean[0].cpu().to(torch.float64).numpy()


_GUIDANCE_TYPES = ("none", "data_consistency", "model_consistency")


def _BuildGuidance(settings: Kratos.Parameters, protocol, torch, scheduler,
                   observation, mask, observation_operator, device, dtype):
    """Builds a DPS guidance term from the "guidance" block, or None.

    Diffusion posterior sampling steers an ALREADY TRAINED denoiser toward
    measurements at sampling time: no retraining, and the constraint can be
    anything differentiable - sparse sensor readings ("data_consistency")
    or a forward operator such as the exact discrete FEM residual
    ("model_consistency").
    """
    default_settings = Kratos.Parameters("""{
        "type"  : "none",
        "std_y" : 0.1,
        "gamma" : 0.0,
        "norm"  : 2
    }""")
    settings.ValidateAndAssignDefaults(default_settings)
    guidance_type = settings["type"].GetString()
    if guidance_type not in _GUIDANCE_TYPES:
        raise ValueError(
            f"Unsupported guidance type \"{guidance_type}\". Supported: "
            f"{', '.join(_GUIDANCE_TYPES)}.")
    if guidance_type == "none":
        return None

    std_y = settings["std_y"].GetDouble()
    if std_y <= 0.0:
        raise ValueError(
            f"\"std_y\" is the observation noise level and must be > 0 [ std_y = "
            f"{std_y} ]. It also sets the guidance strength (the step is ~1/(2 std_y^2)), "
            "so a value far below the data's own scale makes the sampler diverge.")
    gamma = settings["gamma"].GetDouble()
    arguments = {
        "std_y": std_y,
        "gamma": gamma,
        "norm": settings["norm"].GetInt(),
        "alpha_fn": scheduler.alpha,
    }
    if gamma > 0.0:
        arguments["sigma_fn"] = scheduler.sigma  # upstream requires it for gamma > 0

    if guidance_type == "data_consistency":
        if observation is None or mask is None:
            raise ValueError(
                "Guidance type \"data_consistency\" needs both an observation and a "
                "mask (the process reads them from \"observation_fields\" and "
                "\"mask_field\"; a direct caller passes observation= and mask=).")
        y = torch.as_tensor(numpy.asarray(observation)).to(device, dtype)[None]
        mask_tensor = torch.as_tensor(numpy.asarray(mask)).to(device)[None]
        if mask_tensor.shape[1] == 1 and y.shape[1] > 1:
            mask_tensor = mask_tensor.expand_as(y)  # one mask for every channel
        if mask_tensor.shape != y.shape:
            raise ValueError(
                f"The observation mask has shape {tuple(mask_tensor.shape)[1:]} and the "
                f"observations {tuple(y.shape)[1:]}; they must match (or the mask carry "
                "one channel).")
        return protocol.DataConsistencyDPSGuidance(
            mask=mask_tensor.to(torch.bool), y=y, **arguments)

    if observation_operator is None or observation is None:
        raise ValueError(
            "Guidance type \"model_consistency\" needs an observation_operator "
            "A(x0) -> (B, *obs) and the matching observation y.")
    y = torch.as_tensor(numpy.asarray(observation)).to(device, dtype)
    if y.ndim == 0 or y.shape[0] != 1:
        y = y[None]  # the guidance broadcasts t against a batch-first y
    return protocol.ModelConsistencyDPSGuidance(
        observation_operator=observation_operator, y=y, **arguments)


def GenerateEnsemble(model, condition_grid, settings: Kratos.Parameters,
                     observation=None, mask=None, observation_operator=None):
    """Samples an ensemble of fields conditioned on one grid.

    Args:
        model: The trained denoiser, named by "denoiser_interface"
            (default "edm": the net(x, img_lr, sigma) contract).
        condition_grid: (C_in, *spatial) array-like condition (e.g. the
            coarse field sampled by grid_bridge, thin axis squeezed).
        settings: Kratos Parameters; defaults:
            num_samples (8), num_steps (18), solver ("heun"),
            solver_options ({}), output_channels (0 = read the model's
            img_out_channels), seed (-1 = leave the RNG alone),
            api ("auto"), denoiser_interface ("edm"), preconditioner
            ("none"), sigma_min (0.002), sigma_max (80.0), rho (7.0),
            sigma_data (0.5), guidance ({}), patching ({}).

            On the protocol API "solver" is one of "heun", "euler",
            "edm_stochastic_euler" and "edm_stochastic_heun"
            ("solver_options" reaches the solver's constructor; the
            stochastic ones' S_churn lives there). "api" : "legacy" runs
            the deprecated deterministic_sampler instead and accepts
            neither guidance nor patching.
        observation: Measurements in the model's own (normalized) output
            space, for DPS guidance. (C_out, *spatial) for
            "data_consistency", (*obs) for "model_consistency".
        mask: (C_out, *spatial) or (1, *spatial) 0/1 array marking the
            observed entries, for "data_consistency".
        observation_operator: A differentiable A(x0) -> (B, *obs) for
            "model_consistency" - e.g. the exact discrete FEM residual of
            physics.diffusion_residual_operator.

    Returns:
        (num_samples, C_out, *spatial) float64 numpy array.
    """
    torch = _TryImportTorch()

    default_settings = Kratos.Parameters("""{
        "num_samples"        : 8,
        "num_steps"          : 18,
        "solver"             : "heun",
        "solver_options"     : {},
        "output_channels"    : 0,
        "seed"               : -1,
        "api"                : "auto",
        "denoiser_interface" : "edm",
        "preconditioner"     : "none",
        "sigma_min"          : 0.002,
        "sigma_max"          : 80.0,
        "rho"                : 7.0,
        "sigma_data"         : 0.5,
        "guidance"           : {},
        "patching"           : {}
    }""")
    settings.ValidateAndAssignDefaults(default_settings)

    num_samples = settings["num_samples"].GetInt()
    if num_samples < 1:
        raise ValueError(f"\"num_samples\" must be >= 1 [ num_samples = {num_samples} ].")
    num_steps = settings["num_steps"].GetInt()
    output_channels = settings["output_channels"].GetInt()
    if output_channels == 0:
        output_channels = getattr(model, "img_out_channels", 0)
        if not output_channels:
            raise ValueError(
                "\"output_channels\" is 0 and the model exposes no img_out_channels; "
                "set it explicitly.")

    api = _ResolveDiffusionApi(settings["api"].GetString(), "edm_sr")
    guidance_settings = settings["guidance"].Clone()
    guided = (guidance_settings.Has("type")
              and guidance_settings["type"].GetString() != "none")
    patching_settings = settings["patching"].Clone()
    patched = (patching_settings.Has("patch_shape")
               and patching_settings["patch_shape"].size() > 0)
    if api == "legacy" and (guided or patched):
        raise ValueError(
            "DPS guidance and multi-diffusion patching exist only on physicsnemo's "
            "protocol API; they cannot run with \"api\" : \"legacy\".")

    seed = settings["seed"].GetInt()
    if seed >= 0:
        torch.manual_seed(seed)

    parameter = next(model.parameters(), None)
    device = parameter.device if parameter is not None else torch.device("cpu")
    dtype = parameter.dtype if parameter is not None else torch.float32

    condition = torch.as_tensor(numpy.asarray(condition_grid)).to(device, dtype)[None]
    spatial = tuple(condition.shape[2:])

    if api == "legacy":
        _, _, _, deterministic_sampler = _TryImportPhysicsNemoDiffusion()
        samples = []
        with torch.no_grad():
            for _ in range(num_samples):
                latents = torch.randn((1, output_channels) + spatial, device=device, dtype=dtype)
                sample = deterministic_sampler(
                    net=model, latents=latents, img_lr=condition,
                    num_steps=num_steps, solver=settings["solver"].GetString())
                samples.append(sample[0].cpu().to(torch.float64).numpy())
        return numpy.stack(samples)

    if num_steps < 2:
        raise ValueError(
            "\"num_steps\" must be >= 2: the noise schedule spans num_steps - 1 "
            f"intervals, so a single step yields NaN timesteps [ num_steps = {num_steps} ].")

    protocol = _TryImportPhysicsNemoDiffusionProtocol()
    sigma_data = settings["sigma_data"].GetDouble()
    scheduler = protocol.EDMNoiseScheduler(
        sigma_min=settings["sigma_min"].GetDouble(),
        sigma_max=settings["sigma_max"].GetDouble(),
        rho=settings["rho"].GetDouble(), sigma_data=sigma_data)
    protocol_model = _AsProtocolModel(
        model, settings["denoiser_interface"].GetString(), out_channels=output_channels,
        preconditioner=settings["preconditioner"].GetString(), sigma_data=sigma_data)

    solver_options = settings["solver_options"].Clone()
    options = {key: _ReadParameterValue(solver_options[key]) for key in solver_options.keys()}
    if settings["solver"].GetString().startswith("edm_stochastic"):
        # the stochastic solvers scale their churn by their OWN num_steps,
        # which defaults to 18 whatever the sampler is asked for
        options.setdefault("num_steps", num_steps)

    if patched:
        if len(spatial) != 2:
            raise ValueError(
                "Multi-diffusion patching is 2-D upstream; a "
                f"{len(spatial)}-D latent stays whole-grid (the volumetric \"unet3d\" "
                "path always does).")
        multi = _TryImportPhysicsNemoMultiDiffusion()
        patch_shape = tuple(int(n) for n in patching_settings["patch_shape"].GetVector())
        overlap_pixels = (patching_settings["overlap_pix"].GetInt()
                          if patching_settings.Has("overlap_pix") else 0)
        boundary_pixels = (patching_settings["boundary_pix"].GetInt()
                           if patching_settings.Has("boundary_pix") else 0)
        chunk_size = (patching_settings["chunk_size"].GetInt()
                      if patching_settings.Has("chunk_size") else 0)

        def MakePredictor():
            # a MultiDiffusionPredictor MUTATES the wrapper it is given
            # (it turns fusion off inside it), so build a fresh one per call
            wrapper = multi.MultiDiffusionModel2D(
                protocol_model, global_spatial_shape=spatial, condition_patch=True)
            predictor = multi.MultiDiffusionPredictor(
                wrapper, condition=condition, chunk_size=chunk_size or None)
            predictor.set_patching(
                overlap_pixels, boundary_pixels,
                patch_shape=patch_shape, global_shape=spatial)
            return predictor
    else:
        def MakePredictor():
            # samplers.sample never forwards a condition: bind it here
            return functools.partial(protocol_model, condition=condition)

    guidance = _BuildGuidance(
        guidance_settings, protocol, torch, scheduler,
        observation, mask, observation_operator, device, dtype)

    samples = []
    with torch.no_grad():  # the DPS predictors re-enable grad internally
        for _ in range(num_samples):
            time_steps = scheduler.timesteps(num_steps, device=device, dtype=dtype)
            latents = scheduler.init_latents(
                (output_channels,) + spatial, time_steps[0].expand(1),
                device=device, dtype=dtype)
            predictor = MakePredictor()
            if guidance is None:
                denoiser = scheduler.get_denoiser(x0_predictor=predictor)
            else:
                denoiser = scheduler.get_denoiser(
                    score_predictor=protocol.DPSScorePredictor(
                        predictor, scheduler.x0_to_score, guidance))
            sample = protocol.sample(
                denoiser, latents, scheduler, num_steps,
                solver=settings["solver"].GetString(), solver_options=options or None)
            if not bool(torch.isfinite(sample).all()):
                raise ValueError(
                    "The sampler produced non-finite values. "
                    + ("The guidance step scales as 1/(2 std_y^2), so a \"std_y\" far "
                       "below the observations' own scale overshoots and diverges - "
                       "raise it." if guidance is not None else
                       "Check the model's scaling against \"sigma_data\"."))
            samples.append(sample[0].cpu().to(torch.float64).numpy())
    return numpy.stack(samples)


def _ReadParameterValue(entry: Kratos.Parameters):
    """A scalar Parameters entry as its python value."""
    if entry.IsBool():
        return entry.GetBool()
    if entry.IsInt():
        return entry.GetInt()
    if entry.IsNumber():
        return entry.GetDouble()
    if entry.IsString():
        return entry.GetString()
    if entry.IsVector():
        return list(entry.GetVector())
    raise ValueError(
        "\"solver_options\" entries must be scalars, strings or vectors; "
        f"got {entry.PrettyPrintJsonString()}.")
