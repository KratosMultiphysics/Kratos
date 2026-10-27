"""Self-supervised geometry pretraining with AeroJEPA.

physicsnemo.experimental.models.aerojepa.AeroJEPA learns from point sets:
a CONTEXT cloud describing a shape, and QUERY points at which something
about that shape is predicted. The appeal for this application is that the
context cloud costs nothing - it is the surface of any geometry the mesh
bridge can produce - so a model can be pretrained across a family of shapes
BEFORE a single solve exists, and a supervised head fitted afterwards on
the few solves one can afford.

WHAT UPSTREAM DOES NOT SHIP. The JEPA recipe the architecture is named for
- predict the target encoder's tokens, keep the target encoder as an
exponential moving average of the context encoder - has no loss and no EMA
in physicsnemo 2.2; only the forward pass exists. Rather than guess at the
token-level objective, the pretext task here is one the public forward
supports directly: predict the SIGNED DISTANCE at query points from the
context cloud. It is self-supervised in the sense that matters (the labels
come from the geometry, not from a solver), and what transfers is the same
trunk.

torch and physicsnemo are imported lazily.
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
            "PhysicsNeMoApplication.aerojepa_pretraining requires torch, which could not "
            "be imported. Install it with e.g. 'pip install torch'.") from e


def _TryImportAeroJepa():
    import warnings

    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")  # physicsnemo.experimental warns on import
            from physicsnemo.experimental.models.aerojepa import (
                AeroJEPA, AeroJEPATrunk, ContextTransformer, QueryTokenDecoder,
                PrototypeTokenJEPAHead, TargetTransformer)
        return {
            "AeroJEPA": AeroJEPA, "AeroJEPATrunk": AeroJEPATrunk,
            "ContextTransformer": ContextTransformer,
            "TargetTransformer": TargetTransformer,
            "QueryTokenDecoder": QueryTokenDecoder,
            "PrototypeTokenJEPAHead": PrototypeTokenJEPAHead,
        }
    except ImportError as e:
        raise ImportError(
            "AeroJEPA requires physicsnemo >= 2.2 (its experimental models), which could "
            "not be imported. Install it with e.g. "
            "'pip install -U nvidia-physicsnemo'.") from e


def CreateAeroJepaModel(settings: Kratos.Parameters):
    """Builds a small AeroJEPA.

    Args:
        settings: Kratos Parameters; defaults:
            point_input_dim (3), token_dim (32), max_point_tokens (16),
            tokenizer_cluster_size (8), tokenizer_voxel_size (0.25),
            num_heads (2), num_layers (1), neighbor_k (8),
            decoder_hidden_dim (64), decoder_layers (2), out_dim (1),
            condition_dim (1), use_transformer_engine (false).

    Returns:
        An AeroJEPA.

    TRAP: the point tokenizer's default strategy is "voxel_fps_cluster" and
    it REQUIRES a voxel size - without one the forward raises deep inside
    the tokenizer, long after the model was built, so the size is a settings
    key here rather than an optional extra.
    """
    torch = _TryImportTorch()
    parts = _TryImportAeroJepa()

    settings.ValidateAndAssignDefaults(Kratos.Parameters("""{
        "point_input_dim"          : 3,
        "token_dim"                : 32,
        "max_point_tokens"         : 16,
        "tokenizer_cluster_size"   : 8,
        "tokenizer_voxel_size"     : 0.25,
        "num_heads"                : 2,
        "num_layers"               : 1,
        "neighbor_k"               : 8,
        "decoder_hidden_dim"       : 64,
        "decoder_layers"           : 2,
        "out_dim"                  : 1,
        "condition_dim"            : 1,
        "use_transformer_engine"   : false
    }"""))
    token_dim = settings["token_dim"].GetInt()
    use_te = settings["use_transformer_engine"].GetBool()
    encoder_arguments = dict(
        point_input_dim=settings["point_input_dim"].GetInt(),
        token_dim=token_dim,
        max_point_tokens=settings["max_point_tokens"].GetInt(),
        tokenizer_cluster_size=settings["tokenizer_cluster_size"].GetInt(),
        tokenizer_voxel_size=settings["tokenizer_voxel_size"].GetDouble(),
        num_heads=settings["num_heads"].GetInt(),
        num_layers=settings["num_layers"].GetInt(),
        neighbor_k=settings["neighbor_k"].GetInt(),
        use_te=use_te)

    trunk = parts["AeroJEPATrunk"](
        context_encoder=parts["ContextTransformer"](**encoder_arguments),
        target_encoder=parts["TargetTransformer"](**encoder_arguments),
        decoder=parts["QueryTokenDecoder"](
            token_dim=token_dim, query_dim=3,
            hidden_dim=settings["decoder_hidden_dim"].GetInt(),
            num_layers=settings["decoder_layers"].GetInt(),
            out_dim=settings["out_dim"].GetInt(), use_sdf=True,
            cross_attention_heads=settings["num_heads"].GetInt(),
            cross_attention_layers=1,
            cross_attention_k=settings["neighbor_k"].GetInt(), use_te=use_te))
    return parts["AeroJEPA"](
        trunk=trunk,
        predictor=parts["PrototypeTokenJEPAHead"](
            token_dim=token_dim, cond_dim=settings["condition_dim"].GetInt(),
            hidden_dim=settings["decoder_hidden_dim"].GetInt(), depth=1,
            num_heads=settings["num_heads"].GetInt(),
            neighbor_k=settings["neighbor_k"].GetInt(), use_te=use_te))


def CreateGeometryFamilyDataset(count: int = 8, seed: int = 0, n_context: int = 256,
                                n_query: int = 64, resolution: int = 24):
    """A family of implicit shapes as (context, query) point sets.

    Each sample is a sphere, a box or their union at random parameters,
    surfaced through mesh_bridge.generate and sampled by
    mesh_bridge.spatial.ComputeSignedDistance - so both the context cloud
    and its labels come from the geometry alone.

    Returns:
        A list of dicts with context_pos, context_feat, gen_params,
        query_pos, query_sdf and a "kind" label (for probing what the
        embedding learned).
    """
    torch = _TryImportTorch()
    from KratosMultiphysics.PhysicsNeMoApplication.bridges.mesh_bridge import (
        generate, spatial)

    primitives = generate.SdfPrimitives()
    generator = numpy.random.default_rng(seed)
    low, high = (-0.6, -0.6, -0.6), (0.6, 0.6, 0.6)
    axes = [torch.linspace(low[i], high[i], resolution) for i in range(3)]
    lattice = torch.stack(torch.meshgrid(*axes, indexing="ij"), dim=-1).reshape(-1, 3)

    samples = []
    for index in range(count):
        kind = index % 2
        radius = float(generator.uniform(0.25, 0.4))
        if kind == 0:
            phi = primitives["sphere"]((0.0, 0.0, 0.0), radius)
        else:
            half = float(generator.uniform(0.2, 0.35))
            phi = primitives["box"]((-half, -half, -half), (half, half, half))
        # SurfaceFromLevelSet marches a SAMPLED field, not the callable
        field = phi(lattice).reshape(resolution, resolution, resolution)
        mesh = generate.SurfaceFromLevelSet(field, bounding_box=(low, high))

        points = torch.as_tensor(mesh.points, dtype=torch.float32)
        chosen = generator.choice(points.shape[0],
                                  size=min(n_context, points.shape[0]), replace=False)
        context = points[torch.from_numpy(chosen)]
        query = torch.as_tensor(
            generator.uniform(-0.6, 0.6, size=(n_query, 3)), dtype=torch.float32)
        distances = torch.as_tensor(
            numpy.asarray(spatial.ComputeSignedDistance(mesh, query.numpy())),
            dtype=torch.float32).reshape(-1, 1)
        samples.append({
            "context_pos": context,
            "context_feat": context.clone(),
            "gen_params": torch.tensor([radius], dtype=torch.float32),
            "query_pos": query,
            "query_sdf": distances,
            "kind": kind,
        })
    return samples


def PretrainOnGeometry(model, samples, settings: Kratos.Parameters):
    """Fits the model to predict signed distance from the context cloud.

    The pretext task, and deliberately not a solver quantity: the labels are
    the geometry's own, so pretraining costs no solves at all.

    Args:
        model: A CreateAeroJepaModel result.
        samples: CreateGeometryFamilyDataset output.
        settings: Kratos Parameters; defaults: epochs (20),
            learning_rate (1e-3), device ("cpu"), shuffle (true), seed (-1),
            echo_interval (0).

    Returns:
        list[float]: mean loss per epoch.
    """
    torch = _TryImportTorch()

    settings.ValidateAndAssignDefaults(Kratos.Parameters("""{
        "epochs"        : 20,
        "learning_rate" : 1e-3,
        "device"        : "cpu",
        "shuffle"       : true,
        "seed"          : -1,
        "echo_interval" : 0
    }"""))
    if not samples:
        raise ValueError("PretrainOnGeometry needs at least one sample.")
    seed = settings["seed"].GetInt()
    if seed >= 0:
        torch.manual_seed(seed)
    device = model_registry.ResolveDevice(settings["device"].GetString())
    model = model.to(device)
    optimizer = torch.optim.Adam(model.parameters(),
                                 lr=settings["learning_rate"].GetDouble())
    generator = numpy.random.default_rng(None if seed < 0 else seed)
    echo_interval = settings["echo_interval"].GetInt()

    history = []
    model.train()
    for epoch in range(settings["epochs"].GetInt()):
        order = (generator.permutation(len(samples)) if settings["shuffle"].GetBool()
                 else range(len(samples)))
        epoch_loss = 0.0
        for index in order:
            sample = samples[index]
            optimizer.zero_grad()
            predicted = model(
                context_pos=sample["context_pos"].to(device),
                context_feat=sample["context_feat"].to(device),
                gen_params=sample["gen_params"].to(device),
                query_pos=sample["query_pos"].to(device),
                query_sdf=torch.zeros_like(sample["query_sdf"]).to(device))
            loss = torch.nn.functional.mse_loss(
                predicted, sample["query_sdf"].to(device))
            loss.backward()
            optimizer.step()
            epoch_loss += loss.item()
        history.append(epoch_loss / len(samples))
        if echo_interval > 0 and (epoch + 1) % echo_interval == 0:
            Kratos.Logger.PrintInfo(
                "PretrainOnGeometry",
                f"epoch {epoch + 1}: loss = {history[-1]:.6e}")
    model.eval()
    return history


def EncodeGeometry(model, sample, device=None):
    """The trunk's context embedding for one geometry, mean-pooled.

    What pretraining is FOR: a fixed-width descriptor of a shape that a
    downstream head - or a linear probe, as the tests use - can consume.

    TRAP: ``encode_context`` is keyword-only and wants the TARGET point sets
    as well as the context, because the architecture was written for a
    context/target JEPA pairing. With no separate target to offer, the
    context cloud is passed as the target surface and the query points as
    the target volume, which is what the upstream forward does internally
    for a single geometry.
    """
    torch = _TryImportTorch()
    device = device or next(model.parameters()).device
    context_pos = sample["context_pos"].to(device)
    context_feat = sample["context_feat"].to(device)
    query_pos = sample["query_pos"].to(device)
    with torch.no_grad():
        encoded = model.trunk.encode_context(
            context_pos=context_pos, context_feat=context_feat,
            target_surface_pos=context_pos, target_surface_main_feat=context_feat,
            target_volume_pos=query_pos,
            # the target encoder concatenates the surface and volume feature
            # blocks and refuses mismatched widths, so the volume features are
            # zeros of the CONTEXT's width rather than a single column
            target_volume_feat=torch.zeros(
                query_pos.shape[0], context_feat.shape[-1],
                dtype=context_feat.dtype, device=device),
            gen_params=sample["gen_params"].to(device))

    tokens = encoded
    if isinstance(encoded, dict):
        candidates = [value for value in encoded.values()
                      if torch.is_tensor(value) and value.ndim >= 2]
        if not candidates:
            raise ValueError(
                f"encode_context returned no token tensor; it gave {sorted(encoded)}.")
        tokens = candidates[0]
    if isinstance(tokens, (tuple, list)):
        tokens = tokens[0]
    return tokens.reshape(-1, tokens.shape[-1]).mean(dim=0)
