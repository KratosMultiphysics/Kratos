"""torch Dataset over DatasetExportProcess output directories.

Lets users train directly on Kratos-exported .npz samples:

    dataset = CreateNpzDataset("my_dataset",
                               input_keys=["VELOCITY__node_historical"],
                               output_keys=["PRESSURE__node_historical"])
    loader = torch.utils.data.DataLoader(dataset, batch_size=4)

torch is an optional runtime dependency, imported lazily inside the factory.
"""

from pathlib import Path

import numpy


def _TryImportTorch():
    try:
        import torch
        return torch
    except ImportError as e:
        raise ImportError(
            "PhysicsNeMoApplication.torch_dataset requires torch, which could not be "
            "imported. Install it with e.g. 'pip install torch'.") from e


def CreateNpzDataset(directory, input_keys, output_keys):
    """Creates a torch Dataset over a DatasetExportProcess output directory.

    One item per .npz file (sorted by their step suffix). Each item is a
    tuple (inputs, outputs): the named arrays flattened per entity to
    (n_entities, width) float32 tensors and concatenated along the last axis
    — the same layout InferenceProcess feeds to and expects from a model.

    Args:
        directory: The directory holding "<prefix>_<step>.npz" files.
        input_keys: Field keys ("<VARIABLE>__<location>") forming the input.
        output_keys: Field keys forming the training target.

    Returns:
        A torch.utils.data.Dataset instance.
    """
    torch = _TryImportTorch()

    directory = Path(directory)
    files = sorted(directory.glob("*.npz"), key=lambda p: int(p.stem.rsplit("_", 1)[-1]))
    if not files:
        raise FileNotFoundError(f"No .npz sample files found in \"{directory}\".")

    def assemble(data, keys):
        arrays = []
        for key in keys:
            if key not in data:
                raise KeyError(
                    f"Field \"{key}\" not found in sample (available: {sorted(data.files)}).")
            array = numpy.asarray(data[key], dtype=numpy.float32)
            arrays.append(array.reshape(array.shape[0], -1))
        return torch.from_numpy(numpy.concatenate(arrays, axis=-1))

    class KratosNpzDataset(torch.utils.data.Dataset):
        def __len__(self):
            return len(files)

        def __getitem__(self, index):
            with numpy.load(files[index]) as data:
                return assemble(data, input_keys), assemble(data, output_keys)

    return KratosNpzDataset()


def _SortedGridFiles(directory):
    directory = Path(directory)
    files = sorted(directory.glob("*.npz"), key=lambda p: int(p.stem.rsplit("_", 1)[-1]))
    if not files:
        raise FileNotFoundError(f"No .npz grid files found in \"{directory}\".")
    return files


def _LoadGrid(path, squeeze_axis):
    with numpy.load(path) as data:
        grid = numpy.asarray(data["grid"], dtype=numpy.float32)  # (C, D, H, W)
    if squeeze_axis is not None:
        grid = grid.mean(axis=1 + squeeze_axis)  # thin-axis idiom: collapse it
    return grid


def _CheckSqueezeAxis(squeeze_axis):
    if squeeze_axis is not None and squeeze_axis not in (0, 1, 2):
        raise ValueError(f"squeeze_axis must be None, 0, 1 or 2, got {squeeze_axis}.")
    return squeeze_axis


def CreateGridSequenceDataset(directory, nr_tsteps, squeeze_axis=None):
    """Creates a torch Dataset of grid sequences over a GridDatasetExportProcess
    output directory, in the physicsnemo RNN layout.

    Item i pairs the state at step i with the nr_tsteps following states:
    (x0, y) with x0 of shape (C, 1, *spatial) and y of shape
    (C, nr_tsteps, *spatial) - exactly what One2ManyRNN consumes and
    produces ((N, C, 1, ...) -> (N, C, T, ...)), so TrainModel works
    unchanged.

    Args:
        directory: The GridDatasetExportProcess output directory (grids
            sorted by their step suffix; all on the same lattice).
        nr_tsteps: The number of future states per item (the model's
            nr_tsteps).
        squeeze_axis: None (default) keeps the (C, D, H, W) grids for
            dimension=3 models; 0, 1 or 2 collapses that spatial axis by its
            mean for dimension=2 models - the thin-axis idiom for planar
            cases (e.g. grid_shape [16, 16, 2] with squeeze_axis=2).

    Returns:
        A torch.utils.data.Dataset instance.
    """
    torch = _TryImportTorch()
    files = _SortedGridFiles(directory)
    squeeze_axis = _CheckSqueezeAxis(squeeze_axis)
    nr_tsteps = int(nr_tsteps)
    if nr_tsteps < 1:
        raise ValueError(f"nr_tsteps must be >= 1, got {nr_tsteps}.")
    if len(files) <= nr_tsteps:
        raise ValueError(
            f"Need more than nr_tsteps = {nr_tsteps} exported grids, found {len(files)}.")

    class KratosGridSequenceDataset(torch.utils.data.Dataset):
        def __len__(self):
            return len(files) - nr_tsteps

        def __getitem__(self, index):
            initial = _LoadGrid(files[index], squeeze_axis)
            future = numpy.stack(
                [_LoadGrid(files[index + 1 + t], squeeze_axis) for t in range(nr_tsteps)],
                axis=1)  # (C, T, *spatial)
            return torch.from_numpy(initial[:, None]), torch.from_numpy(future)

    return KratosGridSequenceDataset()


def CreateGridPairDataset(input_directory, target_directory, squeeze_axis=None):
    """Creates a torch Dataset of (input_grid, target_grid) pairs from two
    GridDatasetExportProcess output directories matched by step suffix.

    The layout for conditional grid-to-grid training (diffusion downscaling,
    superresolution): item i is (x, y) with x from input_directory (e.g. the
    coarse/condition grid) and y from target_directory (e.g. the fine
    grid), both (C, *spatial) float32.

    Args:
        input_directory / target_directory: The two export directories; they
            must contain the same steps.
        squeeze_axis: See CreateGridSequenceDataset.

    Returns:
        A torch.utils.data.Dataset instance.
    """
    torch = _TryImportTorch()
    input_files = _SortedGridFiles(input_directory)
    target_files = _SortedGridFiles(target_directory)
    squeeze_axis = _CheckSqueezeAxis(squeeze_axis)

    def steps(files):
        return [int(p.stem.rsplit("_", 1)[-1]) for p in files]

    if steps(input_files) != steps(target_files):
        raise ValueError(
            f"Step mismatch between \"{input_directory}\" ({steps(input_files)}) and "
            f"\"{target_directory}\" ({steps(target_files)}).")

    class KratosGridPairDataset(torch.utils.data.Dataset):
        def __len__(self):
            return len(input_files)

        def __getitem__(self, index):
            return (torch.from_numpy(_LoadGrid(input_files[index], squeeze_axis)),
                    torch.from_numpy(_LoadGrid(target_files[index], squeeze_axis)))

    return KratosGridPairDataset()


# Minimal key sets the CAE datapipes need per model type (verified against
# physicsnemo 2.1.1; CaeDatasetExportProcess writes the superset of all of
# them, and CAEDataset ignores keys it is not asked to read).
DOMINO_COMMON_KEYS = ("stl_coordinates", "stl_faces", "stl_centers", "stl_areas",
                      "global_params_values", "global_params_reference")
DOMINO_SURFACE_KEYS = DOMINO_COMMON_KEYS + ("surface_mesh_centers", "surface_normals",
                                            "surface_areas", "surface_fields")
DOMINO_VOLUME_KEYS = DOMINO_COMMON_KEYS + ("volume_mesh_centers", "volume_fields")
TRANSOLVER_SURFACE_KEYS = ("stl_centers", "surface_mesh_centers", "surface_normals",
                           "surface_areas", "surface_fields", "stream_velocity", "air_density")
TRANSOLVER_VOLUME_KEYS = ("stl_centers", "stl_coordinates", "stl_faces",
                          "volume_mesh_centers", "volume_fields", "stream_velocity", "air_density")


def _TryImportCae():
    try:
        from physicsnemo.datapipes.cae.cae_dataset import CAEDataset
        from physicsnemo.datapipes.cae.domino_datapipe import DoMINODataPipe
        from physicsnemo.datapipes.cae.transolver_datapipe import TransolverDataPipe
        return CAEDataset, DoMINODataPipe, TransolverDataPipe
    except ImportError as e:
        raise ImportError(
            "PhysicsNeMoApplication.torch_dataset's CAE datapipe factories require "
            "physicsnemo, which could not be imported. Install it with e.g. "
            "'pip install nvidia-physicsnemo'.") from e


def CreateCaeDataset(directory, keys_to_read, keys_to_read_if_available=None, device=None, **kwargs):
    """Creates a physicsnemo CAEDataset over a CaeDatasetExportProcess output.

    One item per exported "<prefix>_<id>.npz" case; each item is a dict of
    torch tensors holding exactly keys_to_read (extra keys in the files are
    ignored - the exporter writes a superset on purpose).

    Args:
        directory: The CaeDatasetExportProcess output directory.
        keys_to_read: The keys to load (e.g. one of the *_KEYS constants).
        keys_to_read_if_available: Optional {key: default_tensor} extras.
        device: Output device (default "cpu").
        **kwargs: Forwarded to CAEDataset (preload_depth, pin_memory, ...).
    """
    CAEDataset, _, _ = _TryImportCae()
    import torch

    return CAEDataset(
        data_dir=str(directory),
        keys_to_read=list(keys_to_read),
        keys_to_read_if_available=keys_to_read_if_available or {},
        output_device=torch.device(device or "cpu"),
        **kwargs)


class _BoundingBox:
    """The .min/.max attribute object DoMINO's config expects."""

    def __init__(self, low, high):
        self.min = [float(v) for v in low]
        self.max = [float(v) for v in high]


def CreateDoMINODataPipe(directory, model_type, bounding_box, bounding_box_surface=None,
                         keys_to_read=None, phase="train", device=None, **overrides):
    """Creates a DoMINODataPipe over a CaeDatasetExportProcess output.

    Args:
        directory: The export directory (per-case .npz files).
        model_type: "surface", "volume" or "combined".
        bounding_box: ((xmin, ymin, zmin), (xmax, ymax, zmax)) of the volume
            domain - REQUIRED by DoMINO's preprocessing even for
            surface-only models.
        bounding_box_surface: Optional surface bounding box; defaults to
            bounding_box.
        keys_to_read: Defaults to the model_type-appropriate *_KEYS constant
            ("combined" = union of surface + volume sets).
        phase: "train" (default), "val" or "test".
        device: None/"cpu" keeps preprocessing and output on CPU (safe
            default); a CUDA device enables GPU preprocessing.
        **overrides: Forwarded to DoMINODataConfig (grid_resolution,
            sampling, num_surface_neighbors, ...).

    Returns:
        The DoMINODataPipe, dataset already wired (len() / indexable).
    """
    _, DoMINODataPipe, _ = _TryImportCae()
    import torch

    low, high = bounding_box
    surface_low, surface_high = bounding_box_surface or bounding_box
    use_gpu = torch.device(device).type == "cuda" if device is not None else False
    settings = {
        "phase": phase,
        "bounding_box_dims": _BoundingBox(low, high),
        "bounding_box_dims_surf": _BoundingBox(surface_low, surface_high),
        "gpu_preprocessing": use_gpu,
        "gpu_output": use_gpu,
        "sampling": False,
    }
    settings.update(overrides)
    pipe = DoMINODataPipe(input_path=str(directory), model_type=model_type, **settings)

    if keys_to_read is None:
        keys_to_read = {
            "surface": DOMINO_SURFACE_KEYS,
            "volume": DOMINO_VOLUME_KEYS,
            "combined": tuple(dict.fromkeys(DOMINO_SURFACE_KEYS + DOMINO_VOLUME_KEYS)),
        }[model_type]
    pipe.set_dataset(CreateCaeDataset(directory, keys_to_read, device=device))
    return pipe


def CreateTransolverDataPipe(directory, model_type, keys_to_read=None,
                             resolution=None, device=None, **overrides):
    """Creates a TransolverDataPipe over a CaeDatasetExportProcess output.

    Two npz-reader necessities are handled here (both verified):
    volume_sample_from_disk is forced off (the npz reader cannot subsample
    on disk), and resolution defaults to None (use all points) instead of
    Transolver's 200k default, which breaks on smaller meshes.

    Args:
        directory: The export directory (per-case .npz files).
        model_type: "surface" or "volume".
        keys_to_read: Defaults to the model_type-appropriate *_KEYS constant.
        resolution: Optional number of points to sample per case.
        device: Output device (default "cpu").
        **overrides: Forwarded to TransolverDataConfig (include_normals,
            include_sdf, scaling_type, ...).

    Returns:
        The TransolverDataPipe, dataset already wired (len() / indexable).
    """
    _, _, TransolverDataPipe = _TryImportCae()

    settings = {
        "resolution": resolution,
        "volume_sample_from_disk": False,
    }
    settings.update(overrides)
    pipe = TransolverDataPipe(input_path=str(directory), model_type=model_type, **settings)

    if keys_to_read is None:
        keys_to_read = TRANSOLVER_SURFACE_KEYS if model_type == "surface" else TRANSOLVER_VOLUME_KEYS
    pipe.set_dataset(CreateCaeDataset(directory, keys_to_read, device=device))
    return pipe


def CreateMeshDataset(directory, transforms=None, device=None, num_workers=1):
    """Creates a physicsnemo MeshDataset over a MeshExportProcess output.

    One item per saved "*.pmsh" mesh; each item is what physicsnemo's
    MeshReader yields — a (Mesh, metadata) tuple with the exported fields in
    the Mesh's point_data/cell_data. This is the mesh-native counterpart of
    CreateNpzDataset, for training mesh-based models (transforms and device
    streaming are handled by physicsnemo's own datapipe).

    Args:
        directory: The MeshExportProcess output directory.
        transforms: Optional sequence of physicsnemo MeshTransform instances.
        device: Optional device the meshes are moved to.
        num_workers: Reader worker count.

    Returns:
        A physicsnemo.datapipes.mesh_dataset.MeshDataset instance.
    """
    try:
        from physicsnemo.datapipes.mesh_dataset import MeshDataset, MeshReader
    except ImportError as e:
        raise ImportError(
            "PhysicsNeMoApplication.torch_dataset.CreateMeshDataset requires physicsnemo, "
            "which could not be imported. Install it with e.g. 'pip install nvidia-physicsnemo'.") from e

    # MeshReader itself raises a clear ValueError when the directory holds no
    # "*.pmsh" meshes.
    reader = MeshReader(str(directory))
    return MeshDataset(reader, transforms=transforms, device=device, num_workers=num_workers)


def CreateParticleTrajectoryDataset(trajectories, history_size, delta_time,
                                    normalize=False):
    """Windows position trajectories into Learning-to-Simulate samples.

    From each (T, N, 3) position trajectory (numpy array-like; e.g. loaded
    from DatasetExportProcess .npz series), finite-difference velocities
    v_t = (x_t - x_{t-1}) / dt are formed and every valid step t yields one
    (features, target) pair:

        features (N, K*3): the last K velocities, oldest first - exactly
            particle_bridge.BuildKinematicFeatures's layout;
        target (N, 3): the central-difference acceleration
            a_t = (x_{t+1} - 2 x_t + x_{t-1}) / dt**2.

    The per-sample particle positions (for graph building) are exposed as
    dataset.positions[i] (numpy (N, 3)). Normalization statistics over all
    samples are computed either way (dataset.feature_mean/std and
    dataset.target_mean/std, per channel); normalize=True bakes them into
    the returned tensors - FEATURES included, so a model trained this way
    expects standardized velocity histories too. Write BOTH halves into the
    model card (MakeNormalizationCardEntries: the feature statistics as
    "input_normalization", the target statistics as "output_normalization");
    ParticleInferenceProcess reads both, standardizing the history it
    gathers and de-normalizing the acceleration it integrates.

    Args:
        trajectories: One (T, N, 3) array or an iterable of them.
        history_size: K >= 1 velocity states per feature window.
        delta_time: The (constant) time-step of the trajectories.
        normalize: Standardize features and targets with the dataset stats.

    Returns:
        A torch Dataset of (features, target) float64 tensor pairs.
    """
    torch = _TryImportTorch()

    if history_size < 1:
        raise ValueError(f"history_size must be >= 1 [ history_size = {history_size} ].")
    if delta_time <= 0.0:
        raise ValueError(f"delta_time must be > 0 [ delta_time = {delta_time} ].")

    trajectories = numpy.asarray(trajectories, dtype=float) if not isinstance(
        trajectories, (list, tuple)) else [numpy.asarray(t, dtype=float) for t in trajectories]
    if not isinstance(trajectories, list):
        trajectories = [trajectories]

    features_list, targets_list, positions_list = [], [], []
    for trajectory in trajectories:
        if trajectory.ndim != 3 or trajectory.shape[2] != 3:
            raise ValueError(
                f"Each trajectory must have shape (T, N, 3); got {trajectory.shape}.")
        steps = trajectory.shape[0]
        if steps < history_size + 2:
            raise ValueError(
                f"A trajectory needs at least history_size + 2 = {history_size + 2} states; "
                f"got {steps}.")
        velocities = numpy.diff(trajectory, axis=0) / delta_time  # v_t at index t-1
        for t in range(history_size, steps - 1):
            window = [velocities[t - history_size + k] for k in range(history_size)]
            features_list.append(numpy.concatenate(window, axis=1))
            targets_list.append(
                (trajectory[t + 1] - 2.0 * trajectory[t] + trajectory[t - 1]) / delta_time ** 2)
            positions_list.append(trajectory[t])

    features = numpy.stack(features_list)
    targets = numpy.stack(targets_list)

    feature_mean = features.mean(axis=(0, 1))
    feature_std = features.std(axis=(0, 1))
    feature_std[feature_std == 0.0] = 1.0
    target_mean = targets.mean(axis=(0, 1))
    target_std = targets.std(axis=(0, 1))
    target_std[target_std == 0.0] = 1.0
    if normalize:
        features = (features - feature_mean) / feature_std
        targets = (targets - target_mean) / target_std

    dataset = torch.utils.data.TensorDataset(
        torch.from_numpy(features), torch.from_numpy(targets))
    dataset.positions = positions_list
    dataset.feature_mean = feature_mean
    dataset.feature_std = feature_std
    dataset.target_mean = target_mean
    dataset.target_std = target_std
    return dataset


def MakeNormalizationCardEntries(dataset) -> dict:
    """The model-card entries a model trained on a normalize=True particle
    dataset needs: {"input_normalization": ..., "output_normalization": ...}
    from the dataset's feature_mean/std and target_mean/std.

    Both halves of the standardization then travel with the checkpoint -
    ParticleInferenceProcess standardizes the velocity history it gathers
    with the first and de-normalizes the predicted acceleration with the
    second. Merge the result into the card dict handed to SaveTrainedModel.
    """
    from KratosMultiphysics.PhysicsNeMoApplication.deployment import model_registry
    for name in ("feature_mean", "feature_std", "target_mean", "target_std"):
        if not hasattr(dataset, name):
            raise ValueError(
                f"The dataset has no \"{name}\" attribute; pass the dataset returned by "
                "CreateParticleTrajectoryDataset.")
    return {
        "input_normalization": model_registry.MakeMeanStdNormalization(
            dataset.feature_mean, dataset.feature_std),
        "output_normalization": model_registry.MakeMeanStdNormalization(
            dataset.target_mean, dataset.target_std),
    }


def _TryImportMeshTransforms():
    try:
        from physicsnemo.datapipes import (
            MultiDataset, RandomRotateMesh, RandomScaleMesh, RandomTranslateMesh)
        return MultiDataset, RandomRotateMesh, RandomScaleMesh, RandomTranslateMesh
    except ImportError as e:
        raise ImportError(
            "PhysicsNeMoApplication.torch_dataset's augmentation factories require "
            "physicsnemo, which could not be imported. Install it with e.g. "
            "'pip install nvidia-physicsnemo'.") from e


def _MakeCastMeshTransform(dtype=None):
    """A MeshTransform casting a mesh's floating-point data to one dtype.

    Needed because the upstream random augmentations build their rotation
    and scaling factors in torch's GLOBAL default dtype (float32), while
    Kratos exports meshes in float64 - mixing them raises "expected m1 and
    m2 to have the same dtype". Prepending this cast keeps the whole
    pipeline in the dtype models train in.
    """
    torch = _TryImportTorch()
    from physicsnemo.datapipes.transforms.mesh.base import MeshTransform

    target = dtype if dtype is not None else torch.get_default_dtype()
    if isinstance(target, str):
        target = getattr(torch, target)

    class CastMeshDtype(MeshTransform):
        def __call__(self, mesh):
            def cast(tensor):
                return tensor.to(target) if tensor.is_floating_point() else tensor

            return mesh.replace(
                points=cast(mesh.points),
                point_data=mesh.point_data.apply(cast) if mesh.point_data is not None else None,
                cell_data=mesh.cell_data.apply(cast) if mesh.cell_data is not None else None)

    return CastMeshDtype()


def _AsRange(option, default_low, default_high, name):
    """True -> defaults; {"low": .., "high": ..} -> the given range."""
    if option is True:
        return default_low, default_high
    if isinstance(option, dict):
        return option.get("low", default_low), option.get("high", default_high)
    raise ValueError(f"\"{name}\" must be True or a dict with \"low\"/\"high\", got {option!r}.")


def MakeMeshAugmentations(rotation=None, scale=None, translation=None,
                          vector_fields=(), tensor_fields=(), dtype=None):
    """Builds coherent random-augmentation transforms for mesh datasets.

    The transforms rotate/scale the listed VECTOR (N, 3) and rank-2 TENSOR
    (N, 3, 3) point_data fields together with the coordinates - upstream's
    `transform_point_data` defaults to False, which silently leaves vector
    fields in the old frame, and its bare `True` form raises on any
    non-spatial feature field; this factory always passes the per-field dict
    form, so unlisted fields pass through untouched. Translation never
    touches field values (correctly so).

    Args:
        rotation: None, True (uniform SO(3)), or a dict:
            {"mode": "uniform"} - uniform random 3D rotations (3D meshes only);
            {"mode": "axis_aligned", "axes": ["x","y","z"], "low": -pi, "high": pi}
            - one axis picked uniformly, angle drawn from Uniform(low, high).
        scale: None, True (Uniform(0.9, 1.1)) or {"low": .., "high": ..} -
            isotropic random scaling.
        translation: None, True (Uniform(-0.1, 0.1) per axis) or
            {"low": .., "high": ..} with floats or per-axis lists.
        vector_fields: point_data field names transforming as vectors
            (v -> R v), e.g. the "VELOCITY" the mesh export wrote.
        tensor_fields: field names transforming as rank-2 tensors
            (T -> R T R^T), e.g. stress.
        dtype: Floating dtype the meshes are cast to first (default:
            torch's global default, i.e. float32). The upstream transforms
            build their rotation/scale factors in that dtype, so
            float64 Kratos meshes must be cast to match.

    Returns:
        List of physicsnemo MeshTransform instances, in application order
        (pass to CreateMeshDataset(transforms=...) - physicsnemo's Compose
        does NOT accept mesh transforms; the dataset's transform list is
        the composition mechanism).
    """
    torch = _TryImportTorch()
    _, RandomRotateMesh, RandomScaleMesh, RandomTranslateMesh = _TryImportMeshTransforms()

    transform_point_data = {name: True for name in (*vector_fields, *tensor_fields)}
    transforms = [_MakeCastMeshTransform(dtype)]

    if rotation is not None:
        options = {} if rotation is True else dict(rotation)
        mode = options.get("mode", "uniform")
        if mode == "uniform":
            transforms.append(RandomRotateMesh(
                mode="uniform", transform_point_data=transform_point_data))
        elif mode == "axis_aligned":
            low, high = _AsRange(options, -numpy.pi, numpy.pi, "rotation")
            transforms.append(RandomRotateMesh(
                mode="axis_aligned", axes=options.get("axes", ["x", "y", "z"]),
                distribution=torch.distributions.Uniform(float(low), float(high)),
                transform_point_data=transform_point_data))
        else:
            raise ValueError(
                f"Unknown rotation mode \"{mode}\". Use \"uniform\" or \"axis_aligned\".")

    if scale is not None:
        low, high = _AsRange(scale, 0.9, 1.1, "scale")
        transforms.append(RandomScaleMesh(
            distribution=torch.distributions.Uniform(float(low), float(high)),
            transform_point_data=transform_point_data))

    if translation is not None:
        low, high = _AsRange(translation, -0.1, 0.1, "translation")
        distribution = torch.distributions.Uniform(
            torch.as_tensor(low, dtype=torch.float64),
            torch.as_tensor(high, dtype=torch.float64))
        transforms.append(RandomTranslateMesh(distribution=distribution))

    return transforms if len(transforms) > 1 else []


def CreateAugmentedMeshDataset(directory, rotation=None, scale=None, translation=None,
                               vector_fields=(), tensor_fields=(), seed=-1,
                               extra_transforms=(), dtype=None, device=None, num_workers=1):
    """A MeshDataset over a .pmsh series with coherent random augmentations.

    Randomness is redrawn on every __getitem__ (each epoch sees new
    augmentations). With seed >= 0 the dataset's generator is seeded, making
    the draw sequence reproducible; call dataset.set_epoch(e) per epoch to
    reseed deterministically per epoch.

    Args:
        directory: The MeshExportProcess output directory.
        rotation/scale/translation/vector_fields/tensor_fields: see
            MakeMeshAugmentations.
        seed: -1 (default) leaves the RNG alone; >= 0 seeds it.
        extra_transforms: Additional MeshTransforms appended after the
            augmentations (e.g. NormalizeMeshFields, MeshToTensorDict).
        device: Optional device the meshes are moved to.
        num_workers: Reader worker count.
    """
    torch = _TryImportTorch()
    transforms = MakeMeshAugmentations(rotation, scale, translation,
                                       vector_fields, tensor_fields, dtype)
    transforms.extend(extra_transforms)
    dataset = CreateMeshDataset(directory, transforms=transforms,
                                device=device, num_workers=num_workers)
    if seed >= 0:
        dataset.set_generator(torch.Generator().manual_seed(int(seed)))
    return dataset


def CreateMultiMeshDataset(sources, output_strict=True, seed=-1, **augmentation_kwargs):
    """Mixes several .pmsh series (or ready datasets) into one dataset.

    Wraps physicsnemo's MultiDataset: one concatenated index space, each
    item identical to the owning sub-dataset's, with metadata extended by
    "dataset_index". Notes from upstream: with output_strict=True (default)
    sample 0 of EVERY sub-dataset is loaded eagerly at construction to
    check the outputs are stackable, and shuffling does NOT balance
    sub-datasets - larger series dominate proportionally.

    Args:
        sources: Iterable of directories (each becomes a
            CreateAugmentedMeshDataset with the shared augmentation_kwargs)
            and/or ready dataset objects (used as-is).
        output_strict: physicsnemo MultiDataset strictness (see above).
        seed: With seed >= 0, sub-dataset i built from a directory is seeded
            with seed + i (independent but reproducible streams).
        **augmentation_kwargs: Forwarded to CreateAugmentedMeshDataset for
            directory sources (rotation, scale, translation, vector_fields,
            tensor_fields, extra_transforms, device, num_workers).
    """
    MultiDataset, _, _, _ = _TryImportMeshTransforms()

    datasets = []
    for index, source in enumerate(sources):
        if isinstance(source, (str, Path)):
            sub_seed = seed + index if seed >= 0 else -1
            datasets.append(CreateAugmentedMeshDataset(
                source, seed=sub_seed, **augmentation_kwargs))
        else:
            datasets.append(source)
    if not datasets:
        raise ValueError("CreateMultiMeshDataset needs at least one source.")
    return MultiDataset(*datasets, output_strict=output_strict)
