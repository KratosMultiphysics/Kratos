"""Point-cloud token budgets: run a model on a subset of the nodes.

A point-cloud transformer's cost grows with the number of tokens, and a
Kratos mesh routinely has more nodes than such a model was ever trained to
attend over. The DoMINO datapipe subsamples; the point-cloud deployment
path used to feed every node.

This module supplies the two halves of doing that honestly:

- SelectPointSubset picks the subset - a bounding box filter, then
  farthest-point sampling (physicsnemo.nn.functional, which spreads the
  points over the geometry) or a seeded uniform draw.
- ExpandToAllPoints puts the prediction back on EVERY node, giving an
  unselected node its nearest selected node's value, so the field a solver
  reads afterwards is complete whatever budget was used.

Farthest-point sampling is the one to prefer: a uniform draw of a mesh
refined in one corner spends its budget there, while FPS covers the whole
geometry. The difference is measurable as the minimum pairwise distance of
the chosen points, which is what the tests assert.

torch, physicsnemo and scipy are imported lazily.
"""

import numpy

import KratosMultiphysics as Kratos

_SUBSAMPLING_METHODS = ("none", "farthest_point", "uniform")


def _TryImportTorch():
    try:
        import torch
        return torch
    except ImportError as e:
        raise ImportError(
            "PhysicsNeMoApplication.point_subsampling requires torch, which could not "
            "be imported. Install it with e.g. 'pip install torch'.") from e


def _TryImportFarthestPointSampling():
    try:
        from physicsnemo.nn.functional import farthest_point_sampling
        return farthest_point_sampling
    except ImportError as e:
        raise ImportError(
            "The \"farthest_point\" subsampling method requires physicsnemo, which "
            "could not be imported. Install it with e.g. "
            "'pip install -U nvidia-physicsnemo'.") from e


def _TryImportKdTree():
    try:
        from scipy.spatial import cKDTree
        return cKDTree
    except ImportError as e:
        raise ImportError(
            "Expanding a subsampled prediction back onto every node requires scipy, "
            "which could not be imported. Install it with e.g. "
            "'pip install scipy'.") from e


def SubsamplingDefaults() -> Kratos.Parameters:
    """The defaults of a "subsampling" settings block."""
    return Kratos.Parameters("""{
        "method"           : "none",
        "num_points"       : 0,
        "bounding_box_min" : [],
        "bounding_box_max" : [],
        "seed"             : -1
    }""")


def SelectPointSubset(coordinates, settings: Kratos.Parameters):
    """Indices of the points a model should actually see.

    The bounding box is applied first (it is a filter, not a budget), then
    the method reduces whatever survived to "num_points".

    Args:
        coordinates: (N, 3) array-like node coordinates.
        settings: Kratos Parameters; defaults from SubsamplingDefaults() -
            method ("none" | "farthest_point" | "uniform"), num_points (0 =
            no budget), bounding_box_min/max ([] = no filter), seed (-1).

    Returns:
        A sorted (n_selected,) int64 numpy array of indices into
        coordinates, or None when nothing is selected away (so a caller can
        skip the whole expansion path).
    """
    settings.ValidateAndAssignDefaults(SubsamplingDefaults())
    method = settings["method"].GetString()
    if method not in _SUBSAMPLING_METHODS:
        raise ValueError(
            f"Unsupported subsampling method \"{method}\". Supported: "
            f"{', '.join(_SUBSAMPLING_METHODS)}.")

    coordinates = numpy.asarray(coordinates, dtype=numpy.float64)
    n_points = coordinates.shape[0]
    indices = numpy.arange(n_points, dtype=numpy.int64)

    low = list(settings["bounding_box_min"].GetVector())
    high = list(settings["bounding_box_max"].GetVector())
    if low or high:
        if len(low) != 3 or len(high) != 3:
            raise ValueError(
                "\"bounding_box_min\" and \"bounding_box_max\" must both have three "
                f"entries; got {len(low)} and {len(high)}.")
        inside = numpy.all(
            (coordinates >= numpy.array(low)) & (coordinates <= numpy.array(high)), axis=1)
        indices = indices[inside]
        if indices.size == 0:
            raise ValueError(
                "The subsampling bounding box contains no nodes; it is given in the "
                "model part's own coordinates, not normalized ones.")

    budget = settings["num_points"].GetInt()
    if method == "none" or budget <= 0 or budget >= indices.size:
        # nothing was filtered and nothing is over budget: no subset at all
        return None if indices.size == n_points else indices

    seed = settings["seed"].GetInt()
    if method == "uniform":
        generator = numpy.random.default_rng(None if seed < 0 else seed)
        chosen = generator.choice(indices.size, size=budget, replace=False)
    else:  # farthest_point
        torch = _TryImportTorch()
        farthest_point_sampling = _TryImportFarthestPointSampling()
        if seed >= 0:
            torch.manual_seed(seed)
        chosen = farthest_point_sampling(
            torch.from_numpy(coordinates[indices]), budget,
            random_start=seed >= 0).cpu().numpy()
    return numpy.sort(indices[chosen])


def ExpandToAllPoints(values, coordinates, indices):
    """Puts a subset's prediction back onto every point.

    Each unselected point takes its nearest selected point's row, so the
    written field is complete - a solver reading it cannot tell which nodes
    the model actually saw, which is the point of a budget.

    Args:
        values: (n_selected, C) array-like or torch tensor of predictions.
        coordinates: (N, 3) coordinates of every point.
        indices: the (n_selected,) indices SelectPointSubset returned.

    Returns:
        The same type as values, with N rows.
    """
    if values is None:
        return None
    coordinates = numpy.asarray(coordinates, dtype=numpy.float64)
    indices = numpy.asarray(indices, dtype=numpy.int64)
    cKDTree = _TryImportKdTree()
    _, nearest = cKDTree(coordinates[indices]).query(coordinates, k=1)
    nearest = numpy.asarray(nearest, dtype=numpy.int64)

    if hasattr(values, "detach"):  # a torch tensor: stay one
        torch = _TryImportTorch()
        return values[torch.from_numpy(nearest)]
    return numpy.asarray(values)[nearest]
