"""Differentiable shape deformation on Kratos meshes.

The design-parameterization layer: a small set of control parameters maps
differentiably to node coordinates, so gradients of a surrogate objective
flow back to the shape itself (physicsnemo >= 2.2's
physicsnemo.nn.functional deformers, Warp-accelerated on CUDA).

    control displacements --DeformPoints--> deformed coordinates
                                                |
                                     WriteNodeCoordinates
                                                |
                                        a moved Kratos mesh

Pair it with sensitivity_utils.ComputeShapeSensitivities to get dJ/d(control)
by the chain rule through a surrogate, and with RegularizationEnergy to keep
the deformed mesh valid (element inversion is the failure mode that ruins a
shape optimization run).

Two contracts worth stating up front, both easy to get wrong:

- **Control values are DISPLACEMENTS, not destination coordinates.** A zero
  control array is the identity deformation.
- **Coordinate write-back is the only mutating operation here.** Everything
  else is pure: it returns tensors and leaves the model part alone.

Mesh moving: writing coordinates directly is the right thing for a shape
parameterization of the *boundary and volume together* (FFD/RBF move every
node they cover). When only a surface should move and the interior should
follow smoothly, drive MeshMovingApplication with the resulting boundary
displacement instead of writing interior nodes here.

torch/physicsnemo are optional runtime dependencies - imported lazily.
"""

import numpy

import KratosMultiphysics as Kratos

_METHODS = ("ffd", "rbf", "displace", "morph")


def _TryImportTorch():
    try:
        import torch
        return torch
    except ImportError as e:
        raise ImportError(
            "PhysicsNeMoApplication.mesh_bridge.deformation requires torch, which could "
            "not be imported. Install it with e.g. 'pip install torch'.") from e


def _TryImportDeformers():
    try:
        from physicsnemo.nn.functional import (
            displace_points, free_form_deform_points, morph_points,
            radial_basis_function_deform_points)
        return {"ffd": free_form_deform_points,
                "rbf": radial_basis_function_deform_points,
                "displace": displace_points,
                "morph": morph_points}
    except ImportError as e:
        raise ImportError(
            "PhysicsNeMoApplication.mesh_bridge.deformation requires physicsnemo >= 2.2 "
            "(its deformers landed in 2.2), which could not be imported. Install it with "
            "e.g. 'pip install -U nvidia-physicsnemo'.") from e


def _TryImportEnergies():
    try:
        from physicsnemo.mesh import deformation as energies
        return energies
    except ImportError as e:
        raise ImportError(
            "PhysicsNeMoApplication.mesh_bridge.deformation requires physicsnemo >= 2.2 "
            "(its deformation energies landed in 2.2), which could not be imported. "
            "Install it with e.g. 'pip install -U nvidia-physicsnemo'.") from e


def DeformPoints(points, control_displacements, method: str = "ffd", **options):
    """Deforms point coordinates from a differentiable control parameterization.

    Args:
        points: (N, D) coordinates to deform (tensor or array). Pass a tensor
            with requires_grad to differentiate w.r.t. the shape itself.
        control_displacements: The control parameters, per method - and always
            DISPLACEMENTS, never destinations:
            "ffd"      -> (*lattice_resolution, D) lattice displacements;
            "rbf"      -> (C, D) displacements at `control_points`;
            "morph"    -> (C, D) displacements at `control_points`;
            "displace" -> (N, D) one displacement per point.
        method: One of _METHODS.
        **options: Forwarded to the upstream deformer. Required per method:
            "ffd" needs origin/extent (the lattice's bounding box; defaults to
            the points' own bounding box here), "rbf"/"morph" need
            control_points, and "morph" also needs radius.

    Returns:
        (N, D) deformed coordinates, differentiable w.r.t. both `points` and
        `control_displacements`.
    """
    torch = _TryImportTorch()
    deformers = _TryImportDeformers()
    if method not in _METHODS:
        raise ValueError(f"Unknown deformation method \"{method}\". Use one of {_METHODS}.")

    points = torch.as_tensor(points)
    control_displacements = torch.as_tensor(
        control_displacements, dtype=points.dtype, device=points.device)

    if method == "ffd" and "origin" not in options:
        # a lattice must span the points it deforms; the points' own bounding
        # box is the only sane default and makes the identity case exact
        options.setdefault("origin", points.min(dim=0).values.tolist())
        options.setdefault("extent",
                           (points.max(dim=0).values - points.min(dim=0).values).tolist())
    for key in ("control_points",):
        if key in options and options[key] is not None:
            options[key] = torch.as_tensor(options[key], dtype=points.dtype,
                                           device=points.device)

    if method in ("rbf", "morph"):
        if options.get("control_points") is None:
            raise ValueError(f"The \"{method}\" method needs \"control_points\".")
        control_points = options.pop("control_points")
        return deformers[method](points, control_points, control_displacements, **options)
    return deformers[method](points, control_displacements, **options)


# The 2.2 deformers that take a MESH rather than a point cloud: they need
# connectivity (a target surface to project onto, a stiffness matrix to
# filter through), so they cannot join DeformPoints' dispatch.
_MESH_METHODS = ("shrinkwrap", "sobolev")


def _TryImportMeshDeformers():
    try:
        from physicsnemo.mesh.transformations.deform.shrinkwrap import shrinkwrap
        from physicsnemo.mesh.transformations.deform.sobolev import sobolev_deform
        return {"shrinkwrap": shrinkwrap, "sobolev": sobolev_deform}
    except ImportError as e:
        raise ImportError(
            "The mesh deformers require physicsnemo >= 2.2, which could not be "
            "imported. Install it with e.g. 'pip install -U nvidia-physicsnemo'. "
            "(They live in physicsnemo.mesh.transformations.deform, NOT in "
            "physicsnemo.nn.functional with the point deformers.)") from e


def DeformMesh(mesh, method: str, **options):
    """A mesh-level deformation: shrinkwrap or a Sobolev displacement filter.

    - "shrinkwrap" projects every vertex onto the nearest point of a
      ``target`` surface, optionally with an ``offset`` along the target's
      normal. The classic use is fitting a coarse control surface onto a
      detailed one.
    - "sobolev" solves (M + l^2 K) u = M d, filtering a raw per-vertex
      ``displacement`` through the mesh's own stiffness so the deformation
      spreads smoothly instead of denting single nodes. ``fixed_points``
      pins vertices. This is what MeshMovingApplication does for the
      shape-optimization loop, in torch and differentiable.

    Args:
        mesh: The physicsnemo Mesh to deform (it is not modified).
        method: "shrinkwrap" or "sobolev".
        options: The method's own arguments - ``target`` (and optionally
          ``offset``, ``max_distance``, ``point_weights``) for shrinkwrap;
          ``displacement`` and ``length_scale`` (and optionally
          ``fixed_points``, ``max_iterations``, ``tolerance``) for sobolev.

    Returns:
        A new physicsnemo Mesh with the deformed points.
    """
    if method not in _MESH_METHODS:
        raise ValueError(
            f"Unsupported mesh deformation \"{method}\". Supported: "
            f"{', '.join(_MESH_METHODS)}. The point-cloud deformations "
            f"({', '.join(_METHODS)}) go through DeformPoints.")
    deformers = _TryImportMeshDeformers()
    if method == "shrinkwrap":
        if "target" not in options:
            raise ValueError(
                "\"shrinkwrap\" needs a \"target\" surface to project onto.")
        target = options.pop("target")
        return deformers["shrinkwrap"](mesh, target, **options)

    for required in ("displacement", "length_scale"):
        if required not in options:
            raise ValueError(
                f"\"sobolev\" needs \"{required}\"; the filter is "
                "(M + length_scale^2 K) u = M displacement.")
    displacement = options.pop("displacement")
    return deformers["sobolev"](mesh, displacement, **options)


def RegularizationEnergy(reference_mesh, points, energy: str = "strain", **options):
    """Mesh-quality energy of a deformed configuration, as an objective term.

    These are penalties to ADD to a design objective, not constraints:
    "strain" resists distortion away from the reference shape, "inversion"
    blows up as elements approach zero/negative volume (the term that keeps a
    shape optimization from tearing the mesh), "measure" targets a volume
    ratio, and "bending" penalizes surface curvature change.

    Degenerate reference cells deliberately produce NaN upstream rather than
    being silently regularized - treat a NaN here as "your reference mesh is
    broken", not as a numerical hiccup.

    Args:
        reference_mesh: physicsnemo.mesh.Mesh in the REFERENCE configuration.
        points: (N, D) deformed coordinates (differentiable).
        energy: "strain", "inversion", "measure", "total_measure", "bending"
            or "volume".
        **options: Forwarded upstream (lame_lambda/shear_modulus for strain,
            minimum_jacobian for inversion, target_ratio for the measure
            energies, reduction for all of them).

    Returns:
        A scalar tensor by default (reduction="sum"), differentiable w.r.t.
        `points`.
    """
    torch = _TryImportTorch()
    energies = _TryImportEnergies()
    table = {
        "strain": energies.simplex_strain_energy,
        "inversion": energies.simplex_inversion_energy,
        "measure": energies.simplex_measure_energy,
        "total_measure": energies.total_measure_energy,
        "bending": energies.surface_bending_energy,
        "volume": energies.closed_surface_volume_energy,
    }
    if energy not in table:
        raise ValueError(
            f"Unknown energy \"{energy}\". Use one of {tuple(table)}.")
    return table[energy](reference_mesh, torch.as_tensor(points), **options)


def WriteNodeCoordinates(model_part: Kratos.ModelPart, coordinates,
                         update_displacement: bool = False) -> None:
    """Writes deformed coordinates onto a model part's nodes.

    The mutating counterpart of graph_bridge.NodePositions: rows are in
    model_part.Nodes iteration order, matching every gather in this
    application.

    Args:
        model_part: The model part to move.
        coordinates: (N, D) deformed coordinates.
        update_displacement: Also store (X - X0) in DISPLACEMENT, which is
            what MeshMovingApplication and the structural solvers read.
    """
    coordinates = numpy.ascontiguousarray(
        numpy.asarray(coordinates, dtype=numpy.float64))
    if coordinates.shape[0] != model_part.NumberOfNodes():
        raise ValueError(
            f"coordinates has {coordinates.shape[0]} rows but the model part has "
            f"{model_part.NumberOfNodes()} nodes.")

    position_ta = Kratos.TensorAdaptors.NodePositionTensorAdaptor(
        model_part.Nodes, Kratos.Configuration.Current)
    position_ta.CollectData()
    position_ta.data[:] = coordinates
    position_ta.StoreData()

    if update_displacement:
        for node, position in zip(model_part.Nodes, coordinates):
            node.SetSolutionStepValue(
                Kratos.DISPLACEMENT,
                [float(position[0] - node.X0), float(position[1] - node.Y0),
                 float(position[2] - node.Z0)])


def DeformModelPart(model_part: Kratos.ModelPart, control_displacements,
                    method: str = "ffd", update_displacement: bool = False,
                    **options) -> numpy.ndarray:
    """Deforms a model part in place and returns the new coordinates.

    Convenience wrapper: read coordinates, deform, write back. Use
    DeformPoints directly when the deformation must stay in the autograd
    graph (this one detaches, because Kratos coordinates are plain data).
    """
    from KratosMultiphysics.PhysicsNeMoApplication.bridges import graph_bridge
    torch = _TryImportTorch()
    points = torch.as_tensor(graph_bridge.NodePositions(model_part))
    deformed = DeformPoints(points, control_displacements, method, **options)
    coordinates = deformed.detach().cpu().numpy()
    WriteNodeCoordinates(model_part, coordinates, update_displacement)
    return coordinates
