"""Discrete calculus on tessellated Kratos meshes via physicsnemo.mesh.calculus.

Solver-free differential operators (gradient, divergence, curl, Laplacian)
and integrals evaluated directly on the simplex mesh the mesh bridge builds
- for physics-residual features, conserved-quantity monitoring and feature
engineering, without a builder-and-solver assembly. For the physics' own
assembled residual use solver_residuals / differentiable_residual instead.

Backend validity (probed against physicsnemo 2.2 and enforced here):

- LSQ operators are correct on both surface (codim-1) and volume (codim-0)
  meshes; DEC gradient/divergence are silently WRONG on volume (e.g.
  tetrahedral) meshes, so this module refuses them there. The DEC Laplacian
  IS correct on volume meshes - but only at interior points.
- Multi-channel gradients come back derivative-first (N, D, C) from every
  backend; this module transposes them once to its own channel-major
  (N, C, D) contract. physicsnemo < 2.2 returned the LSQ gradient the other
  way round, which is why _TryImportPhysicsNemoCalculus enforces a minimum
  version instead of guessing.
- None of the upstream operators handle boundaries: DEC results are garbage
  at boundary points and LSQ degrades there on coarse meshes. Use
  InteriorPointMask to mask them.
- Surface meshes need the intrinsic (tangent-plane) LSQ gradient; the
  extrinsic default carries an ill-conditioned normal component. Volume
  meshes are insensitive (both coincide). ComputeGradient picks the right
  one from the mesh's codimension unless overridden. (2.2 also fixed two
  upstream intrinsic bugs: multi-channel fields no longer raise, and
  codimension >= 2 meshes no longer return all zeros.)
- All operators are autograd-differentiable (w.r.t. the field values and
  even mesh.points); float32/float64 only (fp16 fails in the LSQ solve).

torch/physicsnemo are optional runtime dependencies - imported lazily.
"""

import numpy

import KratosMultiphysics as Kratos

from KratosMultiphysics.PhysicsNeMoApplication.bridges.mesh_bridge import domain_mesh_builder


def _TryImportTorch():
    try:
        import torch
        return torch
    except ImportError as e:
        raise ImportError(
            "PhysicsNeMoApplication.calculus_bridge requires torch, which could not be "
            "imported. Install it with e.g. 'pip install torch'.") from e


_MINIMUM_PHYSICSNEMO = (2, 2)


def _TryImportPhysicsNemoCalculus():
    try:
        import physicsnemo
        from physicsnemo.mesh import calculus
    except ImportError as e:
        raise ImportError(
            "PhysicsNeMoApplication.calculus_bridge requires physicsnemo, which could "
            "not be imported. Install it with e.g. 'pip install nvidia-physicsnemo'.") from e

    # Version guard, not pedantry: physicsnemo < 2.2 returned the LSQ gradient
    # of a multi-channel field CHANNEL-first while 2.2 returns it
    # derivative-first. Normalizing for the wrong one produces transposed
    # gradients with no error anywhere - so refuse rather than guess.
    version = getattr(physicsnemo, "__version__", "0.0")
    try:
        parsed = tuple(int(part) for part in version.split(".")[:2])
    except ValueError:
        parsed = _MINIMUM_PHYSICSNEMO  # unparsable (dev build): assume current
    if parsed < _MINIMUM_PHYSICSNEMO:
        raise ImportError(
            f"PhysicsNeMoApplication.calculus_bridge requires physicsnemo >= "
            f"{'.'.join(str(p) for p in _MINIMUM_PHYSICSNEMO)}, found {version}. Older "
            "releases return multi-channel gradients in the opposite axis order, which "
            "this bridge would silently mis-normalize. Upgrade with e.g. "
            "'pip install -U nvidia-physicsnemo'.")
    return calculus


def _IsVolumeMesh(mesh) -> bool:
    """codim 0: the manifold fills its ambient space (tets in 3D, tris in 2D)."""
    return int(mesh.points.shape[1]) == int(mesh.cells.shape[1]) - 1


def _CheckDecSupported(mesh, operation: str) -> None:
    if _IsVolumeMesh(mesh):
        raise ValueError(
            f"The DEC {operation} is not valid on volume (codimension-0) meshes - "
            "physicsnemo silently returns wrong values there (still true in 2.2, which "
            "adds no validation of its own). Use method=\"lsq\".")


def ComputeGradient(mesh, values, method: str = "lsq", intrinsic=None):
    """Per-point gradient of a point field.

    Args:
        mesh: physicsnemo.mesh.Mesh (from e.g. domain_mesh_builder.BuildMesh).
        values: (N,) or (N, C) tensor/array of per-point values.
        method: "lsq" (default; valid everywhere) or "dec" (surface/curve
            meshes only - refused on volume meshes where it is silently wrong).
        intrinsic: None (default) picks per codimension: tangent-plane
            (intrinsic) LSQ on surfaces, extrinsic on volume meshes. Pass
            True/False to override; ignored for "dec".

    Returns:
        (N, D) tensor for scalar input, (N, C, D) for multi-channel input -
        i.e. gradient[i, c, d] is d(values[:, c])/dx_d. This CHANNEL-MAJOR
        layout is the bridge's stable contract; upstream returns
        derivative-first (N, D, C) from every backend as of physicsnemo 2.2,
        and is transposed here exactly once so callers never see the
        difference.
    """
    torch = _TryImportTorch()
    calculus = _TryImportPhysicsNemoCalculus()

    values = torch.as_tensor(values)
    if values.dim() not in (1, 2):
        raise ValueError(f"values must be (N,) or (N, C), got shape {tuple(values.shape)}.")

    if method == "dec":
        _CheckDecSupported(mesh, "gradient")
        gradient = calculus.compute_gradient_points_dec(mesh, values)
    elif method == "lsq":
        if intrinsic is None:
            intrinsic = not _IsVolumeMesh(mesh)
        gradient = calculus.compute_gradient_points_lsq(mesh, values, intrinsic=intrinsic)
    else:
        raise ValueError(f"Unknown gradient method \"{method}\". Use \"lsq\" or \"dec\".")

    if values.dim() == 2:
        gradient = gradient.movedim(1, 2)  # (N, D, C) -> (N, C, D)
    return gradient


def ComputeDivergence(mesh, vector_values, method: str = "lsq"):
    """Per-point divergence of a (N, D) point vector field -> (N,) tensor."""
    torch = _TryImportTorch()
    calculus = _TryImportPhysicsNemoCalculus()

    vector_values = torch.as_tensor(vector_values)
    if method == "dec":
        _CheckDecSupported(mesh, "divergence")
        return calculus.compute_divergence_points_dec(mesh, vector_values)
    if method != "lsq":
        raise ValueError(f"Unknown divergence method \"{method}\". Use \"lsq\" or \"dec\".")
    return calculus.compute_divergence_points_lsq(mesh, vector_values)


def ComputeCurl(mesh, vector_values):
    """Per-point curl of a (N, 3) point vector field -> (N, 3) tensor (LSQ)."""
    torch = _TryImportTorch()
    calculus = _TryImportPhysicsNemoCalculus()
    return calculus.compute_curl_points_lsq(mesh, torch.as_tensor(vector_values))


def ComputeLaplacian(mesh, values):
    """Per-point Laplacian of a point field (DEC cotangent formula).

    Valid on surface AND volume meshes, but only at INTERIOR points - the
    formula has no boundary treatment, so mask the result with
    InteriorPointMask before using it.
    """
    torch = _TryImportTorch()
    calculus = _TryImportPhysicsNemoCalculus()
    return calculus.compute_laplacian_points_dec(mesh, torch.as_tensor(values))


def IntegrateField(mesh, values, data_source: str = "points"):
    """Integral of a field over the mesh (float result for scalar fields).

    Args:
        mesh: physicsnemo.mesh.Mesh.
        values: per-point (data_source="points") or per-cell ("cells")
            values, or the name of a point_data/cell_data field.
        data_source: "points" (default) or "cells".
    """
    torch = _TryImportTorch()
    calculus = _TryImportPhysicsNemoCalculus()
    if not isinstance(values, (str, tuple)):
        values = torch.as_tensor(values)
    return calculus.integrate(mesh, values, data_source=data_source)


def InteriorPointMask(mesh):
    """Boolean (N,) tensor: True for interior points, False on the boundary.

    A facet (sub-simplex of one lower dimension) lying in exactly one cell is
    a boundary facet; every point of a boundary facet is a boundary point.
    Pure torch on the connectivity - no optional dependencies beyond torch.
    """
    torch = _TryImportTorch()

    cells = torch.as_tensor(mesh.cells)
    n_points = int(mesh.points.shape[0])
    n_vertices = cells.shape[1]

    facets = []
    for drop in range(n_vertices):
        keep = [v for v in range(n_vertices) if v != drop]
        facets.append(cells[:, keep])
    facets = torch.cat(facets, dim=0).sort(dim=1).values

    unique_facets, counts = torch.unique(facets, dim=0, return_counts=True)
    boundary_points = torch.unique(unique_facets[counts == 1].reshape(-1))

    mask = torch.ones(n_points, dtype=torch.bool, device=cells.device)
    mask[boundary_points] = False
    return mask


_OPERATIONS = ("gradient", "divergence", "curl", "laplacian")


def ComputeNodalDerivatives(model_part: Kratos.ModelPart, settings: Kratos.Parameters):
    """Computes derivative fields of nodal variables and writes them back.

    The feature-engineering entry point: gathers the requested nodal fields
    onto the tessellated mesh, evaluates the requested discrete operators
    and scatters the results back to (non-historical by default) nodal
    variables - physics-consistent derivative features without any solver
    assembly.

    Args:
        model_part: The model part (its Elements are tessellated).
        settings: Kratos Parameters:
            {
                "operations": [
                    { "field"           : "VELOCITY",
                      "field_location"  : "node_historical",
                      "operation"       : "gradient" | "divergence" | "curl" | "laplacian",
                      "output_variable" : "...",
                      "output_location" : "node_non_historical" }
                ],
                "method"        : "lsq",
                "zero_boundary" : false
            }
        With "zero_boundary": true, boundary points (InteriorPointMask) are
        zeroed in every output - recommended for "laplacian" and any DEC op.

    Returns:
        dict {output_variable_name: numpy array written}.
    """
    defaults = Kratos.Parameters("""{
        "operations"    : [],
        "method"        : "lsq",
        "zero_boundary" : false
    }""")
    settings.ValidateAndAssignDefaults(defaults)
    method = settings["method"].GetString()
    zero_boundary = settings["zero_boundary"].GetBool()

    operation_defaults = Kratos.Parameters("""{
        "field"           : "",
        "field_location"  : "node_historical",
        "operation"       : "gradient",
        "output_variable" : "",
        "output_location" : "node_non_historical"
    }""")

    operations = []
    field_specs = []
    for entry in settings["operations"].values():
        entry.ValidateAndAssignDefaults(operation_defaults)
        operation = entry["operation"].GetString()
        if operation not in _OPERATIONS:
            raise ValueError(
                f"Unknown operation \"{operation}\". Use one of {_OPERATIONS}.")
        if not entry["field"].GetString() or not entry["output_variable"].GetString():
            raise ValueError("Each operation needs \"field\" and \"output_variable\".")
        variable = Kratos.KratosGlobals.GetVariable(entry["field"].GetString())
        location = entry["field_location"].GetString()
        if (variable, location) not in field_specs:
            field_specs.append((variable, location))
        operations.append((variable, location, operation,
                           entry["output_variable"].GetString(),
                           entry["output_location"].GetString()))

    mesh, provenance = domain_mesh_builder.BuildMesh(model_part, field_specs)
    mask = InteriorPointMask(mesh) if zero_boundary else None

    written = {}
    for variable, location, operation, output_name, output_location in operations:
        values = mesh.point_data[variable.Name()]
        if operation == "gradient":
            result = ComputeGradient(mesh, values, method=method)
            if result.dim() == 3:  # (N, C, D) vector-field gradient -> flatten
                result = result.reshape(result.shape[0], -1)
        elif operation == "divergence":
            result = ComputeDivergence(mesh, values, method=method)
        elif operation == "curl":
            result = ComputeCurl(mesh, values)
        else:
            result = ComputeLaplacian(mesh, values)
        if mask is not None:
            result = result * mask.reshape(-1, *([1] * (result.dim() - 1)))
        result = result.detach().cpu().numpy()
        domain_mesh_builder.ScatterFieldBack(
            provenance, result, model_part,
            Kratos.KratosGlobals.GetVariable(output_name), output_location)
        written[output_name] = result
    return written


def IntegrateNodalField(model_part: Kratos.ModelPart, variable,
                        data_location: str = "node_historical"):
    """Integral of a nodal field over the tessellated model part.

    Returns a float for scalar variables, a numpy array for vector ones.
    """
    mesh, _ = domain_mesh_builder.BuildMesh(model_part, ((variable, data_location),))
    result = IntegrateField(mesh, mesh.point_data[variable.Name()], data_source="points")
    result = result.detach().cpu().numpy()
    return float(result) if result.ndim == 0 else result
