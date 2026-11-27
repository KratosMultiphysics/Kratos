"""Signed distance fields on Kratos meshes (physicsnemo >= 2.2).

An SDF turns geometry itself into a feature: every point carries its signed
distance to the boundary, which is what lets a surrogate generalize across
shapes instead of memorizing one mesh. This module computes it against a
model part's boundary and, crucially, delivers it the way the rest of the
application already consumes data.

The integration trick: `WriteSignedDistanceField` stores the result into an
ordinary nodal (non-historical) Kratos variable. Every gather in this
application - `grid_bridge.SampleFieldsOnGrid`, `graph_bridge.BuildGraph`,
`inference_process.GatherInputFields` - keys off
`(variable_name, data_location)`, so once the SDF is a variable it flows
into grids, graphs and point clouds through their existing `input_fields`
settings with no signature changes anywhere. `SampleSignedDistanceOnGrid`
covers the one case that cannot work that way: lattice points are not nodes.

Upstream contract, encoded here:

- **Negative inside, positive outside**, zero on the surface.
- It requires a **triangle surface in 3D** and rejects volume meshes, so a
  tetrahedral model part is reduced to its boundary surface first.
- Distances are computed in float32 internally (Warp-backed, CPU and CUDA);
  results are returned as float64 to match this application's convention.
- With `max_dist` set, queries beyond that band return NaN rather than a
  distance - a narrow-band optimization, not an error.

torch/physicsnemo are optional runtime dependencies - imported lazily.
"""

import numpy

import KratosMultiphysics as Kratos

from KratosMultiphysics.PhysicsNeMoApplication.bridges.mesh_bridge import domain_mesh_builder
from KratosMultiphysics.PhysicsNeMoApplication.utilities.tensor_adaptor_dataset_utils import GetTensorAdaptor


def _TryImportTorch():
    try:
        import torch
        return torch
    except ImportError as e:
        raise ImportError(
            "PhysicsNeMoApplication.mesh_bridge.spatial requires torch, which could not "
            "be imported. Install it with e.g. 'pip install torch'.") from e


def _TryImportSignedDistanceField():
    try:
        from physicsnemo.mesh.spatial import signed_distance_field
        return signed_distance_field
    except ImportError as e:
        raise ImportError(
            "PhysicsNeMoApplication.mesh_bridge.spatial requires physicsnemo >= 2.2 "
            "(signed_distance_field landed in 2.2), which could not be imported. Install "
            "it with e.g. 'pip install -U nvidia-physicsnemo'.") from e


def BoundarySurface(mesh):
    """The OUTWARD-oriented triangle surface of a mesh, as the SDF requires.

    A volume (tetrahedral) mesh is rejected by the SDF query, so its boundary
    is extracted here - but extraction alone is not enough: physicsnemo's
    Mesh.get_boundary_mesh() returns triangles with inconsistent winding (the
    closed surface's signed volume comes out 0 instead of the enclosed
    volume), and BOTH sign methods then give meaningless signs, reporting
    interior points as outside.

    So the boundary is rebuilt here with a consistent outward orientation,
    exactly rather than heuristically: a facet belongs to the boundary when it
    appears in a single cell, and its winding is fixed by requiring the normal
    to point AWAY from that cell's opposite vertex. The result satisfies
    signed-volume == enclosed volume, which the tests assert.

    Surface meshes are passed through untouched - their orientation is the
    caller's own.
    """
    torch = _TryImportTorch()
    import physicsnemo.mesh

    if int(mesh.cells.shape[1]) - 1 == 2:
        return mesh

    points = torch.as_tensor(mesh.points)
    cells = torch.as_tensor(mesh.cells)
    n_vertices = int(cells.shape[1])

    facets, opposite = [], []
    for drop in range(n_vertices):
        keep = [v for v in range(n_vertices) if v != drop]
        facets.append(cells[:, keep])
        opposite.append(cells[:, drop])
    facets = torch.cat(facets, dim=0)
    opposite = torch.cat(opposite, dim=0)

    _, inverse, counts = torch.unique(
        facets.sort(dim=1).values, dim=0, return_inverse=True, return_counts=True)
    on_boundary = counts[inverse] == 1
    triangles = facets[on_boundary].clone()
    opposite = opposite[on_boundary]

    corner = points[triangles[:, 0]]
    normal = torch.cross(points[triangles[:, 1]] - corner,
                         points[triangles[:, 2]] - corner, dim=1)
    points_inward = torch.einsum("ij,ij->i", normal, points[opposite] - corner) > 0
    triangles[points_inward] = triangles[points_inward][:, [0, 2, 1]]

    return physicsnemo.mesh.Mesh(points=points, cells=triangles)


def ComputeSignedDistance(mesh, query_points, max_dist=None,
                          use_sign_winding_number: bool = False):
    """Signed distance from arbitrary points to a mesh's boundary surface.

    Args:
        mesh: physicsnemo.mesh.Mesh (volume meshes are reduced to their
            boundary surface automatically).
        query_points: (P, 3) coordinates to evaluate at.
        max_dist: Optional narrow band; queries beyond it return NaN.
        use_sign_winding_number: False (default) uses the angle-weighted
            pseudo-normal sign, which assumes a watertight surface; True uses
            the winding number, which is robust on non-watertight geometry.

    Returns:
        (P,) float64 numpy array - negative inside, positive outside.
    """
    torch = _TryImportTorch()
    signed_distance_field = _TryImportSignedDistanceField()

    surface = BoundarySurface(mesh)
    points = torch.as_tensor(
        numpy.ascontiguousarray(numpy.asarray(query_points, dtype=numpy.float32)))
    points = points.to(surface.points.device)

    keywords = {"use_sign_winding_number": use_sign_winding_number}
    if max_dist is not None:
        keywords["max_dist"] = float(max_dist)
    distance, _, _ = signed_distance_field(surface, points, **keywords)
    return distance.detach().cpu().to(torch.float64).numpy()


def ComputeNodalSignedDistance(model_part: Kratos.ModelPart, max_dist=None,
                               use_sign_winding_number: bool = False,
                               source_container: str = "Elements"):
    """Signed distance at the model part's own nodes.

    Nodes lie ON the boundary surface they are measured against, so boundary
    nodes come out at ~0 and interior nodes negative - the usual "distance to
    the wall" feature.
    """
    from KratosMultiphysics.PhysicsNeMoApplication.bridges import graph_bridge
    mesh, _ = domain_mesh_builder.BuildMesh(
        model_part, (), source_container=source_container)
    return ComputeSignedDistance(
        mesh, graph_bridge.NodePositions(model_part), max_dist, use_sign_winding_number)


def WriteSignedDistanceField(model_part: Kratos.ModelPart, settings: Kratos.Parameters):
    """Computes the nodal SDF and stores it in a Kratos variable.

    Once written, the field is an ordinary nodal variable, so every existing
    gather picks it up through its normal "input_fields"/field_specs settings
    - no plumbing changes needed to feed geometry awareness into a grid,
    graph or point-cloud model.

    Args:
        model_part: The model part to measure.
        settings: Kratos Parameters:
            {
                "output_variable"          : "DISTANCE",
                "output_location"          : "node_non_historical",
                "source_container"         : "Elements",
                "max_dist"                 : 0.0,     // 0 = unbounded
                "use_sign_winding_number"  : false
            }

    Returns:
        The (N,) float64 array that was written.
    """
    defaults = Kratos.Parameters("""{
        "output_variable"         : "DISTANCE",
        "output_location"         : "node_non_historical",
        "source_container"        : "Elements",
        "max_dist"                : 0.0,
        "use_sign_winding_number" : false
    }""")
    settings.ValidateAndAssignDefaults(defaults)

    max_dist = settings["max_dist"].GetDouble()
    distance = ComputeNodalSignedDistance(
        model_part,
        max_dist=(max_dist if max_dist > 0.0 else None),
        use_sign_winding_number=settings["use_sign_winding_number"].GetBool(),
        source_container=settings["source_container"].GetString())

    variable = Kratos.KratosGlobals.GetVariable(settings["output_variable"].GetString())
    tensor_adaptor = GetTensorAdaptor(
        model_part, settings["output_location"].GetString(), variable, collect=True)
    tensor_adaptor.data[:] = distance.reshape(tensor_adaptor.data.shape)
    tensor_adaptor.StoreData()
    return distance


def SampleSignedDistanceOnGrid(model_part: Kratos.ModelPart, grid_shape,
                               bounding_box=None, max_dist=None,
                               use_sign_winding_number: bool = False,
                               source_container: str = "Elements"):
    """SDF sampled on the same lattice grid_bridge.SampleFieldsOnGrid builds.

    The grid case cannot go through a nodal variable: lattice points are not
    nodes. The returned channel concatenates directly onto a sampled grid.

    Returns:
        ((1, *grid_shape) float64 array, bounding_box) - the extra channel
        and the box it was sampled over.
    """
    from KratosMultiphysics.PhysicsNeMoApplication.bridges import grid_bridge
    grid_shape = tuple(int(n) for n in grid_shape)
    if bounding_box is None:
        bounding_box = grid_bridge.ComputeBoundingBox(model_part)
    points = grid_bridge._GridPointCoordinates(grid_shape, bounding_box)

    mesh, _ = domain_mesh_builder.BuildMesh(
        model_part, (), source_container=source_container)
    distance = ComputeSignedDistance(mesh, points, max_dist, use_sign_winding_number)
    return distance.reshape((1,) + grid_shape), bounding_box
