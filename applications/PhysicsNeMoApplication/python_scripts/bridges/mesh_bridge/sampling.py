"""Torch-native point sampling over a mesh's own BVH.

grid_bridge samples fields onto a lattice with Kratos's
BinBasedFastPointLocator, and mapping_bridge transfers between meshes with
MappingApplication. physicsnemo 2.2 offers a third way:
``mesh.sampling.sample_data_at_points`` over a
``mesh.spatial.BVH``, which is GPU-resident and autograd-friendly where the
other two are neither.

Whether it is FASTER is a measurement, not an assumption, so the default
stays where it was: ``benchmarks/benchmark_bridges.py --sampling both``
times the two at the standard sizes and prints the verdict, and
grid_bridge.SampleFieldsOnGrid takes a ``backend`` argument rather than
switching underneath anyone.

What this path uniquely gives is the gradient: a value sampled at a point
is differentiable with respect to the mesh's own data, which a locator
lookup is not.

torch and physicsnemo are imported lazily.
"""

import numpy

import KratosMultiphysics as Kratos


def _TryImportTorch():
    try:
        import torch
        return torch
    except ImportError as e:
        raise ImportError(
            "PhysicsNeMoApplication.mesh_bridge.sampling requires torch, which could "
            "not be imported. Install it with e.g. 'pip install torch'.") from e


def _TryImportSampling():
    try:
        from physicsnemo.mesh.sampling import sample_data_at_points
        from physicsnemo.mesh.spatial import BVH
        return sample_data_at_points, BVH
    except ImportError as e:
        raise ImportError(
            "Torch-native mesh sampling requires physicsnemo >= 2.2, which could not "
            "be imported. Install it with e.g. 'pip install -U nvidia-physicsnemo'.") from e


def BuildBvh(mesh, leaf_size: int = 1):
    """The bounding-volume hierarchy of a mesh, reusable across queries.

    Build it once when sampling the same mesh repeatedly: it is the part
    that costs, and sample_data_at_points rebuilds it every call otherwise.
    """
    _, BVH = _TryImportSampling()
    return BVH.from_mesh(mesh, leaf_size=leaf_size)


def SampleMeshAtPoints(mesh, query_points, data_source: str = "points",
                       bvh=None, **options):
    """Samples a mesh's fields at arbitrary points.

    Args:
        mesh: A physicsnemo Mesh carrying point_data (or cell_data).
        query_points: (n_queries, n_spatial_dims) array-like.
        data_source: "points" (interpolate nodal data) or "cells".
        bvh: A BuildBvh result to reuse; built internally when None.
        options: Forwarded upstream - multiple_cells_strategy,
            project_onto_nearest_cell, tolerance.

    Returns:
        A TensorDict of the sampled fields, keyed as on the mesh.
    """
    sample_data_at_points, _ = _TryImportSampling()
    torch = _TryImportTorch()
    points = torch.as_tensor(numpy.asarray(query_points, dtype=numpy.float64))
    points = points.to(torch.as_tensor(mesh.points).dtype)
    return sample_data_at_points(
        mesh, points, data_source=data_source, bvh=bvh, **options)


def SampleModelPartAtPoints(model_part: Kratos.ModelPart, field_specs, query_points,
                            source_container: str = "Elements", bvh=None,
                            tessellation_mode: str = "smallest_id_diagonal"):
    """Tessellates a model part and samples its fields at arbitrary points.

    The mesh-bridge counterpart of grid_bridge's locator path, for callers
    who want the values as differentiable tensors rather than as a grid.

    Args:
        model_part: The model part.
        field_specs: [(variable_name, data_location)] nodal specs.
        query_points: (n_queries, 3) array-like.
        source_container: Which container to tessellate.
        bvh: An optional prebuilt BVH (see BuildBvh).
        tessellation_mode: Passed to BuildMesh.

    Returns:
        (values, mesh): the (n_queries, total_width) float64 tensor in the
        field specs' column order, and the tessellated mesh so a caller can
        keep its BVH.
    """
    from KratosMultiphysics.PhysicsNeMoApplication.bridges.mesh_bridge import (
        domain_mesh_builder)

    torch = _TryImportTorch()
    mesh_field_specs = [
        (Kratos.KratosGlobals.GetVariable(variable_name)
         if isinstance(variable_name, str) else variable_name, data_location)
        for variable_name, data_location in field_specs]
    mesh, _ = domain_mesh_builder.BuildMesh(
        model_part, field_specs=mesh_field_specs, source_container=source_container,
        tessellation_mode=tessellation_mode)

    sampled = SampleMeshAtPoints(mesh, query_points, data_source="points", bvh=bvh)
    columns = []
    for variable, _ in mesh_field_specs:
        values = sampled[variable.Name()]
        columns.append(values.reshape(values.shape[0], -1))
    return torch.cat(columns, dim=-1).to(torch.float64), mesh
