"""GLOBE: boundary-driven field prediction at arbitrary points.

physicsnemo.experimental.models.globe.GLOBE predicts fields at arbitrary
points from BOUNDARY data alone, with Green's-function-like kernels
evaluated on a dual-tree cluster hierarchy. That is the shape of an
elliptic problem - the boundary conditions determine the interior - and it
is exactly the data BuildDomainMesh already produces from a model part's
sub-model-parts.

What this bridge supplies is the translation:

- BuildGlobeBoundaryMeshes turns named sub-model-parts into the
  ``{name: Mesh}`` dict GLOBE's forward takes, with the boundary values
  living in ``cell_data`` where GLOBE reads them.
- RunGlobeForward calls the model and flattens the returned Mesh's
  ``point_data`` back into the (N, C_out) layout every writer here uses.

VERIFIED ON CPU. GLOBE's cluster tree accepts ``tree_build_device="cpu"``
and both the forward and the backward run there (a 2-D case forwards in
well under a second), so this is not a GPU-only capability.

TRAP: a rank-0 boundary field must be shaped ``(n_cells,)``, NOT
``(n_cells, 1)``. Upstream validates the declared rank spec against the
array's own rank and rejects the column form with a message that names the
field but not the shape.

physicsnemo is imported lazily; the model itself lives in
physicsnemo.experimental, so it carries no API stability guarantee.
"""

import numpy

import KratosMultiphysics as Kratos


def _TryImportTorch():
    try:
        import torch
        return torch
    except ImportError as e:
        raise ImportError(
            "PhysicsNeMoApplication.globe_bridge requires torch, which could not be "
            "imported. Install it with e.g. 'pip install torch'.") from e


def _TryImportGlobe():
    import warnings

    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")  # physicsnemo.experimental warns on import
            from physicsnemo.experimental.models.globe import GLOBE
        return GLOBE
    except ImportError as e:
        raise ImportError(
            "GLOBE requires physicsnemo >= 2.2 (its experimental models), which could "
            "not be imported. Install it with e.g. "
            "'pip install -U nvidia-physicsnemo'.") from e


def _PointDataToCellData(mesh, name, rank: int):
    """Averages a mesh's nodal field onto its cells, at the declared rank.

    GLOBE reads boundary sources from cell_data; the mesh bridge carries
    nodal fields as point_data. A rank-0 field must come out as (n_cells,)
    - the column form (n_cells, 1) is rejected upstream.
    """
    torch = _TryImportTorch()
    values = mesh.point_data[name]
    cells = mesh.cells
    gathered = values[cells.reshape(-1)]
    gathered = gathered.reshape(cells.shape[0], cells.shape[1], *values.shape[1:])
    averaged = gathered.mean(dim=1)
    if rank == 0:
        averaged = averaged.reshape(cells.shape[0])
    return averaged.to(torch.float32)


def BuildGlobeBoundaryMeshes(model_part: Kratos.ModelPart, sub_model_part_names,
                             field_specs, source_container: str = "Elements",
                             ranks=None):
    """The ``{name: Mesh}`` dict GLOBE's forward takes.

    Args:
        model_part: The parent part.
        sub_model_part_names: The named boundaries, e.g. ["Inlet", "Wall"].
        field_specs: [(variable_name, data_location)] nodal specs carrying
            the boundary values, the same shape every gather here uses.
        source_container: Which container the interior is tessellated from.
        ranks: Optional {variable_name: rank} (0 scalar, 1 vector);
            defaults to 0 for every field.

    Returns:
        {sub_model_part_name: physicsnemo Mesh} with the fields in cell_data.
    """
    from KratosMultiphysics.PhysicsNeMoApplication.bridges.mesh_bridge import (
        domain_mesh_builder)

    ranks = ranks or {}
    # the mesh bridge's field specs carry Kratos Variable OBJECTS, unlike
    # grid_bridge's, which carry variable NAMES - accept names here and
    # convert, since every process-facing block in this application is
    # written with names
    mesh_field_specs = [
        (Kratos.KratosGlobals.GetVariable(variable_name)
         if isinstance(variable_name, str) else variable_name, data_location)
        for variable_name, data_location in field_specs]
    for name in sub_model_part_names:
        if not model_part.HasSubModelPart(name):
            available = sorted(part.Name for part in model_part.SubModelParts)
            raise ValueError(
                f"Model part \"{model_part.Name}\" has no sub-model-part "
                f"\"{name}\"; it has {available or 'none'}. GLOBE's boundaries are "
                "named sub-model-parts carrying conditions.")

    domain_mesh, _ = domain_mesh_builder.BuildDomainMesh(
        model_part, field_specs=mesh_field_specs,
        boundary_sub_model_part_names=tuple(sub_model_part_names),
        source_container=source_container)

    boundary_meshes = {}
    for name in sub_model_part_names:
        if name not in domain_mesh.boundaries:
            raise ValueError(
                f"Sub-model-part \"{name}\" produced no boundary mesh; GLOBE's boundary "
                "meshes come from sub-model-parts that carry conditions (or elements).")
        mesh = domain_mesh.boundaries[name]
        for variable, _ in mesh_field_specs:
            field_name = variable.Name()
            mesh.cell_data[field_name] = _PointDataToCellData(
                mesh, field_name, ranks.get(field_name, 0))
        boundary_meshes[name] = mesh
    return boundary_meshes


def RunGlobeForward(model, coordinates, boundary_meshes, reference_lengths,
                    output_names, global_data=None, enable_grad: bool = False):
    """One GLOBE forward pass, flattened to the (N, C_out) layout.

    Args:
        model: A GLOBE instance.
        coordinates: (N, D) query coordinates (torch tensor or array-like).
        boundary_meshes: BuildGlobeBoundaryMeshes output.
        reference_lengths: {name: float or tensor}, matching the model's
            reference_length_names.
        output_names: The model's output field names, in the order the
            process's output fields expect them.
        global_data: Optional TensorDict of case-level data.
        enable_grad: Keep the autograd graph (training).

    Returns:
        (N, C_out) float64 torch tensor.
    """
    torch = _TryImportTorch()

    parameter = next(model.parameters(), None)
    dtype = parameter.dtype if parameter is not None else torch.float32
    points = torch.as_tensor(numpy.asarray(coordinates)).to(dtype)
    lengths = {
        name: (value if torch.is_tensor(value)
               else torch.tensor(float(value), dtype=dtype))
        for name, value in reference_lengths.items()
    }
    # The boundary meshes need the same normalization as the query points,
    # and a Kratos-built one arrives MIXED: the mesh bridge carries its
    # points as Kratos doubles while _PointDataToCellData has already cast
    # the cell data to float32. GLOBE's kernels are float32, so the points
    # alone are enough to fail deep inside a Linear with "mat1 and mat2 must
    # have the same dtype" - far from anything naming the mesh. Promoting
    # the model to float64 instead does not help: its tree scatter is hard
    # float32 and raises from index_add_. Mesh.to returns a NEW mesh, so the
    # caller's meshes (and the inference process's per-node-count cache) are
    # left untouched.
    boundary_meshes = {
        name: (mesh.to(dtype) if hasattr(mesh, "to") else mesh)
        for name, mesh in boundary_meshes.items()
    }

    context = torch.enable_grad() if enable_grad else torch.no_grad()
    with context:
        predicted = model(points, boundary_meshes, lengths, global_data)

    columns = []
    for name in output_names:
        if name not in predicted.point_data:
            raise ValueError(
                f"GLOBE returned no output field \"{name}\"; it predicts "
                f"{sorted(predicted.point_data.keys())}, named by its "
                "output_field_ranks.")
        values = predicted.point_data[name]
        columns.append(values.reshape(values.shape[0], -1))
    return torch.cat(columns, dim=-1).to(torch.float64)
