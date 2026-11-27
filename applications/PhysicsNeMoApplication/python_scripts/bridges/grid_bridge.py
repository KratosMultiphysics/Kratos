"""Bridge between unstructured Kratos meshes and regular voxel grids.

Superresolution (and other image-style) models consume regular
(C, D, H, W) grids, while Kratos fields live on unstructured meshes. This
module samples nodal fields onto a lattice using Kratos's fast point locator
(FE shape-function interpolation — exact for fields the element interpolates
exactly, e.g. linear fields on simplex meshes) and scatters grids back onto
nodes with trilinear interpolation.

Pure Kratos + numpy, with one exception: ComputeGridDerivatives lazily
imports torch/physicsnemo to expose the physics-consistent grid-derivative
operators (everything else works without them).
"""

import numpy

import KratosMultiphysics as Kratos
from KratosMultiphysics.PhysicsNeMoApplication.utilities import array_backend_utils
from KratosMultiphysics.PhysicsNeMoApplication.utilities.tensor_adaptor_dataset_utils import (
    GetTensorAdaptor, RowsOfIds)

_SIMPLEX_TYPES_3D = (
    Kratos.GeometryData.KratosGeometryType.Kratos_Tetrahedra3D4,
)

_NODAL_LOCATIONS = ("node_historical", "node_non_historical")


def ComputeBoundingBox(model_part: Kratos.ModelPart, padding: float = 0.0):
    """Axis-aligned bounding box of the model part's nodes as (low, high)."""
    position_ta = Kratos.TensorAdaptors.NodePositionTensorAdaptor(model_part.Nodes, Kratos.Configuration.Current)
    position_ta.CollectData()
    coordinates = numpy.array(position_ta.data)
    return coordinates.min(axis=0) - padding, coordinates.max(axis=0) + padding


def _GridPointCoordinates(grid_shape, bounding_box):
    low, high = (numpy.asarray(b, dtype=float) for b in bounding_box)
    axes = [numpy.linspace(low[i], high[i], grid_shape[i]) for i in range(3)]
    mesh = numpy.meshgrid(*axes, indexing="ij")
    return numpy.stack([m.ravel() for m in mesh], axis=1)  # (D*H*W, 3)


def _GatherNodalField(model_part, variable_name, data_location):
    if data_location not in _NODAL_LOCATIONS:
        raise ValueError(
            f"Grid sampling supports nodal locations only ({', '.join(_NODAL_LOCATIONS)}), "
            f"got \"{data_location}\".")
    variable = Kratos.KratosGlobals.GetVariable(variable_name)
    tensor_adaptor = GetTensorAdaptor(model_part, data_location, variable)
    data = numpy.array(tensor_adaptor.data)
    return data.reshape(data.shape[0], -1)  # (n_nodes, width)


def SampleFieldsOnGrid(model_part: Kratos.ModelPart,
                       field_specs,
                       grid_shape,
                       bounding_box=None,
                       fill_value: float = 0.0,
                       backend: str = "kratos"):
    """Samples nodal fields of a model part onto a regular voxel grid.

    Args:
        model_part: The source model part (3D).
        field_specs: iterable of (variable_name, data_location) pairs; nodal
            locations only.
        grid_shape: (D, H, W) number of grid points per axis.
        bounding_box: Optional (low, high) arrays; defaults to the model
            part's node bounding box.
        fill_value: Value assigned outside the mesh.
        backend: "kratos" (default) uses BinBasedFastPointLocator3D - the
            vectorized path for all-simplex meshes, a per-point fallback
            for general geometries, fill_value outside the mesh.
            "physicsnemo" tessellates and samples over the mesh's own BVH
            (mesh_bridge.sampling), which is GPU-resident and
            autograd-friendly. The default stays "kratos" because the
            benchmark, not an assumption, decides which is faster - see
            benchmarks/benchmark_bridges.py --sampling both.

    Returns:
        (grid, bounding_box): grid has shape (C, D, H, W) float64, where C is
        the total flattened per-node width across the fields (the same
        channel convention InferenceProcess uses).
    """
    grid_shape = tuple(int(n) for n in grid_shape)
    if backend not in ("kratos", "physicsnemo"):
        raise ValueError(
            f"Unsupported sampling backend \"{backend}\". Use \"kratos\" (the point "
            "locator) or \"physicsnemo\" (the mesh's own BVH).")
    if backend == "physicsnemo":
        return _SampleFieldsOnGridWithBvh(
            model_part, field_specs, grid_shape, bounding_box, fill_value)
    if len(grid_shape) != 3 or any(n < 2 for n in grid_shape):
        raise ValueError(f"grid_shape must be three axis sizes >= 2, got {grid_shape}.")
    if bounding_box is None:
        bounding_box = ComputeBoundingBox(model_part)

    points = _GridPointCoordinates(grid_shape, bounding_box)

    part_ids = numpy.fromiter((node.Id for node in model_part.Nodes),
                              dtype=numpy.int64, count=model_part.NumberOfNodes())
    id_order = numpy.argsort(part_ids, kind="stable")
    sorted_part_ids = part_ids[id_order]
    fields = [_GatherNodalField(model_part, name, location) for name, location in field_specs]
    widths = [field.shape[1] for field in fields]
    stacked = numpy.concatenate(fields, axis=1)  # (n_nodes, C)
    total_width = stacked.shape[1]

    locator = Kratos.BinBasedFastPointLocator3D(model_part)
    locator.UpdateSearchDatabase()

    values = numpy.full((points.shape[0], total_width), float(fill_value))
    all_simplex = all(
        element.GetGeometry().GetGeometryType() in _SIMPLEX_TYPES_3D
        for element in model_part.Elements)

    if all_simplex:
        element_ids, node_ids, shape_values = locator.VectorizedFind(Kratos.Matrix(points))
        element_ids = numpy.asarray(element_ids)
        node_ids = numpy.asarray(node_ids)
        shape_values = numpy.asarray(shape_values)
        found = element_ids != -1  # -1 is the locator's not-found sentinel
        if found.any():
            # searchsorted rather than numpy.vectorize, which is a Python
            # loop over every located point's corner ids.
            rows = id_order[numpy.searchsorted(
                sorted_part_ids, node_ids[found].astype(numpy.int64))]  # (n_found, 4)
            values[found] = numpy.einsum("pk,pkc->pc", shape_values[found], stacked[rows])
    else:
        # per-point locator for general geometries (hexahedra, prisms, ...);
        # the corner rows come from the same sorted-id lookup as above, not
        # from a per-call {id: row} dict
        for i, point in enumerate(points):
            is_found, shape_functions, element = locator.FindPointOnMesh(point)
            if not is_found:
                continue
            corner_ids = numpy.fromiter((node.Id for node in element.GetGeometry()),
                                        dtype=numpy.int64)
            rows = id_order[numpy.searchsorted(sorted_part_ids, corner_ids)]
            values[i] = numpy.asarray(shape_functions) @ stacked[rows]

    grid = values.T.reshape((total_width,) + grid_shape)
    return grid, bounding_box


def InterpolateGridAtPoints(grid, bounding_box, points, backend="numpy"):
    """Trilinear interpolation of a (C, D, H, W) grid at arbitrary points.

    Points are clamped to the bounding box. Exact for linear fields.

    Args:
        backend: "numpy" (default), "cupy" or "auto". Eight strided gathers
            of the whole grid is the shape of work a GPU is built for: the
            CuPy path measured ~14.6x on a 3-channel 64^3 grid at 100k
            points, transfers included. It loses below the size threshold
            (~0.6x on a 1-channel 48^3 grid at 15k points), where it falls
            back to numpy on its own.

    Returns:
        (n_points, C) float64 array.
    """
    n_points = len(numpy.asarray(points))
    xp, _ = array_backend_utils.ResolveArrayModule(
        backend, size_hint=n_points * numpy.asarray(grid).shape[0])
    grid = xp.asarray(grid)
    low, high = (xp.asarray(numpy.asarray(b, dtype=float)) for b in bounding_box)
    shape = xp.asarray(numpy.array(grid.shape[1:]))

    # Fractional lattice coordinates, clamped inside the grid. A FLAT axis
    # (every node at the same Z, i.e. any planar Kratos case) has a zero
    # extent: dividing by it gives NaN, and NaN cast to an index is INT_MIN,
    # which used to crash here. Such an axis carries identical slices, so
    # the whole weight belongs on index 0.
    extent = xp.where(high - low == 0.0, 1.0, high - low)
    fractional = (xp.asarray(numpy.asarray(points, dtype=float)) - low) / extent * (shape - 1)
    fractional = xp.clip(fractional, 0.0, shape - 1)
    base = xp.minimum(fractional.astype(int), shape - 2)
    t = fractional - base  # (n, 3) in [0, 1]

    result = xp.zeros((len(fractional), grid.shape[0]))
    for dx in (0, 1):
        for dy in (0, 1):
            for dz in (0, 1):
                weight = ((t[:, 0] if dx else 1.0 - t[:, 0]) *
                          (t[:, 1] if dy else 1.0 - t[:, 1]) *
                          (t[:, 2] if dz else 1.0 - t[:, 2]))
                corner = grid[:, base[:, 0] + dx, base[:, 1] + dy, base[:, 2] + dz]  # (C, n)
                result += weight[:, None] * corner.T
    return array_backend_utils.ToHost(result)


def _SampleFieldsOnGridWithBvh(model_part, field_specs, grid_shape, bounding_box,
                               fill_value):
    """SampleFieldsOnGrid's physicsnemo backend: tessellate, then sample.

    Same (C, D, H, W) contract as the locator path. Points outside the mesh
    come back NaN from upstream rather than as a sentinel, so they are
    replaced with fill_value here to keep the two backends interchangeable.
    """
    from KratosMultiphysics.PhysicsNeMoApplication.bridges.mesh_bridge import sampling

    if bounding_box is None:
        bounding_box = ComputeBoundingBox(model_part)
    points = _GridPointCoordinates(grid_shape, bounding_box)
    values, _ = sampling.SampleModelPartAtPoints(model_part, field_specs, points)
    values = numpy.asarray(values.detach().cpu().numpy(), dtype=numpy.float64)
    values = numpy.nan_to_num(values, nan=fill_value)
    grid = values.T.reshape((values.shape[1],) + tuple(grid_shape))
    return grid, bounding_box


def _TryImportTorch():
    try:
        import torch
        return torch
    except ImportError as e:
        raise ImportError(
            "PhysicsNeMoApplication.grid_bridge's differentiable interpolation requires "
            "torch, which could not be imported. Install it with e.g. "
            "'pip install torch'.") from e


def InterpolateGridAtPointsTorch(grid, bounding_box, points):
    """Differentiable trilinear interpolation of a (C, D, H, W) grid.

    The numpy/cupy InterpolateGridAtPoints is the fast read-out path and
    breaks the autograd graph. This one is the same arithmetic in torch, so
    a gradient reaches the GRID - which is what a diffusion sample graded
    by a Kratos residual needs (see physics.diffusion_residual_operator).
    Pinned against the numpy version value for value.

    Args:
        grid: (C, D, H, W) torch tensor (any dtype/device; gradients flow).
        bounding_box: (low, high) triples of the grid's physical extent.
        points: (n_points, 3) array-like query coordinates.

    Returns:
        (n_points, C) torch tensor on the grid's device and dtype.
    """
    torch = _TryImportTorch()
    low, high = (torch.as_tensor(numpy.asarray(bound, dtype=float),
                                 device=grid.device, dtype=grid.dtype)
                 for bound in bounding_box)
    shape = torch.as_tensor(numpy.array(grid.shape[1:], dtype=float),
                            device=grid.device, dtype=grid.dtype)

    query = torch.as_tensor(numpy.asarray(points, dtype=float),
                            device=grid.device, dtype=grid.dtype)
    # a flat axis (a planar case) has zero extent - see InterpolateGridAtPoints
    extent = torch.where(high - low == 0.0, torch.ones_like(high), high - low)
    fractional = (query - low) / extent * (shape - 1.0)
    fractional = torch.clamp(fractional, torch.zeros_like(shape), shape - 1.0)
    base = torch.minimum(fractional.to(torch.int64), (shape - 2.0).to(torch.int64))
    t = fractional - base.to(grid.dtype)  # (n, 3) in [0, 1]

    result = torch.zeros((query.shape[0], grid.shape[0]),
                         device=grid.device, dtype=grid.dtype)
    for dx in (0, 1):
        for dy in (0, 1):
            for dz in (0, 1):
                weight = ((t[:, 0] if dx else 1.0 - t[:, 0]) *
                          (t[:, 1] if dy else 1.0 - t[:, 1]) *
                          (t[:, 2] if dz else 1.0 - t[:, 2]))
                corner = grid[:, base[:, 0] + dx, base[:, 1] + dy, base[:, 2] + dz]  # (C, n)
                result = result + weight[:, None] * corner.T
    return result


def ScatterGridToNodes(grid, bounding_box, model_part: Kratos.ModelPart, output_field_specs,
                       backend="numpy") -> None:
    """Writes a (C, D, H, W) grid onto a model part's nodes.

    Channels are split across the output fields by each field's per-node
    width, interpolated trilinearly at the node positions, and stored via the
    tensor adaptors.

    Args:
        backend: Passed through to InterpolateGridAtPoints.
    """
    position_ta = Kratos.TensorAdaptors.NodePositionTensorAdaptor(model_part.Nodes, Kratos.Configuration.Current)
    position_ta.CollectData()
    nodal_values = InterpolateGridAtPoints(grid, bounding_box, numpy.array(position_ta.data),
                                           backend=backend)

    offset = 0
    for variable_name, data_location in output_field_specs:
        if data_location not in _NODAL_LOCATIONS:
            raise ValueError(
                f"Grid scatter supports nodal locations only ({', '.join(_NODAL_LOCATIONS)}), "
                f"got \"{data_location}\".")
        variable = Kratos.KratosGlobals.GetVariable(variable_name)
        tensor_adaptor = GetTensorAdaptor(model_part, data_location, variable)
        width = int(numpy.prod(tensor_adaptor.data.shape[1:], dtype=int))
        chunk = nodal_values[:, offset:offset + width].reshape(tensor_adaptor.data.shape)
        tensor_adaptor.data[:] = chunk
        tensor_adaptor.StoreData()
        offset += width
    if offset != grid.shape[0]:
        raise ValueError(
            f"Grid has {grid.shape[0]} channels but the output fields consume {offset}.")


def _TryImportGridDerivatives():
    try:
        from physicsnemo.nn.functional import derivatives
        import torch
        return torch, derivatives
    except ImportError as e:
        raise ImportError(
            "PhysicsNeMoApplication.grid_bridge.ComputeGridDerivatives requires torch and "
            "physicsnemo, which could not be imported. Install them with e.g. "
            "'pip install torch nvidia-physicsnemo'.") from e


def ComputeGridDerivatives(grid, settings):
    """Physics-consistent spatial derivatives of a (C, *spatial) grid.

    Wraps physicsnemo.nn.functional's uniform/rectilinear/spectral grid
    gradients, which operate on bare 1D-3D scalar fields with PERIODIC
    stencils only. This wrapper adds the channel loop and, for non-periodic
    data, boundary trimming - upstream has no one-sided stencils, so on
    non-periodic fields the outermost layer of every axis wraps around and
    is garbage (the interior is exact for polynomials at the stencil order).
    Spectral derivatives are only meaningful on truly periodic data.

    Args:
        grid: (C, *spatial) array/tensor, 1-3 spatial axes - the layout
            SampleFieldsOnGrid produces. Differentiable when a torch tensor
            with requires_grad is passed.
        settings: Kratos Parameters:
            {
                "operator"          : "uniform" | "rectilinear" | "spectral",
                "spacing"           : [1.0],      // uniform: one per axis (or one value for all)
                "coordinates"       : [[...]],    // rectilinear: per-axis coordinate arrays
                "lengths"           : [1.0],      // spectral: domain lengths per axis
                "order"             : 2,          // uniform only: 2 or 4
                "derivative_orders" : [1],        // e.g. [1, 2]
                "include_mixed"     : false,
                "boundary"          : "periodic"  // or "trim"
            }

    Returns:
        (result, interior_slices): result is (C, num_derivatives, *spatial')
        - first derivatives in axis order, then second derivatives, then
        mixed terms in axis-pair order (x,y), (x,z), (y,z). With
        boundary="trim", spatial' is the interior (each axis loses its first
        and last `trim_width` layers) and interior_slices is the per-axis
        slice tuple that maps the original grid onto it; with "periodic" the
        shape is unchanged and interior_slices covers everything.
    """
    torch, derivatives = _TryImportGridDerivatives()

    defaults = Kratos.Parameters("""{
        "operator"          : "uniform",
        "spacing"           : [1.0],
        "coordinates"       : [],
        "lengths"           : [1.0],
        "order"             : 2,
        "derivative_orders" : [1],
        "include_mixed"     : false,
        "boundary"          : "periodic"
    }""")
    settings.ValidateAndAssignDefaults(defaults)
    operator = settings["operator"].GetString()
    boundary = settings["boundary"].GetString()
    if boundary not in ("periodic", "trim"):
        raise ValueError(f"Unknown boundary mode \"{boundary}\". Use \"periodic\" or \"trim\".")

    grid = torch.as_tensor(grid)
    if grid.dim() < 2 or grid.dim() > 4:
        raise ValueError(
            f"grid must be (C, *spatial) with 1-3 spatial axes, got shape {tuple(grid.shape)}.")
    n_spatial = grid.dim() - 1

    derivative_orders = tuple(settings["derivative_orders"].GetVector())
    derivative_orders = tuple(int(order) for order in derivative_orders)
    include_mixed = settings["include_mixed"].GetBool()

    def per_axis(vector, name):
        values = list(vector)
        if len(values) == 1:
            values = values * n_spatial
        if len(values) != n_spatial:
            raise ValueError(
                f"\"{name}\" needs one entry per spatial axis ({n_spatial}), got {len(values)}.")
        return values

    if operator == "uniform":
        spacing = per_axis(settings["spacing"].GetVector(), "spacing")
        order = settings["order"].GetInt()
        compute = lambda field: derivatives.uniform_grid_gradient(
            field, spacing=spacing, order=order,
            derivative_orders=derivative_orders, include_mixed=include_mixed)
        trim_width = order // 2
    elif operator == "rectilinear":
        if settings["coordinates"].size() != n_spatial:
            raise ValueError(
                f"\"coordinates\" needs one array per spatial axis ({n_spatial}).")
        coordinates = [torch.as_tensor(numpy.array(settings["coordinates"][i].GetVector()),
                                       dtype=grid.dtype, device=grid.device)
                       for i in range(n_spatial)]
        compute = lambda field: derivatives.rectilinear_grid_gradient(
            field, coordinates=coordinates,
            derivative_orders=derivative_orders, include_mixed=include_mixed)
        trim_width = 1
    elif operator == "spectral":
        if boundary == "trim":
            raise ValueError(
                "Spectral derivatives are global: trimming cannot repair non-periodic "
                "data. Use operator \"uniform\" or \"rectilinear\" with boundary \"trim\".")
        lengths = per_axis(settings["lengths"].GetVector(), "lengths")
        compute = lambda field: derivatives.spectral_grid_gradient(
            field, lengths=lengths,
            derivative_orders=derivative_orders, include_mixed=include_mixed)
        trim_width = 0
    else:
        raise ValueError(
            f"Unknown operator \"{operator}\". Use \"uniform\", \"rectilinear\" or \"spectral\".")

    result = torch.stack([compute(grid[c]) for c in range(grid.shape[0])], dim=0)

    interior = tuple(slice(None) for _ in range(n_spatial))
    if boundary == "trim":
        if any(size <= 2 * trim_width for size in grid.shape[1:]):
            raise ValueError(
                f"Grid spatial shape {tuple(grid.shape[1:])} is too small to trim "
                f"{trim_width} boundary layers per side.")
        interior = tuple(slice(trim_width, -trim_width) for _ in range(n_spatial))
        result = result[(slice(None), slice(None)) + interior]
    return result, interior


def _TryImportGridVectorOperators():
    try:
        from physicsnemo.nn.functional import derivatives
        import torch
        for name in ("uniform_grid_divergence", "uniform_grid_curl", "uniform_grid_laplacian"):
            if not hasattr(derivatives, name):
                raise ImportError(
                    f"physicsnemo.nn.functional.derivatives has no {name}; the grid "
                    "divergence/curl/laplacian operators require physicsnemo >= 2.2.")
        return torch, derivatives
    except ImportError as e:
        raise ImportError(
            "PhysicsNeMoApplication.grid_bridge.ComputeGridVectorOperator requires torch "
            "and physicsnemo >= 2.2, which could not be imported. Install them with e.g. "
            "'pip install torch nvidia-physicsnemo'.") from e


_GRID_OPERATIONS = ("divergence", "curl", "laplacian")


def ComputeGridVectorOperator(grid, settings):
    """Divergence, curl or Laplacian of a grid field (physicsnemo >= 2.2).

    The vector counterpart of ComputeGridDerivatives, sharing its operator
    and boundary conventions. Note the two upstream families disagree on
    layout, which this wrapper hides: the gradients take a bare SCALAR field
    and prepend a derivative axis, while divergence/curl take a CHANNEL-FIRST
    vector field whose channel count must equal the number of spatial axes.

    Args:
        grid: For "divergence"/"curl", a (D, *spatial) vector field with
            D == len(spatial) - the layout SampleFieldsOnGrid already
            produces for a D-component field. For "laplacian", either
            (*spatial) or (C, *spatial), applied per channel.
        settings: Kratos Parameters:
            {
                "operation" : "divergence" | "curl" | "laplacian",
                "operator"  : "uniform" | "rectilinear",
                "spacing"   : [1.0],       // uniform
                "coordinates" : [[...]],   // rectilinear
                "order"     : 2,           // uniform only (2 or 4)
                "boundary"  : "periodic" | "trim"
            }

    Returns:
        (result, interior_slices). Shapes: divergence -> (*spatial);
        curl -> (*spatial) in 2D (scalar vorticity) and (3, *spatial) in 3D;
        laplacian -> the input's shape. With boundary="trim" the outermost
        layers are cropped (the stencils are periodic-only, exactly as for
        the gradients) and interior_slices maps the original grid onto the
        result.
    """
    torch, derivatives = _TryImportGridVectorOperators()

    defaults = Kratos.Parameters("""{
        "operation"       : "divergence",
        "operator"        : "uniform",
        "spacing"         : [1.0],
        "coordinates"     : [],
        "order"           : 2,
        "boundary"        : "periodic",
        "has_channel_axis": false,
        "implementation"  : "auto"
    }""")
    settings.ValidateAndAssignDefaults(defaults)

    operation = settings["operation"].GetString()
    if operation not in _GRID_OPERATIONS:
        raise ValueError(
            f"Unknown operation \"{operation}\". Use one of {_GRID_OPERATIONS}.")
    operator = settings["operator"].GetString()
    if operator not in ("uniform", "rectilinear"):
        raise ValueError(
            f"Unknown operator \"{operator}\". Use \"uniform\" or \"rectilinear\" "
            "(spectral divergence/curl/laplacian do not exist upstream).")
    boundary = settings["boundary"].GetString()
    if boundary not in ("periodic", "trim"):
        raise ValueError(f"Unknown boundary mode \"{boundary}\". Use \"periodic\" or \"trim\".")

    grid = torch.as_tensor(grid)
    # (C, *spatial) and (*spatial) are indistinguishable from the shape alone
    # for a laplacian, so the channel axis is declared rather than guessed
    has_channel_axis = (operation in ("divergence", "curl")
                        or settings["has_channel_axis"].GetBool())
    n_spatial = grid.dim() - 1 if has_channel_axis else grid.dim()
    if operation in ("divergence", "curl"):
        if int(grid.shape[0]) != n_spatial:
            raise ValueError(
                f"\"{operation}\" needs a channel-first vector field whose channel count "
                f"equals the number of spatial axes; got {tuple(grid.shape)}.")
    if not 1 <= n_spatial <= 3:
        raise ValueError(f"grid must have 1-3 spatial axes, got shape {tuple(grid.shape)}.")

    def per_axis(vector, name):
        values = list(vector)
        if len(values) == 1:
            values = values * n_spatial
        if len(values) != n_spatial:
            raise ValueError(
                f"\"{name}\" needs one entry per spatial axis ({n_spatial}), got {len(values)}.")
        return values

    # The Warp backend computes in float32 and is picked automatically when a
    # CUDA device exists, silently costing ~7 digits on float64 input; keep the
    # torch backend for float64 unless the caller overrides it.
    implementation = settings["implementation"].GetString()
    if implementation == "auto":
        implementation = "torch" if grid.dtype == torch.float64 else None
    elif implementation == "default":
        implementation = None

    order = settings["order"].GetInt()
    if operator == "uniform":
        keywords = {"spacing": per_axis(settings["spacing"].GetVector(), "spacing"),
                    "order": order, "implementation": implementation}
        trim_width = order // 2
        table = {"divergence": derivatives.uniform_grid_divergence,
                 "curl": derivatives.uniform_grid_curl,
                 "laplacian": derivatives.uniform_grid_laplacian}
    else:
        if settings["coordinates"].size() != n_spatial:
            raise ValueError(
                f"\"coordinates\" needs one array per spatial axis ({n_spatial}).")
        keywords = {"coordinates": [
            torch.as_tensor(numpy.array(settings["coordinates"][i].GetVector()),
                            dtype=grid.dtype, device=grid.device)
            for i in range(n_spatial)], "implementation": implementation}
        trim_width = 1
        table = {"divergence": derivatives.rectilinear_grid_divergence,
                 "curl": derivatives.rectilinear_grid_curl,
                 "laplacian": derivatives.rectilinear_grid_laplacian}

    if operation == "laplacian" and has_channel_axis:
        result = torch.stack([table[operation](grid[c], **keywords)
                              for c in range(grid.shape[0])], dim=0)
    else:
        result = table[operation](grid, **keywords)

    interior = tuple(slice(None) for _ in range(n_spatial))
    if boundary == "trim":
        spatial_shape = grid.shape[-n_spatial:]
        if any(int(size) <= 2 * trim_width for size in spatial_shape):
            raise ValueError(
                f"Grid spatial shape {tuple(int(s) for s in spatial_shape)} is too small "
                f"to trim {trim_width} boundary layers per side.")
        interior = tuple(slice(trim_width, -trim_width) for _ in range(n_spatial))
        leading = result.dim() - n_spatial
        result = result[(slice(None),) * leading + interior]
    return result, interior
