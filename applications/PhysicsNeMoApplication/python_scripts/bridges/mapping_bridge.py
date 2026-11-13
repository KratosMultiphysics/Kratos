"""Field transfer between non-matching meshes via MappingApplication.

The grid bridge (grid_bridge.py) interpolates through the FE shape functions
of the source mesh - exact, but tied to meshes the point locator handles and
to lattices inside the source domain. For ML grids (or auxiliary ML meshes)
that match no tessellation of the Kratos mesh, MappingApplication's mappers
do the transfer instead: nearest-neighbor/nearest-element/barycentric
searches that are robust to partial overlap, non-simplex geometries and MPI
distribution.

The typical pattern: CreateBackgroundGridModelPart builds a structured
"ML grid" model part over a bounding box, MappingBridge maps solver fields
onto it (and inversely maps predictions back), and GatherGridArray hands
the mapped fields to the grid-model stack as a (C, D, H, W) array.

MappingApplication is a compiled optional dependency: it is imported lazily
so this module can always be imported, mirroring the torch/physicsnemo
policy.
"""

import numpy

import KratosMultiphysics as Kratos


def _TryImportMappingApplication():
    try:
        import KratosMultiphysics.MappingApplication  # noqa: F401 - registers the mappers
        from KratosMultiphysics.MappingApplication import python_mapper_factory
        return python_mapper_factory
    except ImportError as e:
        raise ImportError(
            "PhysicsNeMoApplication.mapping_bridge requires MappingApplication, which is "
            "not available. Compile Kratos with "
            "'add_app ${KRATOS_APP_DIR}/MappingApplication'.") from e


def CreateBackgroundGridModelPart(model: Kratos.Model,
                                  name: str,
                                  bounding_box,
                                  divisions,
                                  dimension: int = 3,
                                  historical_variables=()):
    """Creates a structured background-grid model part over a bounding box.

    The "ML grid that matches no tessellation": a structured simplicial
    mesh whose nodes form the regular (divisions + 1)-per-axis lattice the
    mapped fields are read from (GatherGridArray). The elements only exist
    so mappers and locators can search the part - the lattice is the point.

    Args:
        model: The Kratos Model to create the part in.
        name: Name of the new model part.
        bounding_box: (low, high) arrays; in 2D only x/y are used.
        divisions: cells per axis (the lattice has divisions + 1 points per
            axis).
        dimension: 2 or 3.
        historical_variables: nodal solution-step variables to allocate
            (mapping targets).

    Returns:
        The new model part.
    """
    low, high = (numpy.asarray(b, dtype=float) for b in bounding_box)
    divisions = int(divisions)
    if divisions < 1:
        raise ValueError(f"divisions must be >= 1, got {divisions}.")
    if dimension not in (2, 3):
        raise ValueError(f"dimension must be 2 or 3, got {dimension}.")

    model_part = model.CreateModelPart(name)
    model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = dimension
    for variable in historical_variables:
        model_part.AddNodalSolutionStepVariable(variable)

    if dimension == 2:
        generator_geometry = Kratos.Quadrilateral2D4(
            Kratos.Node(1, low[0], low[1], 0.0), Kratos.Node(2, low[0], high[1], 0.0),
            Kratos.Node(3, high[0], high[1], 0.0), Kratos.Node(4, high[0], low[1], 0.0))
        element_name = "Element2D3N"
    else:
        generator_geometry = Kratos.Hexahedra3D8(
            Kratos.Node(1, low[0], low[1], low[2]), Kratos.Node(2, high[0], low[1], low[2]),
            Kratos.Node(3, high[0], high[1], low[2]), Kratos.Node(4, low[0], high[1], low[2]),
            Kratos.Node(5, low[0], low[1], high[2]), Kratos.Node(6, high[0], low[1], high[2]),
            Kratos.Node(7, high[0], high[1], high[2]), Kratos.Node(8, low[0], high[1], high[2]))
        element_name = "Element3D4N"

    mesh_parameters = Kratos.Parameters("""{
        "number_of_divisions"        : %d,
        "element_name"               : "%s",
        "condition_name"             : "LineCondition",
        "create_skin_sub_model_part" : false,
        "elements_properties_id"     : 1
    }""" % (divisions, element_name))
    Kratos.StructuredMeshGeneratorProcess(generator_geometry, model_part, mesh_parameters).Execute()
    return model_part


class MappingBridge:
    """Maps fields between an origin and a destination model part."""

    def __init__(self, origin_model_part: Kratos.ModelPart,
                 destination_model_part: Kratos.ModelPart,
                 settings: Kratos.Parameters = None) -> None:
        python_mapper_factory = _TryImportMappingApplication()

        if settings is None:
            settings = Kratos.Parameters("""{}""")
        default_settings = Kratos.Parameters("""{
            "mapper_type" : "nearest_element",
            "echo_level"  : 0
        }""")
        settings.ValidateAndAssignDefaults(default_settings)

        self.origin_model_part = origin_model_part
        self.destination_model_part = destination_model_part
        if origin_model_part.IsDistributed():
            self._mapper = python_mapper_factory.CreateMPIMapper(
                origin_model_part, destination_model_part, settings.Clone())
        else:
            self._mapper = python_mapper_factory.CreateMapper(
                origin_model_part, destination_model_part, settings.Clone())

    def MapFields(self, field_pairs) -> None:
        """Maps origin fields onto destination fields.

        Args:
            field_pairs: iterable of (origin_variable_name,
                destination_variable_name) pairs (historical nodal
                variables, the mappers' native data location).
        """
        for origin_name, destination_name in field_pairs:
            self._mapper.Map(
                Kratos.KratosGlobals.GetVariable(origin_name),
                Kratos.KratosGlobals.GetVariable(destination_name))

    def InverseMapFields(self, field_pairs) -> None:
        """Maps destination fields back onto origin fields (same pair order
        as MapFields: (origin_variable_name, destination_variable_name))."""
        for origin_name, destination_name in field_pairs:
            self._mapper.InverseMap(
                Kratos.KratosGlobals.GetVariable(origin_name),
                Kratos.KratosGlobals.GetVariable(destination_name))


def GatherGridArray(model_part: Kratos.ModelPart, variable_names, grid_shape) -> numpy.ndarray:
    """Reads mapped nodal fields off a background-grid part as a grid array.

    Nodes are ordered by their coordinates (x outermost, z innermost), so no
    assumption on the mesh generator's numbering is needed - only that the
    nodes form a full (D, H, W) lattice.

    Args:
        model_part: A background-grid model part.
        variable_names: historical nodal variable names; each contributes
            its flattened per-node width to the channel axis.
        grid_shape: (D, H, W) lattice points per axis (divisions + 1); 2D
            grids use (H, W, 1)-style shapes with a size-1 trailing axis.

    Returns:
        (C, *grid_shape) float64 array with axes (channel, x, y, z).
    """
    grid_shape = tuple(int(n) for n in grid_shape)
    expected = int(numpy.prod(grid_shape))
    if model_part.NumberOfNodes() != expected:
        raise ValueError(
            f"Model part \"{model_part.FullName()}\" has {model_part.NumberOfNodes()} nodes "
            f"but grid_shape {grid_shape} implies {expected}.")

    coordinates = numpy.array([[node.X, node.Y, node.Z] for node in model_part.Nodes])
    # snap to integer lattice indices per axis (robust against the floating
    # noise of generated coordinates), then order lexicographically
    indices = numpy.zeros_like(coordinates, dtype=numpy.int64)
    for axis in range(3):
        low, high = coordinates[:, axis].min(), coordinates[:, axis].max()
        if grid_shape[axis] > 1 and high > low:
            spacing = (high - low) / (grid_shape[axis] - 1)
            indices[:, axis] = numpy.rint((coordinates[:, axis] - low) / spacing)
    order = numpy.lexsort((indices[:, 2], indices[:, 1], indices[:, 0]))

    channels = []
    for name in variable_names:
        variable = Kratos.KratosGlobals.GetVariable(name)
        values = numpy.array(
            [numpy.atleast_1d(numpy.asarray(node.GetSolutionStepValue(variable)))
             for node in model_part.Nodes])[order]  # (N, width), lattice order
        for component in range(values.shape[1]):
            channels.append(values[:, component].reshape(grid_shape))
    return numpy.stack(channels)
