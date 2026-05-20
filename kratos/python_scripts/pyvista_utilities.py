# Import Kratos Multiphysics
import KratosMultiphysics as KM

# Import PyVista, along with NumPy for array handling
import numpy as np
import pyvista as pv

# Mapping of Kratos Geometry types to VTK Cell Types (using their integer values)
# NOTE: Ref Kratos VtkDefinitions mapping in vtk_definition.cpp
GEOMETRY_TYPE_TO_VTK_CELL_TYPE = {
    KM.GeometryData.KratosGeometryType.Kratos_Point2D: 1,
    KM.GeometryData.KratosGeometryType.Kratos_Point3D: 1,
    KM.GeometryData.KratosGeometryType.Kratos_Line2D2: 3,
    KM.GeometryData.KratosGeometryType.Kratos_Line3D2: 3,
    KM.GeometryData.KratosGeometryType.Kratos_Triangle2D3: 5,
    KM.GeometryData.KratosGeometryType.Kratos_Triangle3D3: 5,
    KM.GeometryData.KratosGeometryType.Kratos_Quadrilateral2D4: 9,
    KM.GeometryData.KratosGeometryType.Kratos_Quadrilateral3D4: 9,
    KM.GeometryData.KratosGeometryType.Kratos_Tetrahedra3D4: 10,
    KM.GeometryData.KratosGeometryType.Kratos_Hexahedra3D8: 12,
    KM.GeometryData.KratosGeometryType.Kratos_Prism3D6: 13,
    KM.GeometryData.KratosGeometryType.Kratos_Pyramid3D5: 14,
    KM.GeometryData.KratosGeometryType.Kratos_Line2D3: 21,
    KM.GeometryData.KratosGeometryType.Kratos_Line3D3: 21,
    KM.GeometryData.KratosGeometryType.Kratos_Triangle2D6: 22,
    KM.GeometryData.KratosGeometryType.Kratos_Triangle3D6: 22,
    KM.GeometryData.KratosGeometryType.Kratos_Quadrilateral2D8: 23,
    KM.GeometryData.KratosGeometryType.Kratos_Quadrilateral3D8: 23,
    KM.GeometryData.KratosGeometryType.Kratos_Tetrahedra3D10: 24,
    KM.GeometryData.KratosGeometryType.Kratos_Hexahedra3D20: 25,
    KM.GeometryData.KratosGeometryType.Kratos_Prism3D15: 26,
    KM.GeometryData.KratosGeometryType.Kratos_Pyramid3D13: 27,
    KM.GeometryData.KratosGeometryType.Kratos_Quadrilateral2D9: 28,
    KM.GeometryData.KratosGeometryType.Kratos_Quadrilateral3D9: 28,
    KM.GeometryData.KratosGeometryType.Kratos_Hexahedra3D27: 29
}


def _CreateUnstructuredGrid(modelPart, entities, useDeformedConfiguration, nodalVariables, entityVariables):
    """Helper function to convert a collection of entities (Elements or Conditions) into a pyvista.UnstructuredGrid."""
    if len(entities) == 0:
        return pv.UnstructuredGrid()

    numNodes = len(modelPart.Nodes)
    if numNodes == 0:
        return pv.UnstructuredGrid()

    # 1. Build points array, id→index dict, and flat node-ID array for vectorised remapping
    points = np.zeros((numNodes, 3), dtype=np.float64)
    nodeIdToIdx = {}
    nodeIdArr = np.empty(numNodes, dtype=np.int64)
    for idx, node in enumerate(modelPart.Nodes):
        nodeIdArr[idx] = node.Id
        nodeIdToIdx[node.Id] = idx
        if useDeformedConfiguration:
            points[idx, 0] = node.X
            points[idx, 1] = node.Y
            points[idx, 2] = node.Z
        else:
            points[idx, 0] = node.X0
            points[idx, 1] = node.Y0
            points[idx, 2] = node.Z0

    # idToIdx[kratos_node_id] = 0-based point index, for O(1) vectorised lookup
    maxNodeId = int(nodeIdArr.max())
    idToIdx = np.empty(maxNodeId + 1, dtype=np.int32)
    idToIdx[nodeIdArr] = np.arange(numNodes, dtype=np.int32)

    # 2. Build cells and cell types
    # Fast path: ConnectivityIdsTensorAdaptor fills a (N, K) numpy array of node IDs
    # in parallel C++. Check() verifies uniform geometry; throws for mixed meshes.
    cells = None
    cellTypes = None

    firstGeomType = next(iter(entities)).GetGeometry().GetGeometryType()
    firstVtkType = GEOMETRY_TYPE_TO_VTK_CELL_TYPE.get(firstGeomType)
    if firstVtkType is not None:
        try:
            adaptor = KM.TensorAdaptors.ConnectivityIdsTensorAdaptor(entities)
            adaptor.Check()
            adaptor.CollectData()
            nodeIds = adaptor.data  # shape (N, K): Kratos node IDs (1-based, int32)
            N, K = nodeIds.shape
            nodeIndices = idToIdx[nodeIds]  # vectorised remap to 0-based point indices
            countCol = np.full((N, 1), K, dtype=np.int32)
            cells = np.hstack([countCol, nodeIndices]).ravel().astype(np.int32)
            cellTypes = np.full(N, firstVtkType, dtype=np.uint8)
        except Exception:
            pass  # mixed geometry or unsupported type; fall through to slow path

    if cells is None:
        # Slow path: mixed or unsupported geometry types
        cellsList = []
        cellTypesList = []
        for entity in entities:
            geom = entity.GetGeometry()
            geomType = geom.GetGeometryType()
            vtkType = GEOMETRY_TYPE_TO_VTK_CELL_TYPE.get(geomType)
            if vtkType is None:
                continue
            cellNodes = [nodeIdToIdx[node.Id] for node in geom]
            cellsList.append(len(cellNodes))
            cellsList.extend(cellNodes)
            cellTypesList.append(vtkType)
        if len(cellsList) == 0:
            return pv.UnstructuredGrid()
        cells = np.array(cellsList, dtype=np.int32)
        cellTypes = np.array(cellTypesList, dtype=np.uint8)

    if len(cellTypes) == 0:
        return pv.UnstructuredGrid()

    grid = pv.UnstructuredGrid(cells, cellTypes, points)

    # 3. Add point (nodal) variables
    for var in nodalVariables:
        if isinstance(var, str):
            var = KM.KratosGlobals.GetVariable(var)

        firstNode = next(iter(modelPart.Nodes))
        isHistorical = firstNode.HasSolutionStepValue(var)
        isNonHistorical = firstNode.Has(var)

        if not (isHistorical or isNonHistorical):
            continue

        val = firstNode.GetSolutionStepValue(var) if isHistorical else firstNode.GetValue(var)

        if isinstance(val, (float, int)):
            arr = np.zeros(numNodes, dtype=np.float64)
            for idx, node in enumerate(modelPart.Nodes):
                arr[idx] = node.GetSolutionStepValue(var) if isHistorical else node.GetValue(var)
            grid.point_data[var.Name()] = arr
        elif isinstance(val, (KM.Array3, KM.Vector, list)) or type(val).__name__ == "Vector":
            arr = np.zeros((numNodes, 3), dtype=np.float64)
            for idx, node in enumerate(modelPart.Nodes):
                v = node.GetSolutionStepValue(var) if isHistorical else node.GetValue(var)
                arr[idx, 0] = v[0]
                arr[idx, 1] = v[1]
                arr[idx, 2] = v[2]
            grid.point_data[var.Name()] = arr

    # 4. Add cell (entity) variables
    for var in entityVariables:
        if isinstance(var, str):
            var = KM.KratosGlobals.GetVariable(var)

        valType = None
        for entity in entities:
            geomType = entity.GetGeometry().GetGeometryType()
            if geomType in GEOMETRY_TYPE_TO_VTK_CELL_TYPE and entity.Has(var):
                valType = type(entity.GetValue(var))
                break

        if valType is None:
            continue

        if issubclass(valType, (float, int)):
            arr = np.zeros(len(cellTypes), dtype=np.float64)
            cellIdx = 0
            for entity in entities:
                geomType = entity.GetGeometry().GetGeometryType()
                if geomType not in GEOMETRY_TYPE_TO_VTK_CELL_TYPE:
                    continue
                if entity.Has(var):
                    arr[cellIdx] = entity.GetValue(var)
                cellIdx += 1
            grid.cell_data[var.Name()] = arr
        elif issubclass(valType, (KM.Array3, KM.Vector, list)) or valType.__name__ == "Vector":
            arr = np.zeros((len(cellTypes), 3), dtype=np.float64)
            cellIdx = 0
            for entity in entities:
                geomType = entity.GetGeometry().GetGeometryType()
                if geomType not in GEOMETRY_TYPE_TO_VTK_CELL_TYPE:
                    continue
                if entity.Has(var):
                    v = entity.GetValue(var)
                    arr[cellIdx, 0] = v[0]
                    arr[cellIdx, 1] = v[1]
                    arr[cellIdx, 2] = v[2]
                cellIdx += 1
            grid.cell_data[var.Name()] = arr

    return grid


def ModelPartToPyVista(modelPart, useDeformedConfiguration=False, nodalVariables=[], elementVariables=[], conditionVariables=[], exportElements=True, exportConditions=False):
    """Converts a Kratos ModelPart into a pyvista.UnstructuredGrid or pyvista.MultiBlock mesh.

    Parameters:
    modelPart (KM.ModelPart): The Kratos ModelPart to convert.
    useDeformedConfiguration (bool): If True, use current coordinates (X, Y, Z). If False, use initial coordinates (X0, Y0, Z0).
    nodalVariables (list): List of KM.Variable or string names for nodal data to export.
    elementVariables (list): List of KM.Variable or string names for element (cell) data to export.
    conditionVariables (list): List of KM.Variable or string names for condition (cell) data to export.
    exportElements (bool): Whether to export elements.
    exportConditions (bool): Whether to export conditions.

    Returns:
    pyvista.UnstructuredGrid or pyvista.MultiBlock: The converted PyVista mesh.
    """
    if exportElements and exportConditions:
        elemGrid = _CreateUnstructuredGrid(
            modelPart,
            modelPart.Elements,
            useDeformedConfiguration,
            nodalVariables,
            elementVariables
        )
        condGrid = _CreateUnstructuredGrid(
            modelPart,
            modelPart.Conditions,
            useDeformedConfiguration,
            nodalVariables,
            conditionVariables
        )
        blocks = pv.MultiBlock()
        blocks["elements"] = elemGrid
        blocks["conditions"] = condGrid
        return blocks
    elif exportElements:
        return _CreateUnstructuredGrid(
            modelPart,
            modelPart.Elements,
            useDeformedConfiguration,
            nodalVariables,
            elementVariables
        )
    elif exportConditions:
        return _CreateUnstructuredGrid(
            modelPart,
            modelPart.Conditions,
            useDeformedConfiguration,
            nodalVariables,
            conditionVariables
        )
    else:
        return pv.UnstructuredGrid()


def PlotModelPart(modelPart, scalars=None, useDeformedConfiguration=False, **kwargs):
    """Plot a Kratos ModelPart using PyVista's interactive plotting engine."""
    nodalVars = [scalars] if scalars is not None else []
    grid = ModelPartToPyVista(
        modelPart,
        useDeformedConfiguration=useDeformedConfiguration,
        nodalVariables=nodalVars,
        exportElements=len(modelPart.Elements) > 0,
        exportConditions=len(modelPart.Conditions) > 0 and len(modelPart.Elements) == 0
    )

    if scalars is not None:
        scalarsName = scalars if isinstance(scalars, str) else scalars.Name()
        if isinstance(grid, pv.MultiBlock):
            for block in grid:
                if scalarsName in block.point_data:
                    block.set_active_scalars(scalarsName)
        else:
            if scalarsName in grid.point_data:
                grid.set_active_scalars(scalarsName)

    return grid.plot(**kwargs)


def SaveModelPart(modelPart, filename, useDeformedConfiguration=False, nodalVariables=[], elementVariables=[], conditionVariables=[]):
    """Save a Kratos ModelPart to a VTK/VTU/VTM file using PyVista's save mechanisms."""
    exportElements = len(modelPart.Elements) > 0
    exportConditions = len(modelPart.Conditions) > 0

    grid = ModelPartToPyVista(
        modelPart,
        useDeformedConfiguration=useDeformedConfiguration,
        nodalVariables=nodalVariables,
        elementVariables=elementVariables,
        conditionVariables=conditionVariables,
        exportElements=exportElements,
        exportConditions=exportConditions
    )
    grid.save(filename)


# ==============================================================================
# ADVANCED PYVISTA POST-PROCESSING UTILITIES (Kratos Convention)
# ==============================================================================

def ComputeWarpedMesh(modelPart, vectorVariable, factor=1.0, nodalVariables=[], elementVariables=[], conditionVariables=[]):
    """Generates a deformed mesh representation by warping the points based on a vector variable (e.g. DISPLACEMENT)."""
    grid = ModelPartToPyVista(
        modelPart,
        useDeformedConfiguration=False,  # We deforms points manually with the vector field
        nodalVariables=nodalVariables + [vectorVariable],
        elementVariables=elementVariables,
        conditionVariables=conditionVariables,
        exportElements=len(modelPart.Elements) > 0,
        exportConditions=len(modelPart.Conditions) > 0
    )
    varName = vectorVariable if isinstance(vectorVariable, str) else vectorVariable.Name()

    if isinstance(grid, pv.MultiBlock):
        for key in grid.keys():
            if grid[key].n_points > 0 and varName in grid[key].point_data:
                grid[key] = grid[key].warp_by_vector(varName, factor=factor)
        return grid
    else:
        if grid.n_points > 0 and varName in grid.point_data:
            return grid.warp_by_vector(varName, factor=factor)
        return grid


def CreateOrthogonalSlices(modelPart, x=None, y=None, z=None, nodalVariables=[], elementVariables=[]):
    """Extracts orthogonal slices from a 3D Kratos ModelPart domain using interactive cutting planes."""
    grid = ModelPartToPyVista(
        modelPart,
        useDeformedConfiguration=False,
        nodalVariables=nodalVariables,
        elementVariables=elementVariables,
        exportElements=True,
        exportConditions=False
    )
    if isinstance(grid, pv.MultiBlock):
        grid = grid["elements"]

    if grid.n_points == 0:
        return grid

    bounds = grid.bounds
    if x is None:
        x = 0.5 * (bounds[0] + bounds[1])
    if y is None:
        y = 0.5 * (bounds[2] + bounds[3])
    if z is None:
        z = 0.5 * (bounds[4] + bounds[5])

    return grid.slice_orthogonal(x=x, y=y, z=z)


def CreateIsosurfaces(modelPart, variable, valuesOrNumber=5, nodalVariables=[]):
    """Generates contour surfaces (isosurfaces) from a Kratos ModelPart based on a nodal scalar variable."""
    varName = variable if isinstance(variable, str) else variable.Name()
    grid = ModelPartToPyVista(
        modelPart,
        useDeformedConfiguration=False,
        nodalVariables=nodalVariables + [variable],
        exportElements=True,
        exportConditions=False
    )
    if isinstance(grid, pv.MultiBlock):
        grid = grid["elements"]

    if grid.n_points == 0:
        return grid

    return grid.contour(isosurfaces=valuesOrNumber, scalars=varName)


def CreateStreamlines(modelPart, velocityVariable="VELOCITY", sourceCenter=None, sourceRadius=None, nPoints=100, integrationDirection="both"):
    """Generates fluid streamlines through a vector velocity field starting from a spherical seed region."""
    varName = velocityVariable if isinstance(velocityVariable, str) else velocityVariable.Name()
    grid = ModelPartToPyVista(
        modelPart,
        useDeformedConfiguration=False,
        nodalVariables=[velocityVariable],
        exportElements=True,
        exportConditions=False
    )
    if isinstance(grid, pv.MultiBlock):
        grid = grid["elements"]

    if grid.n_points == 0:
        return grid

    bounds = grid.bounds
    if sourceCenter is None:
        sourceCenter = [
            0.5 * (bounds[0] + bounds[1]),
            0.5 * (bounds[2] + bounds[3]),
            0.5 * (bounds[4] + bounds[5])
        ]
    if sourceRadius is None:
        diag = np.sqrt((bounds[1]-bounds[0])**2 + (bounds[3]-bounds[2])**2 + (bounds[5]-bounds[4])**2)
        sourceRadius = 0.1 * diag

    return grid.streamlines(
        vectors=varName,
        source_center=sourceCenter,
        source_radius=sourceRadius,
        n_points=nPoints,
        integration_direction=integrationDirection
    )


def CreateThresholdedMesh(modelPart, variable, thresholdValue, thresholdType="above", nodalVariables=[], elementVariables=[]):
    """Extracts a subset of the mesh where a scalar variable satisfies a threshold condition.

    Parameters:
    modelPart (KM.ModelPart): The Kratos ModelPart to convert and filter.
    variable (KM.Variable or str): Nodal or cell scalar variable to threshold by.
    thresholdValue (float or tuple): The threshold value (float) or range (tuple/list of two floats).
    thresholdType (str): Type of threshold: 'above', 'below', or 'between'.
    nodalVariables (list): Additional nodal variables to export.
    elementVariables (list): Additional elemental variables to export.

    Returns:
    pyvista.UnstructuredGrid: The thresholded mesh subset.
    """
    varName = variable if isinstance(variable, str) else variable.Name()

    grid = ModelPartToPyVista(
        modelPart,
        useDeformedConfiguration=False,
        nodalVariables=nodalVariables + [variable],
        elementVariables=elementVariables + [variable],
        exportElements=True,
        exportConditions=False
    )
    if isinstance(grid, pv.MultiBlock):
        grid = grid["elements"]

    if grid.n_points == 0:
        return grid

    # Determine if variable is in point_data or cell_data
    if varName in grid.point_data:
        preference = 'point'
    elif varName in grid.cell_data:
        preference = 'cell'
    else:
        # Fallback to point
        preference = 'point'

    if thresholdType == "above":
        return grid.threshold(value=thresholdValue, scalars=varName, method='upper', preference=preference)
    elif thresholdType == "below":
        return grid.threshold(value=thresholdValue, scalars=varName, method='lower', preference=preference)
    elif thresholdType == "between":
        if not isinstance(thresholdValue, (list, tuple)) or len(thresholdValue) != 2:
            raise ValueError("For 'between' thresholdType, thresholdValue must be a tuple/list of (lower_limit, upper_limit).")
        return grid.threshold(value=thresholdValue, scalars=varName, preference=preference)
    else:
        raise ValueError(f"Unknown thresholdType '{thresholdType}'. Supported types are 'above', 'below', 'between'.")


def CreateVectorGlyphs(modelPart, vectorVariable, scaleFactor=1.0, glyphType="arrow", nodalVariables=[]):
    """Generates 3D glyphs (e.g. arrows) at nodes of a Kratos ModelPart oriented and scaled by a vector variable.

    Parameters:
    modelPart (KM.ModelPart): The Kratos ModelPart to convert and generate glyphs for.
    vectorVariable (KM.Variable or str): Nodal vector variable to orient and scale glyphs by.
    scaleFactor (float): Scaling factor for the glyphs.
    glyphType (str): Type of glyph: 'arrow', 'cone', 'sphere', 'cylinder'.
    nodalVariables (list): Additional nodal variables to export.

    Returns:
    pyvista.PolyData: The generated glyphs.
    """
    varName = vectorVariable if isinstance(vectorVariable, str) else vectorVariable.Name()

    grid = ModelPartToPyVista(
        modelPart,
        useDeformedConfiguration=False,
        nodalVariables=nodalVariables + [vectorVariable],
        exportElements=True,
        exportConditions=False
    )
    if isinstance(grid, pv.MultiBlock):
        grid = grid["elements"]

    if grid.n_points == 0:
        return pv.PolyData()

    if varName not in grid.point_data:
        raise ValueError(f"Vector variable '{varName}' not found in the converted mesh point data.")

    # Determine glyph geometry
    if glyphType == "arrow":
        geom = pv.Arrow()
    elif glyphType == "cone":
        geom = pv.Cone()
    elif glyphType == "sphere":
        geom = pv.Sphere()
    elif glyphType == "cylinder":
        geom = pv.Cylinder()
    else:
        geom = pv.Arrow()

    glyphs = grid.glyph(orient=varName, scale=varName, factor=scaleFactor, geom=geom)
    if 'GlyphVector' in glyphs.point_data:
        glyphs.point_data[varName] = glyphs.point_data['GlyphVector']
    if 'GlyphScale' in glyphs.point_data:
        glyphs.point_data[varName + "_magnitude"] = glyphs.point_data['GlyphScale']
    return glyphs


def CreateClippedMesh(modelPart, normal=[0.0, 0.0, 1.0], origin=None, nodalVariables=[], elementVariables=[]):
    """Clips a Kratos ModelPart domain by a plane defined by an origin and normal vector.

    Parameters:
    modelPart (KM.ModelPart): The Kratos ModelPart to convert and clip.
    normal (tuple or list): 3D normal vector of the clipping plane.
    origin (tuple or list): 3D point defining the origin of the clipping plane. If None, uses the center of the mesh bounds.
    nodalVariables (list): Nodal variables to export.
    elementVariables (list): Elemental variables to export.

    Returns:
    pyvista.UnstructuredGrid: The clipped mesh subset.
    """
    grid = ModelPartToPyVista(
        modelPart,
        useDeformedConfiguration=False,
        nodalVariables=nodalVariables,
        elementVariables=elementVariables,
        exportElements=True,
        exportConditions=False
    )
    if isinstance(grid, pv.MultiBlock):
        grid = grid["elements"]

    if grid.n_points == 0:
        return grid

    if origin is None:
        bounds = grid.bounds
        origin = [
            0.5 * (bounds[0] + bounds[1]),
            0.5 * (bounds[2] + bounds[3]),
            0.5 * (bounds[4] + bounds[5])
        ]

    return grid.clip(normal=normal, origin=origin)

