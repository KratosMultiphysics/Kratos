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

    # 1. Build points array
    points = np.zeros((numNodes, 3), dtype=np.float64)
    nodeIdToIdx = {}
    for idx, node in enumerate(modelPart.Nodes):
        nodeIdToIdx[node.Id] = idx
        if useDeformedConfiguration:
            points[idx, 0] = node.X
            points[idx, 1] = node.Y
            points[idx, 2] = node.Z
        else:
            points[idx, 0] = node.X0
            points[idx, 1] = node.Y0
            points[idx, 2] = node.Z0

    # 2. Build cells and cell types
    cells = []
    cellTypes = []

    for entity in entities:
        geom = entity.GetGeometry()
        geomType = geom.GetGeometryType()
        vtkType = GEOMETRY_TYPE_TO_VTK_CELL_TYPE.get(geomType)

        if vtkType is None:
            continue

        cellNodes = [nodeIdToIdx[node.Id] for node in geom]
        cells.append(len(cellNodes))
        cells.extend(cellNodes)
        cellTypes.append(vtkType)

    if len(cells) == 0:
        return pv.UnstructuredGrid()

    grid = pv.UnstructuredGrid(np.array(cells, dtype=np.int32), np.array(cellTypes, dtype=np.uint8), points)

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
