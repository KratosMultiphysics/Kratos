# Import Kratos Multiphysics
import KratosMultiphysics as KM

# Import NumPy for array handling
import numpy as np

# NOTE: pyvista is an optional dependency. Every function that needs it imports
# it locally so the rest of Kratos can load successfully when pyvista is absent.

# Mapping of Kratos Geometry types to VTK Cell Types (integer values).
# Ref: Kratos VtkDefinitions mapping in vtk_definition.cpp
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

# Component suffix tables used by _BuildScalarLabel.
_VECTOR_SUFFIXES  = ["X", "Y", "Z"]
_VOIGT_SUFFIXES   = ["XX", "YY", "ZZ", "XY", "YZ", "XZ"]
_TENSOR9_SUFFIXES = ["XX", "XY", "XZ", "YX", "YY", "YZ", "ZX", "ZY", "ZZ"]
_PRINCIPAL_LABELS = ["(Max. Principal)", "(Int. Principal)", "(Min. Principal)"]


# ==============================================================================
# Internal helpers
# ==============================================================================

def _GetVarName(variable):
    """Return the variable name string, accepting either a KM.Variable or a plain string."""
    return variable if isinstance(variable, str) else variable.Name()


def _BuildScalarLabel(name, component, ndim):
    """Build a scalar bar title from array *name*, *component* index, and array *ndim*.

    Handles vector (3), Voigt (6), full-tensor (9), and principal-value arrays
    using the same label conventions as FElupe's Scene.plot().

    Parameters:
    name (str): Variable name (may contain 'Principal Values of ' prefix).
    component (int or None): Component index; None means magnitude.
    ndim (int): Number of components in the original data array.

    Returns:
    str: Human-readable scalar bar label.
    """
    is_principal = "Principal Values of " in name
    display_name = name.replace("Principal Values of ", "") if is_principal else name

    if component is None:
        suffix = "(Magnitude)" if is_principal else "Magnitude"
    elif ndim <= 1:
        suffix = ""
    elif is_principal:
        suffix = _PRINCIPAL_LABELS[min(component, 2)]
    elif ndim <= 3:
        suffix = _VECTOR_SUFFIXES[min(component, 2)]
    elif ndim == 6:
        suffix = _VOIGT_SUFFIXES[min(component, 5)]
    elif ndim == 9:
        suffix = _TENSOR9_SUFFIXES[min(component, 8)]
    else:
        suffix = str(component)

    return f"{display_name} {suffix}".strip()


def _ExtractComponent(arr, component, name):
    """Extract a 1-D scalar array from a possibly multi-component NumPy array.

    For arrays whose *name* contains 'Principal Values of', the columns are
    sorted in descending order (largest eigenvalue first) before extraction,
    matching FElupe's convention.

    Parameters:
    arr (np.ndarray): Shape (N,) for scalars or (N, K) for multi-component.
    component (int or None): Component to extract; None returns the L2 norm.
    name (str): Variable name used for principal-value detection.

    Returns:
    (np.ndarray, int): (1-D scalar array, original ndim).
    """
    if arr.ndim == 1:
        return arr.copy(), 1

    ndim = arr.shape[-1]
    work = arr.copy()

    if "Principal Values of " in name:
        work = np.flip(np.sort(work, axis=-1), axis=-1)

    if component is None:
        return np.linalg.norm(work, axis=-1), ndim

    c = min(int(component), ndim - 1)
    return work[:, c].copy(), ndim


# ==============================================================================
# Core conversion
# ==============================================================================

def _CreateUnstructuredGrid(modelPart, entities, useDeformedConfiguration, nodalVariables, entityVariables):
    """Convert a collection of Kratos entities (Elements or Conditions) to a pyvista.UnstructuredGrid.

    Handles scalar, 3-D vector (Array3 / Vector), and Matrix nodal and entity
    variables. Matrix values are stored as row-major flat arrays of shape
    (N, rows * cols).
    """
    import pyvista as pv

    if len(entities) == 0:
        return pv.UnstructuredGrid()

    numNodes = len(modelPart.Nodes)
    if numNodes == 0:
        return pv.UnstructuredGrid()

    # 1. Build points array and node-ID → index mappings
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

    maxNodeId = int(nodeIdArr.max())
    idToIdx = np.empty(maxNodeId + 1, dtype=np.int32)
    idToIdx[nodeIdArr] = np.arange(numNodes, dtype=np.int32)

    # 2. Build cell connectivity and cell-type arrays
    cells = None
    cellTypes = None

    firstGeomType = next(iter(entities)).GetGeometry().GetGeometryType()
    firstVtkType = GEOMETRY_TYPE_TO_VTK_CELL_TYPE.get(firstGeomType)
    if firstVtkType is not None:
        try:
            # Fast path: uniform geometry — fill via vectorised C++ adaptor
            adaptor = KM.TensorAdaptors.ConnectivityIdsTensorAdaptor(entities)
            adaptor.Check()
            adaptor.CollectData()
            nodeIds = adaptor.data  # (N, K) int32 Kratos node IDs (1-based)
            N, K = nodeIds.shape
            nodeIndices = idToIdx[nodeIds]
            countCol = np.full((N, 1), K, dtype=np.int32)
            cells = np.hstack([countCol, nodeIndices]).ravel().astype(np.int32)
            cellTypes = np.full(N, firstVtkType, dtype=np.uint8)
        except Exception:
            pass  # mixed geometry or unsupported type — fall through to slow path

    if cells is None:
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

    # 3. Add nodal (point) variables
    for var in nodalVariables:
        if isinstance(var, str):
            var = KM.KratosGlobals.GetVariable(var)

        firstNode = next(iter(modelPart.Nodes))
        isHistorical   = firstNode.HasSolutionStepValue(var)
        isNonHistorical = firstNode.Has(var)

        if not (isHistorical or isNonHistorical):
            continue

        val = firstNode.GetSolutionStepValue(var) if isHistorical else firstNode.GetValue(var)
        valTypeName = type(val).__name__

        def _get(node):
            return node.GetSolutionStepValue(var) if isHistorical else node.GetValue(var)

        if isinstance(val, (float, int, bool)):
            arr = np.array([_get(n) for n in modelPart.Nodes], dtype=np.float64)
            grid.point_data[var.Name()] = arr

        elif isinstance(val, KM.Array3) or valTypeName in ("Array3", "Vector"):
            try:
                sz = len(val)
            except Exception:
                sz = 3
            ncols = max(sz, 3)
            arr = np.zeros((numNodes, ncols), dtype=np.float64)
            for idx, node in enumerate(modelPart.Nodes):
                v = _get(node)
                for k in range(min(len(v), ncols)):
                    arr[idx, k] = v[k]
            if ncols < 3:
                arr = np.pad(arr, ((0, 0), (0, 3 - ncols)))
            grid.point_data[var.Name()] = arr

        elif valTypeName == "Matrix":
            r, c = val.Size1(), val.Size2()
            arr = np.zeros((numNodes, r * c), dtype=np.float64)
            for idx, node in enumerate(modelPart.Nodes):
                m = _get(node)
                for i in range(r):
                    for j in range(c):
                        arr[idx, i * c + j] = m[i, j]
            grid.point_data[var.Name()] = arr

    # 4. Add entity (cell) variables
    nCells = len(cellTypes)
    for var in entityVariables:
        if isinstance(var, str):
            var = KM.KratosGlobals.GetVariable(var)

        sampleVal = None
        for entity in entities:
            geomType = entity.GetGeometry().GetGeometryType()
            if geomType in GEOMETRY_TYPE_TO_VTK_CELL_TYPE and entity.Has(var):
                sampleVal = entity.GetValue(var)
                break

        if sampleVal is None:
            continue

        valType     = type(sampleVal)
        valTypeName = valType.__name__

        if isinstance(sampleVal, (float, int, bool)):
            arr = np.zeros(nCells, dtype=np.float64)
            cellIdx = 0
            for entity in entities:
                if entity.GetGeometry().GetGeometryType() not in GEOMETRY_TYPE_TO_VTK_CELL_TYPE:
                    continue
                if entity.Has(var):
                    arr[cellIdx] = entity.GetValue(var)
                cellIdx += 1
            grid.cell_data[var.Name()] = arr

        elif isinstance(sampleVal, KM.Array3) or valTypeName in ("Array3", "Vector"):
            try:
                sz = len(sampleVal)
            except Exception:
                sz = 3
            ncols = max(sz, 3)
            arr = np.zeros((nCells, ncols), dtype=np.float64)
            cellIdx = 0
            for entity in entities:
                if entity.GetGeometry().GetGeometryType() not in GEOMETRY_TYPE_TO_VTK_CELL_TYPE:
                    continue
                if entity.Has(var):
                    v = entity.GetValue(var)
                    for k in range(min(len(v), ncols)):
                        arr[cellIdx, k] = v[k]
                cellIdx += 1
            if ncols < 3:
                arr = np.pad(arr, ((0, 0), (0, 3 - ncols)))
            grid.cell_data[var.Name()] = arr

        elif valTypeName == "Matrix":
            r, c = sampleVal.Size1(), sampleVal.Size2()
            arr = np.zeros((nCells, r * c), dtype=np.float64)
            cellIdx = 0
            for entity in entities:
                if entity.GetGeometry().GetGeometryType() not in GEOMETRY_TYPE_TO_VTK_CELL_TYPE:
                    continue
                if entity.Has(var):
                    m = entity.GetValue(var)
                    for i in range(r):
                        for j in range(c):
                            arr[cellIdx, i * c + j] = m[i, j]
                cellIdx += 1
            grid.cell_data[var.Name()] = arr

    return grid


def ModelPartToPyVista(
    modelPart,
    useDeformedConfiguration=False,
    nodalVariables=[],
    elementVariables=[],
    conditionVariables=[],
    exportElements=True,
    exportConditions=False
):
    """Convert a Kratos ModelPart to a pyvista.UnstructuredGrid or pyvista.MultiBlock mesh.

    Parameters:
    modelPart (KM.ModelPart): The Kratos ModelPart to convert.
    useDeformedConfiguration (bool): Use current (X/Y/Z) vs initial (X0/Y0/Z0) coordinates.
    nodalVariables (list): KM.Variable or str names for nodal data to export.
    elementVariables (list): KM.Variable or str names for element cell data.
    conditionVariables (list): KM.Variable or str names for condition cell data.
    exportElements (bool): Export elements.
    exportConditions (bool): Export conditions.

    Returns:
    pyvista.UnstructuredGrid or pyvista.MultiBlock: Converted mesh.
    """
    import pyvista as pv

    if exportElements and exportConditions:
        elemGrid = _CreateUnstructuredGrid(
            modelPart, modelPart.Elements,
            useDeformedConfiguration, nodalVariables, elementVariables
        )
        condGrid = _CreateUnstructuredGrid(
            modelPart, modelPart.Conditions,
            useDeformedConfiguration, nodalVariables, conditionVariables
        )
        blocks = pv.MultiBlock()
        blocks["elements"]   = elemGrid
        blocks["conditions"] = condGrid
        return blocks
    elif exportElements:
        return _CreateUnstructuredGrid(
            modelPart, modelPart.Elements,
            useDeformedConfiguration, nodalVariables, elementVariables
        )
    elif exportConditions:
        return _CreateUnstructuredGrid(
            modelPart, modelPart.Conditions,
            useDeformedConfiguration, nodalVariables, conditionVariables
        )
    else:
        return pv.UnstructuredGrid()


# ==============================================================================
# Plotting
# ==============================================================================

def PlotModelPart(
    modelPart,
    name=None,
    component=0,
    label=None,
    factor=1.0,
    warpByVector=None,
    showEdges=True,
    showUndeformed=True,
    cmap="turbo",
    view="default",
    theme=None,
    scalarBarArgs=None,
    scalarBarVertical=False,
    addAxes=True,
    offScreen=False,
    plotter=None,
    notebook=False,
    extractSurface=False,
    nonlinearSubdivision=1,
    smoothShading=True,
    splitSharpEdges=True,
    edgeColor="black",
    lineWidth=1.0,
    useDeformedConfiguration=False,
    nodalVariables=None,
    elementVariables=None,
    conditionVariables=None,
    exportElements=True,
    exportConditions=False,
    **kwargs
):
    """Plot a Kratos ModelPart using PyVista's rendering engine.

    The function converts the ModelPart, optionally warps it by a vector field,
    selects a scalar for coloring (with component extraction and principal-value
    sorting), and returns a configured pyvista.Plotter.  The caller decides
    whether to call .show() interactively or .screenshot() for off-screen images.

    Parameters:
    modelPart (KM.ModelPart): ModelPart to render.
    name (KM.Variable or str or None): Variable to color by.
    component (int or None): Component index; None shows the L2 magnitude.
        Suffix labels: vectors→X/Y/Z, Voigt-6→XX/YY/ZZ/XY/YZ/XZ,
        tensor-9→XX/XY/.../ZZ, principal-3→(Max./Int./Min. Principal).
    label (str or None): Override scalar bar title (auto-generated when None).
    factor (float): Warp scaling factor applied when warpByVector is set.
    warpByVector (KM.Variable or str or None): Nodal vector variable used to
        warp the mesh (e.g. KM.DISPLACEMENT). When set the undeformed mesh is
        shown as a 20%-opacity ghost (controlled by showUndeformed).
    showEdges (bool): Show cell edges.
    showUndeformed (bool): Show undeformed ghost mesh when warping is active.
    cmap (str): Matplotlib colormap name.
    view (str or None): Camera preset: 'xy', 'xz', 'yz', 'iso', or 'default'.
        'default' auto-selects 'xy' + parallel projection for flat 2-D meshes.
    theme (str or None): PyVista theme ('default', 'document', 'dark', …).
    scalarBarArgs (dict or None): Extra kwargs forwarded to the scalar bar.
    scalarBarVertical (bool): Show scalar bar vertically on the right side.
    addAxes (bool): Add coordinate axes widget.
    offScreen (bool): Initialize plotter in off-screen mode (for screenshots).
    plotter (pyvista.Plotter or None): Reuse an existing plotter to compose scenes.
    notebook (bool): Enable Jupyter inline rendering (implies off-screen).
    extractSurface (bool): Extract surface before rendering (hides internal edges
        of quadratic cells).
    nonlinearSubdivision (int): >1 enables smooth curved surface of quadratic cells
        (implies extractSurface).
    smoothShading (bool): Smooth shading when nonlinearSubdivision > 1.
    splitSharpEdges (bool): Split sharp edges for smooth shading.
    edgeColor (str): Edge line colour.
    lineWidth (float): Edge line width.
    useDeformedConfiguration (bool): Use current (X/Y/Z) coordinates.
    nodalVariables (list or None): Additional nodal variables to export.
    elementVariables (list or None): Elemental variables to export.
    conditionVariables (list or None): Condition variables to export.
    exportElements (bool): Export elements.
    exportConditions (bool): Export conditions.
    **kwargs: Passed directly to plotter.add_mesh().

    Returns:
    pyvista.Plotter: Configured plotter — call .show() or .screenshot() next.
    """
    import pyvista as pv

    if nodalVariables is None:
        nodalVariables = []
    if elementVariables is None:
        elementVariables = []
    if conditionVariables is None:
        conditionVariables = []

    # Ensure the display variable and warp variable are included in the export
    extraNodal = []
    if name is not None and name not in nodalVariables:
        extraNodal.append(name)
    if (warpByVector is not None
            and warpByVector not in nodalVariables
            and warpByVector not in extraNodal):
        extraNodal.append(warpByVector)
    allNodal = list(nodalVariables) + extraNodal

    grid = ModelPartToPyVista(
        modelPart,
        useDeformedConfiguration=useDeformedConfiguration,
        nodalVariables=allNodal,
        elementVariables=list(elementVariables),
        conditionVariables=list(conditionVariables),
        exportElements=exportElements,
        exportConditions=exportConditions
    )

    if isinstance(grid, pv.MultiBlock):
        grid = grid["elements"] if "elements" in grid.keys() else grid[0]

    if plotter is None:
        if theme is not None:
            pv.set_plot_theme(theme)
        plotter = pv.Plotter(off_screen=offScreen or notebook)

    if grid is None or grid.n_points == 0:
        return plotter

    # Warp by vector field
    warpName = _GetVarName(warpByVector) if warpByVector is not None else None
    mesh = grid
    if warpName is not None and warpName in grid.point_data:
        if showUndeformed:
            plotter.add_mesh(grid, show_edges=False, opacity=0.2, line_width=lineWidth)
        mesh = grid.warp_by_vector(warpName, factor=factor)

    # Identify and extract the scalar to display
    varName = _GetVarName(name) if name is not None else None
    scalarsName = None
    ndim = 1

    if varName is not None:
        arr = None
        dataType = None
        if varName in mesh.point_data:
            arr      = mesh.point_data[varName]
            dataType = "point"
        elif varName in mesh.cell_data:
            arr      = mesh.cell_data[varName]
            dataType = "cell"

        if arr is not None:
            scalar1d, ndim = _ExtractComponent(arr, component, varName)

            # For scalar arrays (1-D) with component 0, use the array in-place;
            # for multi-component or magnitude, store a derived scalar under a temp key.
            if arr.ndim == 1 and (component == 0 or component is None):
                scalarsName = varName
            else:
                scalarsName = f"_kratos_{varName}_{component}"
                if dataType == "point":
                    mesh.point_data[scalarsName] = scalar1d
                else:
                    mesh.cell_data[scalarsName] = scalar1d

            if label is None:
                label = _BuildScalarLabel(varName, component, ndim)

    # Build scalar bar args dict
    sbArgs = None
    if scalarsName is not None and label is not None:
        sbArgs = {"title": label, "vertical": scalarBarVertical}
        if scalarBarArgs:
            sbArgs.update(scalarBarArgs)

    # Surface extraction (for quadratic / Lagrange cells)
    renderMesh = mesh
    if nonlinearSubdivision > 1 or extractSurface:
        nlSubd = max(1, nonlinearSubdivision)
        try:
            renderMesh = mesh.extract_surface(
                algorithm=None, nonlinear_subdivision=nlSubd
            )
            edges = (
                mesh.separate_cells()
                .extract_surface(algorithm=None, nonlinear_subdivision=1)
                .extract_feature_edges()
            )
            edgeActor = plotter.add_mesh(edges, color=edgeColor, line_width=lineWidth)
            edgeActor.mapper.SetResolveCoincidentTopologyToPolygonOffset()
        except Exception:
            renderMesh = mesh

        addKw = dict(scalars=scalarsName, cmap=cmap, show_edges=False, line_width=lineWidth)
        if nlSubd > 1:
            addKw["smooth_shading"]   = smoothShading
            addKw["split_sharp_edges"] = splitSharpEdges
        if sbArgs is not None:
            addKw["scalar_bar_args"] = sbArgs
        addKw.update(kwargs)
        plotter.add_mesh(renderMesh, **addKw)
    else:
        addKw = dict(
            scalars=scalarsName,
            cmap=cmap,
            show_edges=showEdges,
            line_width=lineWidth,
            edge_color=edgeColor,
        )
        if sbArgs is not None:
            addKw["scalar_bar_args"] = sbArgs
        addKw.update(kwargs)
        plotter.add_mesh(renderMesh, **addKw)

    # Camera preset
    if view == "default":
        if grid.n_points > 0 and np.allclose(grid.points[:, 2], 0.0):
            view = "xy"
            plotter.enable_parallel_projection()
        else:
            view = None
            plotter.camera.elevation = -15
            plotter.camera.azimuth   = -100
    if view is not None:
        plotter.camera_position = view

    if addAxes:
        plotter.add_axes()

    return plotter


def ScreenshotModelPart(
    modelPart,
    filename="modelpart.png",
    transparentBackground=None,
    scale=None,
    **kwargs
):
    """Render a Kratos ModelPart off-screen and save a screenshot to *filename*.

    All keyword arguments are forwarded to PlotModelPart (name, component,
    cmap, warpByVector, nodalVariables, …).

    Parameters:
    modelPart (KM.ModelPart): ModelPart to render.
    filename (str or None): Output image path. Pass None to return a numpy
        image array instead of writing a file.
    transparentBackground (bool or None): Transparent background flag.
    scale (int or None): Resolution scale factor.

    Returns:
    numpy.ndarray or None: Image array when filename is None, otherwise None.
    """
    p = PlotModelPart(modelPart, offScreen=True, **kwargs)
    return p.screenshot(filename, transparent_background=transparentBackground, scale=scale)


# ==============================================================================
# File I/O
# ==============================================================================

def SaveModelPart(
    modelPart,
    filename,
    useDeformedConfiguration=False,
    nodalVariables=[],
    elementVariables=[],
    conditionVariables=[]
):
    """Save a Kratos ModelPart to a VTK/VTU/VTM file via PyVista."""
    grid = ModelPartToPyVista(
        modelPart,
        useDeformedConfiguration=useDeformedConfiguration,
        nodalVariables=nodalVariables,
        elementVariables=elementVariables,
        conditionVariables=conditionVariables,
        exportElements=len(modelPart.Elements) > 0,
        exportConditions=len(modelPart.Conditions) > 0
    )
    grid.save(filename)


# ==============================================================================
# Advanced post-processing utilities
# ==============================================================================

def ComputeWarpedMesh(
    modelPart,
    vectorVariable,
    factor=1.0,
    nodalVariables=[],
    elementVariables=[],
    conditionVariables=[]
):
    """Generate a deformed mesh by warping nodal positions by a vector variable.

    Parameters:
    modelPart (KM.ModelPart): Source model part.
    vectorVariable (KM.Variable or str): Nodal vector variable (e.g. KM.DISPLACEMENT).
    factor (float): Warp scaling factor.
    nodalVariables (list): Additional nodal variables to carry through.
    elementVariables (list): Elemental variables to include.
    conditionVariables (list): Condition variables to include.

    Returns:
    pyvista.UnstructuredGrid or pyvista.MultiBlock: Warped mesh.
    """
    import pyvista as pv

    grid = ModelPartToPyVista(
        modelPart,
        useDeformedConfiguration=False,
        nodalVariables=nodalVariables + [vectorVariable],
        elementVariables=elementVariables,
        conditionVariables=conditionVariables,
        exportElements=len(modelPart.Elements) > 0,
        exportConditions=len(modelPart.Conditions) > 0
    )
    varName = _GetVarName(vectorVariable)

    if isinstance(grid, pv.MultiBlock):
        for key in grid.keys():
            if grid[key].n_points > 0 and varName in grid[key].point_data:
                grid[key] = grid[key].warp_by_vector(varName, factor=factor)
        return grid
    else:
        if grid.n_points > 0 and varName in grid.point_data:
            return grid.warp_by_vector(varName, factor=factor)
        return grid


def CreateExtractedSurface(modelPart, nodalVariables=[], elementVariables=[]):
    """Extract the outer surface of a 3-D Kratos ModelPart as a PyVista PolyData mesh.

    Useful for rendering volumetric meshes as smooth closed surfaces and for
    surface-based operations (feature-edge highlighting, smooth shading, etc.).

    Parameters:
    modelPart (KM.ModelPart): Source model part.
    nodalVariables (list): Nodal variables to transfer to the surface mesh.
    elementVariables (list): Elemental variables to transfer.

    Returns:
    pyvista.PolyData: Extracted surface mesh.
    """
    import pyvista as pv

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
        return pv.PolyData()

    return grid.extract_surface()


def CreateOrthogonalSlices(modelPart, x=None, y=None, z=None, nodalVariables=[], elementVariables=[]):
    """Extract orthogonal slices from a 3-D Kratos ModelPart domain.

    Parameters:
    modelPart (KM.ModelPart): Source model part.
    x (float or None): X-coordinate of the YZ-plane slice (midpoint if None).
    y (float or None): Y-coordinate of the XZ-plane slice (midpoint if None).
    z (float or None): Z-coordinate of the XY-plane slice (midpoint if None).
    nodalVariables (list): Nodal variables to include.
    elementVariables (list): Elemental variables to include.

    Returns:
    pyvista.MultiBlock: Three orthogonal slice surfaces.
    """
    import pyvista as pv

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
    """Generate contour surfaces (isosurfaces) from a nodal scalar variable.

    Parameters:
    modelPart (KM.ModelPart): Source model part.
    variable (KM.Variable or str): Nodal scalar variable to contour.
    valuesOrNumber (int or list): Number of evenly-spaced contour levels, or
        an explicit list of iso-values.
    nodalVariables (list): Additional nodal variables to include.

    Returns:
    pyvista.PolyData: Isosurface mesh.
    """
    import pyvista as pv

    varName = _GetVarName(variable)
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
        return pv.PolyData()

    return grid.contour(isosurfaces=valuesOrNumber, scalars=varName)


def CreateStreamlines(
    modelPart,
    velocityVariable="VELOCITY",
    sourceCenter=None,
    sourceRadius=None,
    nPoints=100,
    integrationDirection="both"
):
    """Generate fluid streamlines through a nodal vector velocity field.

    Parameters:
    modelPart (KM.ModelPart): Source model part.
    velocityVariable (KM.Variable or str): Nodal vector variable for the flow field.
    sourceCenter (list or None): [x, y, z] seed sphere centre (mesh centroid if None).
    sourceRadius (float or None): Seed sphere radius (10% of diagonal if None).
    nPoints (int): Number of seed points on the sphere.
    integrationDirection (str): 'forward', 'backward', or 'both'.

    Returns:
    pyvista.PolyData: Streamline tubes.
    """
    import pyvista as pv

    varName = _GetVarName(velocityVariable)
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
        return pv.PolyData()

    bounds = grid.bounds
    if sourceCenter is None:
        sourceCenter = [
            0.5 * (bounds[0] + bounds[1]),
            0.5 * (bounds[2] + bounds[3]),
            0.5 * (bounds[4] + bounds[5])
        ]
    if sourceRadius is None:
        diag = np.sqrt(
            (bounds[1] - bounds[0]) ** 2
            + (bounds[3] - bounds[2]) ** 2
            + (bounds[5] - bounds[4]) ** 2
        )
        sourceRadius = 0.1 * diag

    return grid.streamlines(
        vectors=varName,
        source_center=sourceCenter,
        source_radius=sourceRadius,
        n_points=nPoints,
        integration_direction=integrationDirection
    )


def CreateThresholdedMesh(
    modelPart,
    variable,
    thresholdValue,
    thresholdType="above",
    nodalVariables=[],
    elementVariables=[]
):
    """Extract mesh cells that satisfy a scalar threshold condition.

    Parameters:
    modelPart (KM.ModelPart): Source model part.
    variable (KM.Variable or str): Scalar variable to threshold.
    thresholdValue (float or tuple): Threshold value for 'above'/'below', or
        (lower, upper) pair for 'between'.
    thresholdType (str): 'above', 'below', or 'between'.
    nodalVariables (list): Additional nodal variables to include.
    elementVariables (list): Elemental variables to include.

    Returns:
    pyvista.UnstructuredGrid: Filtered mesh subset.
    """
    import pyvista as pv

    varName = _GetVarName(variable)
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

    preference = "point" if varName in grid.point_data else "cell"

    if thresholdType == "above":
        return grid.threshold(value=thresholdValue, scalars=varName, method="upper",
                              preference=preference)
    elif thresholdType == "below":
        return grid.threshold(value=thresholdValue, scalars=varName, method="lower",
                              preference=preference)
    elif thresholdType == "between":
        if not isinstance(thresholdValue, (list, tuple)) or len(thresholdValue) != 2:
            raise ValueError(
                "For 'between' thresholdType, thresholdValue must be a (lower, upper) pair."
            )
        return grid.threshold(value=thresholdValue, scalars=varName, preference=preference)
    else:
        raise ValueError(
            f"Unknown thresholdType '{thresholdType}'. Supported: 'above', 'below', 'between'."
        )


def CreateVectorGlyphs(
    modelPart,
    vectorVariable,
    scaleFactor=1.0,
    glyphType="arrow",
    nodalVariables=[]
):
    """Generate oriented 3-D glyphs at mesh nodes scaled by a vector variable.

    Parameters:
    modelPart (KM.ModelPart): Source model part.
    vectorVariable (KM.Variable or str): Nodal vector variable for orientation and scale.
    scaleFactor (float): Global glyph size scale factor.
    glyphType (str): Glyph geometry: 'arrow', 'cone', 'sphere', or 'cylinder'.
    nodalVariables (list): Additional nodal variables to include.

    Returns:
    pyvista.PolyData: Glyph mesh.
    """
    import pyvista as pv

    varName = _GetVarName(vectorVariable)
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

    _glyph_map = {
        "arrow":    pv.Arrow(),
        "cone":     pv.Cone(),
        "sphere":   pv.Sphere(),
        "cylinder": pv.Cylinder(),
    }
    geom = _glyph_map.get(glyphType, pv.Arrow())

    glyphs = grid.glyph(orient=varName, scale=varName, factor=scaleFactor, geom=geom)
    if "GlyphVector" in glyphs.point_data:
        glyphs.point_data[varName] = glyphs.point_data["GlyphVector"]
    if "GlyphScale" in glyphs.point_data:
        glyphs.point_data[varName + "_magnitude"] = glyphs.point_data["GlyphScale"]
    return glyphs


def CreateClippedMesh(
    modelPart,
    normal=[0.0, 0.0, 1.0],
    origin=None,
    nodalVariables=[],
    elementVariables=[]
):
    """Clip a Kratos ModelPart by a plane defined by an origin and normal vector.

    Parameters:
    modelPart (KM.ModelPart): Source model part.
    normal (list): 3-D normal vector of the clipping plane.
    origin (list or None): 3-D origin of the clipping plane (mesh centre if None).
    nodalVariables (list): Nodal variables to include.
    elementVariables (list): Elemental variables to include.

    Returns:
    pyvista.UnstructuredGrid: Clipped mesh.
    """
    import pyvista as pv

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
