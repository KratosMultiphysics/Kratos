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

def _create_unstructured_grid(model_part, entities, use_deformed_configuration, nodal_variables, entity_variables):
    """Helper function to convert a collection of entities (Elements or Conditions) into a pyvista.UnstructuredGrid."""
    if len(entities) == 0:
        return pv.UnstructuredGrid()

    num_nodes = len(model_part.Nodes)
    if num_nodes == 0:
        return pv.UnstructuredGrid()

    # 1. Build points array
    points = np.zeros((num_nodes, 3), dtype=np.float64)
    node_id_to_idx = {}
    for idx, node in enumerate(model_part.Nodes):
        node_id_to_idx[node.Id] = idx
        if use_deformed_configuration:
            points[idx, 0] = node.X
            points[idx, 1] = node.Y
            points[idx, 2] = node.Z
        else:
            points[idx, 0] = node.X0
            points[idx, 1] = node.Y0
            points[idx, 2] = node.Z0

    # 2. Build cells and cell types
    cells = []
    cell_types = []

    for entity in entities:
        geom = entity.GetGeometry()
        geom_type = geom.GetGeometryType()
        vtk_type = GEOMETRY_TYPE_TO_VTK_CELL_TYPE.get(geom_type)

        if vtk_type is None:
            continue

        cell_nodes = [node_id_to_idx[node.Id] for node in geom]
        cells.append(len(cell_nodes))
        cells.extend(cell_nodes)
        cell_types.append(vtk_type)

    if len(cells) == 0:
        return pv.UnstructuredGrid()

    grid = pv.UnstructuredGrid(np.array(cells, dtype=np.int32), np.array(cell_types, dtype=np.uint8), points)

    # 3. Add point (nodal) variables
    for var in nodal_variables:
        if isinstance(var, str):
            var = KM.KratosGlobals.GetVariable(var)

        first_node = next(iter(model_part.Nodes))
        is_historical = first_node.HasSolutionStepValue(var)
        is_non_historical = first_node.Has(var)

        if not (is_historical or is_non_historical):
            continue

        val = first_node.GetSolutionStepValue(var) if is_historical else first_node.GetValue(var)

        if isinstance(val, (float, int)):
            arr = np.zeros(num_nodes, dtype=np.float64)
            for idx, node in enumerate(model_part.Nodes):
                arr[idx] = node.GetSolutionStepValue(var) if is_historical else node.GetValue(var)
            grid.point_data[var.Name()] = arr
        elif isinstance(val, (KM.Array3, KM.Vector, list)) or type(val).__name__ == "Vector":
            arr = np.zeros((num_nodes, 3), dtype=np.float64)
            for idx, node in enumerate(model_part.Nodes):
                v = node.GetSolutionStepValue(var) if is_historical else node.GetValue(var)
                arr[idx, 0] = v[0]
                arr[idx, 1] = v[1]
                arr[idx, 2] = v[2]
            grid.point_data[var.Name()] = arr

    # 4. Add cell (entity) variables
    for var in entity_variables:
        if isinstance(var, str):
            var = KM.KratosGlobals.GetVariable(var)

        val_type = None
        for entity in entities:
            geom_type = entity.GetGeometry().GetGeometryType()
            if geom_type in GEOMETRY_TYPE_TO_VTK_CELL_TYPE and entity.Has(var):
                val_type = type(entity.GetValue(var))
                break

        if val_type is None:
            continue

        if issubclass(val_type, (float, int)):
            arr = np.zeros(len(cell_types), dtype=np.float64)
            cell_idx = 0
            for entity in entities:
                geom_type = entity.GetGeometry().GetGeometryType()
                if geom_type not in GEOMETRY_TYPE_TO_VTK_CELL_TYPE:
                    continue
                if entity.Has(var):
                    arr[cell_idx] = entity.GetValue(var)
                cell_idx += 1
            grid.cell_data[var.Name()] = arr
        elif issubclass(val_type, (KM.Array3, KM.Vector, list)) or val_type.__name__ == "Vector":
            arr = np.zeros((len(cell_types), 3), dtype=np.float64)
            cell_idx = 0
            for entity in entities:
                geom_type = entity.GetGeometry().GetGeometryType()
                if geom_type not in GEOMETRY_TYPE_TO_VTK_CELL_TYPE:
                    continue
                if entity.Has(var):
                    v = entity.GetValue(var)
                    arr[cell_idx, 0] = v[0]
                    arr[cell_idx, 1] = v[1]
                    arr[cell_idx, 2] = v[2]
                cell_idx += 1
            grid.cell_data[var.Name()] = arr

    return grid


def model_part_to_pyvista(model_part, use_deformed_configuration=False, nodal_variables=[], element_variables=[], condition_variables=[], export_elements=True, export_conditions=False):
    """Converts a Kratos ModelPart into a pyvista.UnstructuredGrid or pyvista.MultiBlock mesh.

    Parameters:
    model_part (KM.ModelPart): The Kratos ModelPart to convert.
    use_deformed_configuration (bool): If True, use current coordinates (X, Y, Z). If False, use initial coordinates (X0, Y0, Z0).
    nodal_variables (list): List of KM.Variable or string names for nodal data to export.
    element_variables (list): List of KM.Variable or string names for element (cell) data to export.
    condition_variables (list): List of KM.Variable or string names for condition (cell) data to export.
    export_elements (bool): Whether to export elements.
    export_conditions (bool): Whether to export conditions.

    Returns:
    pyvista.UnstructuredGrid or pyvista.MultiBlock: The converted PyVista mesh.
    """
    if export_elements and export_conditions:
        elem_grid = _create_unstructured_grid(
            model_part,
            model_part.Elements,
            use_deformed_configuration,
            nodal_variables,
            element_variables
        )
        cond_grid = _create_unstructured_grid(
            model_part,
            model_part.Conditions,
            use_deformed_configuration,
            nodal_variables,
            condition_variables
        )
        blocks = pv.MultiBlock()
        blocks["elements"] = elem_grid
        blocks["conditions"] = cond_grid
        return blocks
    elif export_elements:
        return _create_unstructured_grid(
            model_part,
            model_part.Elements,
            use_deformed_configuration,
            nodal_variables,
            element_variables
        )
    elif export_conditions:
        return _create_unstructured_grid(
            model_part,
            model_part.Conditions,
            use_deformed_configuration,
            nodal_variables,
            condition_variables
        )
    else:
        return pv.UnstructuredGrid()


def plot_model_part(model_part, scalars=None, use_deformed_configuration=False, **kwargs):
    """Plot a Kratos ModelPart using PyVista's interactive plotting engine."""
    nodal_vars = [scalars] if scalars is not None else []
    grid = model_part_to_pyvista(
        model_part,
        use_deformed_configuration=use_deformed_configuration,
        nodal_variables=nodal_vars,
        export_elements=len(model_part.Elements) > 0,
        export_conditions=len(model_part.Conditions) > 0 and len(model_part.Elements) == 0
    )

    if scalars is not None:
        scalars_name = scalars if isinstance(scalars, str) else scalars.Name()
        if isinstance(grid, pv.MultiBlock):
            for block in grid:
                if scalars_name in block.point_data:
                    block.set_active_scalars(scalars_name)
        else:
            if scalars_name in grid.point_data:
                grid.set_active_scalars(scalars_name)

    return grid.plot(**kwargs)


def save_model_part(model_part, filename, use_deformed_configuration=False, nodal_variables=[], element_variables=[], condition_variables=[]):
    """Save a Kratos ModelPart to a VTK/VTU/VTM file using PyVista's save mechanisms."""
    export_elements = len(model_part.Elements) > 0
    export_conditions = len(model_part.Conditions) > 0

    grid = model_part_to_pyvista(
        model_part,
        use_deformed_configuration=use_deformed_configuration,
        nodal_variables=nodal_variables,
        element_variables=element_variables,
        condition_variables=condition_variables,
        export_elements=export_elements,
        export_conditions=export_conditions
    )
    grid.save(filename)
