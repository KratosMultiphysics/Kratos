"""Create a file containing xdmf metadata for results stored in HDF5.
"""
import KratosMultiphysics
import KratosMultiphysics.HDF5Application as KratosHDF5
import os, sys, h5py, xdmf


def GenerateXdmfConnectivities(file_name):
    with h5py.File(file_name, "r") as h5py_file:
        has_xdmf = ("Xdmf" in h5py_file.get('/ModelData').keys())
    if not has_xdmf:
        KratosHDF5.HDF5XdmfConnectivitiesWriterProcess(file_name, "/ModelData").Execute()


def GetSpatialGrid(h5py_file):
    elements_path = "/ModelData/Xdmf/Elements/"
    coordinates_path = "/ModelData/Nodes/Local/Coordinates"
    spatial_grid = xdmf.SpatialGrid()
    coords_data = xdmf.HDF5UniformDataItem(h5py_file.get(coordinates_path))
    geom = xdmf.Geometry(coords_data)
    elems_group = h5py_file.get(elements_path)
    for name in elems_group.keys():
        if isinstance(elems_group[name], h5py.Group):
            single_elem_group = elems_group.get(name)
            dim = single_elem_group.attrs["Dimension"]
            num_points = single_elem_group.attrs["NumberOfNodes"]
            cell_type = xdmf.TopologyCellType(dim, num_points)
            connectivity_data = xdmf.HDF5UniformDataItem(h5py_file.get(elements_path + '/' + name + '/Connectivities'))
            topology = xdmf.UniformMeshTopology(cell_type, connectivity_data)
            spatial_grid.add_grid(xdmf.UniformGrid(name, geom, topology))
    return spatial_grid


def GetNodalResults(h5py_file):
    nodal_results_path = "/ResultsData/NodalResults"
    results = []
    results_group = h5py_file.get(nodal_results_path)
    for variable_name in results_group.keys():
        if isinstance(results_group[variable_name], h5py.Dataset):
            data = xdmf.HDF5UniformDataItem(results_group.get(variable_name))
            results.append(xdmf.NodalSolutionStepData(variable_name, data))
    return results


def GetListOfTimeLabels(file_name):
    list_of_file_names = []
    time_prefix = file_name.replace(".h5", "-")
    for name in os.listdir():
        if name.find(time_prefix) == 0:
            list_of_file_names.append(name)
    list_of_time_labels = []
    for name in list_of_file_names:
        list_of_time_labels.append(name.replace(".h5", "")[len(time_prefix):])
    list_of_time_labels.sort(key=float)
    return list_of_time_labels


def main():
    file_name = sys.argv[1]
    temporal_grid = xdmf.TemporalGrid()
    GenerateXdmfConnectivities(file_name)
    # Get the initial spatial grid from the base file.
    with h5py.File(file_name, "r") as h5py_file:
        current_spatial_grid = GetSpatialGrid(h5py_file)
    for current_time in GetListOfTimeLabels(file_name):
        current_file_name = file_name.replace(".h5", "-" + current_time + ".h5")
        # Check if the current file has mesh information.
        with h5py.File(current_file_name, "r") as h5py_file:
            has_mesh = ("ModelData" in h5py_file.keys())
        if has_mesh:
            GenerateXdmfConnectivities(current_file_name)
        with h5py.File(current_file_name, "r") as h5py_file:
            if has_mesh:
                # Update current spatial grid
                current_spatial_grid = GetSpatialGrid(h5py_file)
            # Initialize the current grid with the spatial grid.
            current_grid = xdmf.SpatialGrid()
            for grid in current_spatial_grid.grids:
                current_grid.add_grid(xdmf.UniformGrid(grid.name, grid.geometry, grid.topology))
            # Add the (time-dependent) results.
            for nodal_result in GetNodalResults(h5py_file):
                current_grid.add_attribute(nodal_result)
        # Add the current grid to the temporal grid.
        temporal_grid.add_grid(xdmf.Time(current_time), current_grid)
    # Create the domain.
    domain = xdmf.Domain(temporal_grid)
    # Write.
    xdmf_file_name = file_name.replace(".h5", ".xdmf")
    xdmf.ET.ElementTree(xdmf.Xdmf(domain).create_xml_element()).write(xdmf_file_name)


if __name__ == '__main__':
    main()
