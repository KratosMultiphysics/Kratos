"""Create a file containing xdmf metadata for results stored in HDF5."""

import KratosMultiphysics
import KratosMultiphysics.HDF5Application as KratosHDF5
import os, sys, xdmf
import warnings
with warnings.catch_warnings():
    # suppressing an import-related warningfrom h5py
    # problem appears when using it in a test with python >=3.6
    warnings.simplefilter('ignore', category=ImportWarning)
    import h5py

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
    results = {}
    if "/ResultsData/NodalSolutionStepData" in h5py_file:
        AddNodalData(h5py_file.get("/ResultsData/NodalSolutionStepData"), results)
    if "/ResultsData/NodalDataValues" in h5py_file:
        AddNodalData(h5py_file.get("/ResultsData/NodalDataValues"), results)

    return list(results.values())

def AddNodalData(results_group, results):
    for variable_name in results_group.keys():
        if isinstance(results_group[variable_name], h5py.Dataset):
            if variable_name in results:
                raise ValueError('Nodal result "' + variable_name + '" is already defined.')
            data = xdmf.HDF5UniformDataItem(results_group.get(variable_name))
            results[variable_name] = xdmf.NodalData(variable_name, data)

def GetElementResults(h5py_file):
    element_results_path = "/ResultsData/ElementDataValues"
    results = []
    if not element_results_path in h5py_file:
        return results
    results_group = h5py_file.get(element_results_path)
    for variable_name in results_group.keys():
        if isinstance(results_group[variable_name], h5py.Dataset):
            data = xdmf.HDF5UniformDataItem(results_group.get(variable_name))
            results.append(xdmf.ElementSolutionStepData(variable_name, data))
    return results


def GetListOfTimeLabels(file_name):
    list_of_file_names = []
    path, file_name = os.path.split(file_name)
    if path == "": path = "." # os.listdir fails with empty path
    time_prefix = file_name.replace(".h5", "-")
    for name in os.listdir(path):
        if name.find(time_prefix) == 0:
            list_of_file_names.append(name)
    list_of_time_labels = []
    for name in list_of_file_names:
        list_of_time_labels.append(name.replace(".h5", "")[len(time_prefix):])
    list_of_time_labels.sort(key=float)
    return list_of_time_labels


def WriteXdmfFile(file_name):
    #todo(msandre): generalize to WriteXdmfFile(xdmf_file_name, list_of_h5_file_paths):
    temporal_grid = xdmf.TemporalGrid()
    GenerateXdmfConnectivities(file_name)
    # Get the initial spatial grid from the base file.
    with h5py.File(file_name, "r") as h5py_file:
        current_spatial_grid = GetSpatialGrid(h5py_file)
    for current_time in GetListOfTimeLabels(file_name):
        current_file_name = file_name.replace(".h5", "-" + current_time + ".h5")
        try:
            # Check if the current file has mesh information.
            with h5py.File(current_file_name, "r") as h5py_file:
                has_mesh = ("ModelData" in h5py_file.keys())
                has_data = ("/ResultsData" in h5py_file.keys())
        except OSError:
            # in case this file cannot be opened skip it
            # this can be the case if the file is already opened
            warn_msg  = 'No xdmf-data was written for file:\n"'
            warn_msg += current_file_name + '"'
            KratosMultiphysics.Logger.PrintWarning("XDMF-Writing", warn_msg)
            continue
        if not has_data:
            continue
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
            for element_result in GetElementResults(h5py_file):
                current_grid.add_attribute(element_result)
        # Add the current grid to the temporal grid.
        temporal_grid.add_grid(xdmf.Time(current_time), current_grid)
    # Create the domain.
    domain = xdmf.Domain(temporal_grid)
    # Write.
    raw_file_name = os.path.split(file_name)[1]
    xdmf_file_name = raw_file_name.replace(".h5", ".xdmf")
    xdmf.ET.ElementTree(xdmf.Xdmf(domain).create_xml_element()).write(xdmf_file_name)


if __name__ == '__main__':
    file_name = sys.argv[1]
    WriteXdmfFile(file_name)
