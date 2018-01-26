"""Create a file containing xdmf metadata for results stored in HDF5.

This module should be used when the simulation mesh is not changing in time.
"""
import KratosMultiphysics
import KratosMultiphysics.HDF5Application as KratosHDF5
import os, sys, h5py, xdmf

class CreateFixedMeshXdmfFileProcess(KratosMultiphysics.Process):

    _elements_path = "/ModelData/Xdmf/Elements/"

    _coordinates_path = "/ModelData/Nodes/Local/Coordinates"

    _nodal_results_path = "/ResultsData/NodalResults"

    def __init__(self, file_name):
        self._file_name = file_name

    def Execute(self):
        self._generate_xdmf_connectivities()
        spatial_grid = self._get_spatial_grid()
        temporal_grid = xdmf.TemporalGrid()
        list_of_times = self._get_list_of_time_labels()
        for current_time in list_of_times:
            # Build the current grid from the reference grid.
            current_grid = xdmf.SpatialGrid()
            for grid in spatial_grid.grids:
                current_grid.add_grid(xdmf.UniformGrid(grid.name, grid.geometry, grid.topology))
            # Add the (time-dependent) results.
            for nodal_result in self._get_nodal_results(current_time):
                current_grid.add_attribute(nodal_result)
            # Add the current grid to the temporal grid.
            temporal_grid.add_grid(xdmf.Time(current_time), current_grid)
        # Create the domain.
        domain = xdmf.Domain(temporal_grid)
        # Write.
        xdmf_file_name = self._file_name.replace(".h5", ".xdmf")
        xdmf.ET.ElementTree(xdmf.Xdmf(domain).create_xml_element()).write(xdmf_file_name)

    def _generate_xdmf_connectivities(self):
        with h5py.File(self._file_name, "r") as h5py_file:
            if not "Xdmf" in h5py_file.get('/ModelData').keys():
                KratosHDF5.HDF5XdmfConnectivitiesWriterProcess(self._file_name, "/ModelData").Execute()
        # Create the reference mesh's xml hierarchy.

    def _get_spatial_grid(self):
        spatial_grid = xdmf.SpatialGrid()
        with h5py.File(self._file_name, "r") as h5py_file:
            data_set = h5py_file.get(self._coordinates_path)
            coords_data = xdmf.HDF5UniformDataItem(data_set)
            geom = xdmf.Geometry(coords_data)
            elems_group = h5py_file.get(self._elements_path)
            for elem_name in elems_group.keys():
                single_elem_group = elems_group.get(elem_name)
                dim = single_elem_group.attrs["Dimension"]
                num_points = single_elem_group.attrs["NumberOfNodes"]
                cell_type = xdmf.TopologyCellType(dim, num_points)
                connectivity_data = xdmf.HDF5UniformDataItem(h5py_file.get(self._elements_path + '/' + elem_name + '/Connectivities'))
                topology = xdmf.UniformMeshTopology(cell_type, connectivity_data)
                spatial_grid.add_grid(xdmf.UniformGrid(elem_name, geom, topology))
        return spatial_grid

    def _get_nodal_results(self, current_time):
        results = []
        results_file_name = self._file_name.replace(".h5", "-" + current_time + ".h5")
        with h5py.File(results_file_name, "r") as h5py_file:
            results_group = h5py_file.get(self._nodal_results_path)
            for variable_name in results_group.keys():
                if variable_name != "Partition":
                    data = xdmf.HDF5UniformDataItem(results_group.get(variable_name))
                    results.append(xdmf.NodalSolutionStepData(variable_name, data))
        return results

    def _get_list_of_time_labels(self):
        list_of_file_names = []
        time_prefix = self._file_name.replace(".h5", "") + "-"
        for file_name in os.listdir():
            if file_name.find(time_prefix) == 0:
                list_of_file_names.append(file_name)
        list_of_time_labels = []
        for file_name in list_of_file_names:
            list_of_time_labels.append(file_name.replace(".h5", "")[len(time_prefix):])
        list_of_time_labels.sort(key=float)
        return list_of_time_labels

if __name__ == '__main__':
    file_name = sys.argv[1]
    CreateFixedMeshXdmfFileProcess(file_name).Execute()
