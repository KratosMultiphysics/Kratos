"""Generate xdmf metadata file from results stored in HDF5.

After running a Kratos simulation with HDF5 output on a fixed mesh, an Xdmf
metadata file can be generated for post-processing.
"""
import os, sys, h5py, copy, xdmf_utils

class XdmfFixedMeshPostProcess(object):
    """Constructs the XDMF hierarchy for fixed mesh results."""

    def __init__(self, model_part_name):
        self.model_part_name = model_part_name

    def Execute(self):
        xdmf_utils.WriteSortedCoordinates(self.model_part_name + ".h5", "/ModelData/Nodes/Local")
        # Create the reference mesh's xml hierarchy.
        with h5py.File(self.model_part_name + ".h5", "r") as h5py_file:
            fixed_mesh = xdmf_utils.KratosCollectionGrid(h5py_file, "/ModelData")
        temporal_grid = xdmf_utils.XdmfElement("Grid")
        temporal_grid.root.set("GridType", "Collection")
        temporal_grid.root.set("CollectionType", "Temporal")
        list_of_times = self._get_list_of_time_labels()
        for current_time in list_of_times:
            # Create the current mesh.
            current_mesh = copy.deepcopy(fixed_mesh)
            results_file_name = self.model_part_name + "-" + current_time + ".h5"
            self._add_results_to_children(results_file_name, current_mesh)
            grid_time = xdmf_utils.KratosTime(current_time)
            current_mesh.root.append(grid_time.root)
            temporal_grid.root.append(current_mesh.root)
        domain = xdmf_utils.XdmfElement("Domain")
        domain.root.append(temporal_grid.root)
        xdmf = xdmf_utils.XdmfElement("Xdmf")
        xdmf.root.set("Version", "3.0")
        xdmf.root.append(domain.root)
        xdmf.Write(self.model_part_name + ".xdmf")

    def _add_results_to_children(self, results_file_name, parent):
        with h5py.File(results_file_name, "r") as h5py_file:
            nodal_results = h5py_file.get('/ResultsData/NodalResults')
            for variable_name in nodal_results.keys():
                if variable_name == "Partition":
                    continue
                variable_path = '/ResultsData/NodalResults/' + variable_name
                variable_data = xdmf_utils.KratosNodalSolutionStepDataAttribute(h5py_file, variable_path)
                for child in parent.root:
                    child.append(variable_data.root)

    def _get_list_of_time_labels(self):
        list_of_file_names = []
        time_prefix = self.model_part_name + "-"
        for file_name in os.listdir():
            if file_name.find(time_prefix) == 0:
                list_of_file_names.append(file_name)
        list_of_time_labels = []
        for file_name in list_of_file_names:
            list_of_time_labels.append(file_name[len(time_prefix):-3])
        return list_of_time_labels

if __name__ == '__main__':
    model_part_file_name = sys.argv[1]
    model_part_name = model_part_file_name.split(sep=".")[0]
    XdmfFixedMeshPostProcess(model_part_name).Execute()
