import xml.etree.ElementTree as ET
import h5py
from pathlib import Path

import KratosMultiphysics
import KratosMultiphysics.HDF5Application as KratosHDF5
from KratosMultiphysics.HDF5Application.core.pattern import IdentifyPattern
from KratosMultiphysics.HDF5Application.core.xdmf import SpatialGrid
from KratosMultiphysics.HDF5Application.core.xdmf import HDF5UniformDataItem
from KratosMultiphysics.HDF5Application.core.xdmf import Geometry
from KratosMultiphysics.HDF5Application.core.xdmf import TopologyCellType
from KratosMultiphysics.HDF5Application.core.xdmf import UniformMeshTopology
from KratosMultiphysics.HDF5Application.core.xdmf import UniformGrid
from KratosMultiphysics.HDF5Application.core.xdmf import TemporalGrid
from KratosMultiphysics.HDF5Application.core.xdmf import Time
from KratosMultiphysics.HDF5Application.core.xdmf import Domain
from KratosMultiphysics.HDF5Application.core.xdmf import Xdmf
from KratosMultiphysics.HDF5Application.core.xdmf import EntityData
from KratosMultiphysics.HDF5Application.core.dataset_generator import DataSetGenerator

def GetPatternDetailsFromFileName(file_name_path: Path) -> str:
    file_relative_path = str(file_name_path.relative_to(Path(".")))
    if file_relative_path.endswith(".h5"):
        file_relative_path = file_relative_path[:-1]
        pattern, tag_type_dict = IdentifyPattern(file_relative_path)
        return pattern + "5", tag_type_dict
    else:
        return IdentifyPattern(file_relative_path)

def GetSpatialGridPaths(h5_file: h5py.File, list_of_model_data_paths: 'list[str]', current_path: str) -> None:
    current_data = h5_file[current_path]
    if isinstance(current_data, h5py.Group):
        if "__model_part_name" in current_data.attrs.keys():
            list_of_model_data_paths.append(current_path)
        else:
            for itr in current_data:
                GetSpatialGridPaths(h5_file, list_of_model_data_paths, f"{current_path}/{str(itr)}")

def __GetModelPartUniformGrids(h5_group: h5py.Group, grid_name: str, coordinates: Geometry, is_root: bool, container_type: KratosMultiphysics.Globals.DataLocation) -> 'list[UniformGrid]':
    grids: 'list[UniformGrid]' = []
    for name, entity_group in h5_group.items():
        if isinstance(entity_group, h5py.Group):
            if "NumberOfNodes" in entity_group.attrs and "WorkingSpaceDimension" in entity_group.attrs and "Connectivities" in entity_group.keys():
                cell_type = TopologyCellType(entity_group.attrs["WorkingSpaceDimension"], entity_group.attrs["NumberOfNodes"])
                connectivities = HDF5UniformDataItem(entity_group["Connectivities"])
                topology = UniformMeshTopology(cell_type, connectivities)
                grids.append(UniformGrid(f"{grid_name}.{name}", coordinates, topology, is_root, container_type))
                KratosMultiphysics.Logger.PrintInfo("XDMF", f"Added \"{grid_name}.{name}\" spatial grid from \"{h5_group.file.filename}:{h5_group.name}\".")
            else:
                raise RuntimeError(f"Either \"Connectivities\", \"NumberOfNodes\" or \"WorkingSpaceDimension\" are not defined for {h5_group.name}.")
    return grids

def GenerateUniformGridsForSubModelParts(h5_file: h5py.File, grid_path: str, grid_name: str, coordinates: Geometry, spatial_grid: SpatialGrid):
    grid_data = h5_file[grid_path]

    for name, h5_object in grid_data.items():
        if isinstance(h5_object, h5py.Dataset) and name == "Points":
            cell_type = TopologyCellType(3, 1)
            points = HDF5UniformDataItem(h5_object)
            topology = UniformMeshTopology(cell_type, points)
            spatial_grid.AddGrid(UniformGrid(f"{grid_name}.Points", coordinates, topology, False, KratosMultiphysics.Globals.DataLocation.NodeNonHistorical))
            KratosMultiphysics.Logger.PrintInfo("XDMF", f"Added \"{grid_name}.{name}\" spatial grid from \"{h5_object.file.filename}:{h5_object.name}\".")
        elif isinstance(h5_object, h5py.Group):
            if name == "Conditions":
                list(map(spatial_grid.AddGrid, __GetModelPartUniformGrids(h5_object, grid_name, coordinates, False, KratosMultiphysics.Globals.DataLocation.Condition)))
            elif name == "Elements":
                list(map(spatial_grid.AddGrid, __GetModelPartUniformGrids(h5_object, grid_name, coordinates, False, KratosMultiphysics.Globals.DataLocation.Element)))
            else:
                GenerateUniformGridsForSubModelParts(h5_file, f"{grid_path}/{name}", f"{grid_name}.{name}", coordinates, spatial_grid)

def CreateXdmfSpatialGrid(h5_file_name: str, model_data_path: str) -> SpatialGrid:
    """Generate a spatial grid from data in h4 file name and model data path.

    Args:
        h5_file_name (str): hdf5 file name.
        model_data_path (str): model data prefix.

    Raises:
        RuntimeError: If the model data path is not found.
        RuntimeError: IF the specified model_data_path does not contain model data.

    Returns:
        SpatialGrid: Spatial grid containing root model part nodes, conditions, elements and sub-model parts.
    """
    h5_file = h5py.File(h5_file_name, "r")

    if model_data_path not in h5_file.keys():
        raise RuntimeError(f"The model_data_path = \"{model_data_path}\" not found in {h5_file_name}.")

    model_data = h5_file[model_data_path]

    if not "Xdmf" in model_data.keys():
        # temporarily close existing file so following operation can open
        # and add Xdmf content
        h5_file.close()
        # try generating the xdmf data.
        KratosHDF5.HDF5XdmfConnectivitiesWriterOperation(h5_file_name, [model_data_path]).Execute()
        # again open the h5_file
        h5_file = h5py.File(h5_file_name, "r")

    spatial_grid = SpatialGrid()
    model_data = h5_file[model_data_path]

    if "__model_part_name" not in model_data.attrs.keys():
        raise RuntimeError(f"The attribute \"__model_part_name\" not found in \"{h5_file_name}:{model_data_path}\".")

    model_part_name = ''.join([chr(i) for i in model_data.attrs["__model_part_name"]])

    coordinates = Geometry(HDF5UniformDataItem(h5_file[f"{model_data_path}/Nodes/Local/Coordinates"]))

    if "Xdmf" in model_data.keys():
        xdmf_data = model_data["Xdmf"]

        # add root model part entities
        if "Conditions" in xdmf_data.keys():
            list(map(spatial_grid.AddGrid, __GetModelPartUniformGrids(xdmf_data["Conditions"], model_part_name, coordinates, True, KratosMultiphysics.Globals.DataLocation.Condition)))

        if "Elements" in xdmf_data.keys():
            list(map(spatial_grid.AddGrid, __GetModelPartUniformGrids(xdmf_data["Elements"], model_part_name, coordinates, True, KratosMultiphysics.Globals.DataLocation.Element)))

        # now add sub model part entitites
        if "SubModelParts" in xdmf_data.keys():
            GenerateUniformGridsForSubModelParts(h5_file, f"{model_data_path}/Xdmf/SubModelParts", model_part_name, coordinates, spatial_grid)
    else:
        KratosMultiphysics.Logger.PrintInfo("XDMF", f"No Xdmf data found in {h5_file_name}:{model_data_path}.")

    return spatial_grid

def WriteMeshToXdmf(h5_file_name: str, output_file_name: str, model_data_prefix: str = "/") -> None:
    list_of_model_part_prefixes = []
    if model_data_prefix == "/":
        with h5py.File(h5_file_name, "r") as h5_file:
            GetSpatialGridPaths(h5_file, list_of_model_part_prefixes, "/")
    else:
        list_of_model_part_prefixes = [model_data_prefix]

    temporal_grid = TemporalGrid()
    for spatial_grid_path in list_of_model_part_prefixes:
        spatial_grid = CreateXdmfSpatialGrid(h5_file_name, spatial_grid_path)
        temporal_grid.AddGrid(Time(0.0), spatial_grid)

    domain = Domain(temporal_grid)
    xdmf = Xdmf(domain)
    ET.ElementTree(xdmf.CreateXmlElement()).write(output_file_name)

def WriteDataSetsToXdmf(dataset_generator: DataSetGenerator, output_file_name: str) -> None:
    """Write the given datasets by the generator to the specified output file name.

    This method can write given datasets in the generator to xdmf. The generator does
    not have to be yielding datasets in a specific order. The datasets will
    be ordered based on the temporal value each yielding provides in this method.

    Since the HDF5 dataset file paths are written relative to the point where
    this method is called, the @ref output_file_name should not have any relative
    paths in it.

    Args:
        dataset_generator (DataSetGenerator): tuple[float, EntityData] generator
        output_file_name (str): output xdmf filename

    Raises:
        RuntimeError: If output_file_name includes a relative path.
    """
    if "/" in output_file_name or "\\"  in output_file_name:
        raise RuntimeError(f"Output file name should not include a relative path [ output_file_name = \"{output_file_name}\" ].")

    temporal_information: 'dict[float, list[EntityData]]' = {}

    for temporal_value, dataset in dataset_generator.Iterate():
        if temporal_value not in temporal_information.keys():
            temporal_information[temporal_value]: 'list[EntityData]' = []
        temporal_information[temporal_value].append(dataset)

    # now we are done with the dataset_generator. So we can safely generate the spatial grids.
    # otherwise, if the dataset_generator also keeps the same file open, it may conflict with
    # the xdmf data generation which is required when generating spatial grids.
    temporal_grid = TemporalGrid()
    spatial_grids: 'dict[str, SpatialGrid]' = {}

    for temporal_value in sorted(list(temporal_information.keys())):
        # there is no sense in haiving mutiple meshes in the same time step. Hence mesh data
        # is retrieved from the first dataset.
        datasets = temporal_information[temporal_value]

        mesh_location = datasets[0].mesh_location
        if mesh_location not in spatial_grids.keys():
            mesh_location_data = mesh_location.split(":")
            spatial_grids[mesh_location] = CreateXdmfSpatialGrid(mesh_location_data[0], mesh_location_data[1])

        # create a copy of the spatial grid
        current_spatial_grid = SpatialGrid()
        for grid in spatial_grids[mesh_location].grids:
            current_spatial_grid.AddGrid(UniformGrid(grid.name, grid.geometry, grid.topology, grid.is_root, grid.container_type))

        for dataset in datasets:
            current_spatial_grid.AddAttribute(dataset)

        temporal_grid.AddGrid(Time(temporal_value), current_spatial_grid)
        KratosMultiphysics.Logger.PrintInfo("XDMF", f"Written data for temporal value = {temporal_value}.")

    domain = Domain(temporal_grid)
    xdmf = Xdmf(domain)
    ET.ElementTree(xdmf.CreateXmlElement()).write(output_file_name)