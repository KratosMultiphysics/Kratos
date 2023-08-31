import xml.etree.ElementTree as ET
import typing
import h5py
import abc
from pathlib import Path

import KratosMultiphysics
from KratosMultiphysics.HDF5Application.core.pattern import PathPatternEntity
from KratosMultiphysics.HDF5Application.core.pattern import GetMachingEntitiesWithTagData
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

def GetPatternDetailsFromFileName(file_name_path: Path) -> str:
    file_relative_path = str(file_name_path.relative_to(Path(".")))
    if file_relative_path.endswith(".h5"):
        file_relative_path = file_relative_path[:-1]
        pattern, tag_type_dict = IdentifyPattern(file_relative_path)
        return pattern + "5", tag_type_dict
    else:
        return IdentifyPattern(file_relative_path)

def GetSortedFilesListFromPattern(path: Path,
                                  pattern: str,
                                  tag_type_dict: {"<time>": float, "<step>": int},
                                  sorting_functor=lambda _, *args: tuple(args)) -> 'list[tuple[typing.Any]]':
    return GetMachingEntitiesWithTagData(PathPatternEntity(path), pattern, tag_type_dict, sorting_functor)

def GetSpatialGridPaths(h5_file: h5py.File, list_of_model_data_paths: 'list[str]', current_path: str) -> None:
    current_data = h5_file[current_path]
    if isinstance(current_data, h5py.Group):
        if "__model_part_name" in current_data.attrs.keys():
            list_of_model_data_paths.append(current_path)
        else:
            for itr in current_data:
                GetSpatialGridPaths(h5_file, list_of_model_data_paths, f"{current_path}/{str(itr)}")

def __GetModelPartUniformGrids(h5_group: h5py.Group, grid_name: str, coordinates: Geometry) -> 'list[UniformGrid]':
    grids: 'list[UniformGrid]' = []
    for name, entity_group in h5_group.items():
        if isinstance(entity_group, h5py.Group):
            if "NumberOfNodes" in entity_group.attrs and "WorkingSpaceDimension" in entity_group.attrs and "Connectivities" in entity_group.keys():
                cell_type = TopologyCellType(entity_group.attrs["WorkingSpaceDimension"], entity_group.attrs["NumberOfNodes"])
                connectivities = HDF5UniformDataItem(entity_group["Connectivities"])
                topology = UniformMeshTopology(cell_type, connectivities)
                grids.append(UniformGrid(f"{grid_name}.{name}", coordinates, topology))
                KratosMultiphysics.Logger.PrintInfo("XDMF", f"Added {grid_name}.{name} spatial grid.")
            else:
                raise RuntimeError(f"Either \"Connectivities\", \"NumberOfNodes\" or \"WorkingSpaceDimension\" are not defined for {h5_group.name}.")
    return grids

def GenerateUniformGrids(h5_file: h5py.File, grid_path: str, grid_name: str, coordinates: Geometry, spatial_grid: SpatialGrid):
    grid_data = h5_file[grid_path]

    for name, h5_object in grid_data.items():
        if isinstance(h5_object, h5py.Dataset) and name == "Points":
            cell_type = TopologyCellType(3, 1)
            points = HDF5UniformDataItem(h5_object)
            topology = UniformMeshTopology(cell_type, points)
            spatial_grid.add_grid(UniformGrid(f"{grid_name}.Points", coordinates, topology))
            KratosMultiphysics.Logger.PrintInfo("XDMF", f"Added {grid_name}.Points spatial grid.")
        elif isinstance(h5_object, h5py.Group):
            if name in ["Conditions", "Elements"]:
                list(map(spatial_grid.add_grid, __GetModelPartUniformGrids(h5_object, grid_name, coordinates)))
            else:
                GenerateUniformGrids(h5_file, f"{grid_path}/{name}", f"{grid_name}.{name}", coordinates, spatial_grid)

def GenerateSpatialGrid(h5_file: h5py.File, model_data_path: str) -> 'tuple[SpatialGrid, list[UniformGrid], list[UniformGrid]]':
    spatial_grid = SpatialGrid()
    coordinates = Geometry(HDF5UniformDataItem(h5_file[f"{model_data_path}/Nodes/Local/Coordinates"]))

    model_data = h5_file[model_data_path]
    model_part_name = ''.join([chr(i) for i in h5_file[model_data_path].attrs["__model_part_name"]])

    if "Xdmf" in model_data.keys():
        xdmf_data = model_data["Xdmf"]

        # add root model part entities
        condition_grids: 'list[UniformGrid]' = []
        if "Conditions" in xdmf_data.keys():
            condition_grids = __GetModelPartUniformGrids(xdmf_data["Conditions"], model_part_name, coordinates)

        element_grids: 'list[UniformGrid]' = []
        if "Elements" in xdmf_data.keys():
            element_grids = __GetModelPartUniformGrids(xdmf_data["Elements"], model_part_name, coordinates)

        # now add sub model part entitites
        if "SubModelParts" in xdmf_data.keys():
            GenerateUniformGrids(h5_file, f"{model_data_path}/Xdmf/SubModelParts", model_part_name, coordinates, spatial_grid)
    else:
        KratosMultiphysics.Logger.PrintInfo("XDMF", f"No Xdmf data found in {h5_file.filename}:{model_data_path}.")

    return spatial_grid, condition_grids, element_grids

def GetAvailableDataSets(h5_file: h5py.File, group_paths_list: 'list[str]') -> 'list[EntityData]':
    # data_sets holds a dictionary of data sets found in the specified group_paths_list.
    # in data_sets, first we store the mesh_location, then container_type, then data_names list
    data_sets: 'list[EntityData]' = []

    def __check_item(_: str, h5_object: typing.Any) -> None:
        if isinstance(h5_object, h5py.Dataset):
            attribs = h5_object.attrs

            container_type = None
            if "__container_type" in attribs.keys():
                container_type = ''.join(chr(i) for i in attribs["__container_type"])
            mesh_location = None
            if "__mesh_location" in attribs.keys():
                mesh_location = ''.join(chr(i) for i in attribs["__mesh_location"])
            data_name = None
            if "__data_name" in attribs.keys():
                data_name = ''.join(chr(i) for i in attribs["__data_name"])

            if container_type and mesh_location and data_name:
                data_sets.append(EntityData(h5_object))

    for group_path in group_paths_list:
        h5_file[group_path].visititems(__check_item)

    return data_sets

class OrderedDataSets(abc.ABC):
    @abc.abstractmethod
    def Iterate(self) -> 'typing.Generator[float, str, list[EntityData]]':
        pass

class SingleMeshMultiFileSameOrderedDataSets(OrderedDataSets):
    def __init__(self,
                 path: Path,
                 file_name_pattern: str,
                 temporal_tag_position = 0,
                 tag_type_dict = {"<time>": float, "<step>": int},
                 sorting_functor=lambda _, *args: tuple(args)) -> None:
        # generated the sorted file data lists
        self.file_data_list =  GetSortedFilesListFromPattern(path, file_name_pattern, tag_type_dict, sorting_functor)
        self.temporal_tag_position = temporal_tag_position

        if not self.file_data_list:
            raise RuntimeError(f"No files found matching the file name pattern at {str(path)} [ file_name_pattern = {file_name_pattern} ].")

        # since all files should have the same results. We only open one file and get the data sets.
        with h5py.File(self.file_data_list[0][0].Get(), "r") as h5_file:
            self.data_sets = GetAvailableDataSets(h5_file, ["/"])

        if not self.data_sets:
            raise RuntimeError("No data sets found.")

    def Iterate(self) -> 'typing.Generator[float, str, list[EntityData]]':
        for file_data in self.file_data_list:
            yield float(file_data[self.temporal_tag_position + 1]), self.data_sets[0].mesh_location, self.data_sets

def WriteOrderedDataSetsToXdmf(ordered_datasets: OrderedDataSets, output_file_name: str) -> None:
    temporal_grid = TemporalGrid()

    spatial_grids: 'dict[str, tuple[SpatialGrid, list[UniformGrid], list[UniformGrid]]]' = {}

    for control_value, mesh_location, datasets in ordered_datasets.Iterate():
        # get the spatial grid
        if not mesh_location in spatial_grids.keys():
            mesh_data = mesh_location.split(":")
            with h5py.File(mesh_data[0], "r") as h5_file:
                spatial_grids[mesh_location] = GenerateSpatialGrid(h5_file, mesh_data[1])

        # now create a duplicate of the found spatial grid
        current_spatial_grid = SpatialGrid()
        spatial_grid, condition_grids, element_grids = spatial_grids[mesh_location]

        # add the condition and element grids
        current_condition_grids: 'list[UniformGrid]' = []
        for grid in condition_grids:
            condition_grid = UniformGrid(grid.name, grid.geometry, grid.topology)
            current_condition_grids.append(condition_grid)
            current_spatial_grid.add_grid(condition_grid)

        current_element_grids: 'list[UniformGrid]' = []
        for grid in element_grids:
            condition_grid = UniformGrid(grid.name, grid.geometry, grid.topology)
            current_element_grids.append(condition_grid)
            current_spatial_grid.add_grid(condition_grid)

        # now add submodel part grids
        for grid in spatial_grid.grids:
            current_spatial_grid.add_grid(UniformGrid(grid.name, grid.geometry, grid.topology))

        # now add the results
        for dataset in datasets:
            if dataset.container_type == "NODES":
                # add nodal results to all the grids
                for grid in current_spatial_grid.grids:
                    grid.add_attribute(dataset)
            else:
                if dataset.container_type == "CONDITIONS":
                    # now we have to check whether we have multiple condition types.
                    # because, the condition data is stored as one dataset, but for different conditons
                    # are stored in different groups with different xdmf indices.
                    if len(current_condition_grids) != 1:
                        raise RuntimeError(f"The condition data_set {dataset.name} can be written to only one condition data grid. [ number of condition grids = {len(current_condition_grids)} ].")
                    current_condition_grids[0].add_attribute(dataset)
                elif dataset.container_type == "ELEMENTS":
                    # now we have to check whether we have multiple condition types.
                    # because, the condition data is stored as one dataset, but for different conditons
                    # are stored in different groups with different xdmf indices.
                    if len(current_element_grids) != 1:
                        raise RuntimeError(f"The condition data_set {dataset.name} can be written to only one condition data grid. [ number of condition grids = {len(current_element_grids)} ].")
                    current_element_grids[0].add_attribute(dataset)
                else:
                    raise RuntimeError(f"Unsupported container_type = \"{dataset.container_type} ].")

        temporal_grid.add_grid(Time(control_value), current_spatial_grid)

    domain = Domain(temporal_grid)
    xdmf = Xdmf(domain)
    ET.ElementTree(xdmf.create_xml_element()).write(output_file_name)