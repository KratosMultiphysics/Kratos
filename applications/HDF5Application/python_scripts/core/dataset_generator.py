import abc
import typing
import h5py
import copy
from pathlib import Path
from operator import itemgetter

from KratosMultiphysics.HDF5Application.core.xdmf import EntityData
from KratosMultiphysics.HDF5Application.core.pattern import PatternEntity
from KratosMultiphysics.HDF5Application.core.pattern import PathPatternEntity
from KratosMultiphysics.HDF5Application.core.pattern import GetMachingEntities

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

class HDF5PatternEntity(PatternEntity):
    def __init__(self, name: str, hdf5_object: 'typing.Union[h5py.Dataset, h5py.Group]') -> None:
        self.__hdf5_object = hdf5_object
        self.__name = name

    def Iterate(self) -> 'typing.Generator[HDF5PatternEntity, None, None]':
        for name, itr in self.__hdf5_object.items():
            yield HDF5PatternEntity(name, itr)

    def IsLeaf(self) -> bool:
        return isinstance(self.__hdf5_object, h5py.Dataset)

    def Name(self) -> str:
        return self.__name

    def Get(self) -> 'typing.Union[h5py.Dataset, h5py.Group]':
        return self.__hdf5_object

def GetDataSetPatterns(dataset_pattern: str) -> 'tuple[str, str]':
    if dataset_pattern == "":
        raise RuntimeError(f"Dataset pattern cannot be empty.")

    dataset_pattern_data = dataset_pattern.split(":")

    if len(dataset_pattern_data) == 2:
        if not dataset_pattern_data[1].startswith("/"):
            raise RuntimeError(f"Dataset prefix pattern must always start with \"/\" [ dataset_pattern = {dataset_pattern} ].")
        return dataset_pattern_data[0], dataset_pattern_data[1]
    elif len(dataset_pattern_data) == 1:
        return dataset_pattern_data[0], "/"
    else:
        raise RuntimeError(f"Dataset pattern should only have one \":\" delimeter seperating fileanme and dataset prefix [ dataset_pattern = \"{dataset_pattern}\" ].")

def HasTags(pattern: str, tag_type_dict: 'dict[str, typing.Any]') -> bool:
    for tag in tag_type_dict.keys():
        if pattern.find(tag) != -1:
            return True

    return False

class DatasetGenerator(abc.ABC):
    @abc.abstractmethod
    def Iterate(self) -> 'typing.Generator[tuple[typing.Union[int, float], EntityData], None, None]':
        pass

class SingleMeshMultiFileSameDatasetGenerator(DatasetGenerator):
    def __init__(self, dataset_pattern: str, temporal_value_tag_position: int = 0, tag_type_dict: 'dict[str, typing.Any]' = {"<time>": float, "<step>": int}) -> None:
        self.hdf5_file_name_pattern, self.dataset_prefix = GetDataSetPatterns(dataset_pattern)

        if HasTags(self.dataset_prefix, tag_type_dict):
            raise RuntimeError(f"Dataset prefix tags are not supported in SingleMeshMultiFileSameDatasetGenerator [ dataset_pattern: {dataset_pattern} ].")

        if not HasTags(self.hdf5_file_name_pattern, tag_type_dict):
            raise RuntimeError(f"File name is required to have tags in SingleMeshMultiFileSameDatasetGenerator [ dataset_pattern = {dataset_pattern} ].")

        if temporal_value_tag_position < 0:
            raise RuntimeError(f"Temporal value tag position should be positive [ temporal_value_tag_position = {temporal_value_tag_position} ].")

        self.temporal_value_tag_position = temporal_value_tag_position
        self.tag_type_dict = tag_type_dict

    def Iterate(self) -> 'typing.Generator[tuple[typing.Union[int, float], EntityData], None, None]':
        # get the generator
        generator: 'typing.Generator[tuple[PatternEntity, ...], None, None]' = GetMachingEntities(PathPatternEntity(Path(".")), self.hdf5_file_name_pattern, self.tag_type_dict)

        # get the first file
        first_file_data = next(generator, None)
        if first_file_data is None:
            raise RuntimeError(f"No matching dataset files found at \"{Path('.').absolute()}\" for file_name_pattern = \"{self.hdf5_file_name_pattern}\".")

        with h5py.File(str(first_file_data[0].Get())) as h5_file:
            datasets = GetAvailableDataSets(h5_file, [self.dataset_prefix])

        # yield the first finding
        for dataset in datasets:
            yield itemgetter(self.temporal_value_tag_position + 1)(first_file_data), dataset

        # now yield the rest of the files
        for file_data in generator:
            temporal_value = itemgetter(self.temporal_value_tag_position + 1)(file_data)
            for dataset in datasets:
                file_path: Path = file_data[0].Get()
                if file_path.is_file():
                    file_name = str(file_path.relative_to(Path(".")))
                    copied_dataset = copy.deepcopy(dataset)
                    copied_dataset.data.file_name = file_name
                    yield temporal_value, copied_dataset

class GenericDatasetGenerator(DatasetGenerator):
    def __init__(self, dataset_pattern: str, temporal_value_tag_position: int = 0, tag_type_dict: 'dict[str, typing.Any]' = {"<time>": float, "<step>": int}) -> None:
        self.hdf5_file_name_pattern, self.dataset_prefix_pattern = GetDataSetPatterns(dataset_pattern)

        self.temporal_value_tag_position = temporal_value_tag_position
        self.tag_type_dict = tag_type_dict

        if not self.dataset_prefix_pattern.startswith("/"):
            raise RuntimeError(f"The dataset prefix pattern should always start with \"\\\" [ dataset_pattern = {dataset_pattern} ].")

        # remove the first "/" character for pattern mathing.
        self.dataset_prefix_pattern = self.dataset_prefix_pattern[1:]

    def Iterate(self) -> 'typing.Generator[tuple[typing.Union[int, float], EntityData], None, None]':
        # get the generator
        file_generator: 'typing.Generator[tuple[PatternEntity, ...], None, None]' = GetMachingEntities(PathPatternEntity(Path(".")), self.hdf5_file_name_pattern, self.tag_type_dict)

        try:
            first_file = next(file_generator)
            if len(first_file) - 1 > self.temporal_value_tag_position:
                temporal_value_getter = lambda h5_file_data, _: itemgetter(self.temporal_value_tag_position + 1)(h5_file_data)
            else:
                temporal_value_getter = lambda _, dataset_data: itemgetter(self.temporal_value_tag_position + 2 - len(first_file))(dataset_data)

            # now yield the data of the first file
            for data in self.__GetFileDataSets(first_file, temporal_value_getter):
                yield data

            # now yeild data of the rest of the files
            for file_data in file_generator:
                for data in self.__GetFileDataSets(file_data, temporal_value_getter):
                    yield data
        except StopIteration:
            raise RuntimeError(f"No matching files were found for the pattern \"{self.hdf5_file_name_pattern}:{self.dataset_prefix_pattern}\".")

    def __GetFileDataSets(self, file_data: 'tuple[PatternEntity, ...]', temporal_value_getter: typing.Any) -> 'typing.Generator[tuple[float, EntityData], None, None]':
        with h5py.File(file_data[0].Get()) as h5_file:
            dataset_generator: 'typing.Generator[tuple[PatternEntity, ...], None, None]' = GetMachingEntities(HDF5PatternEntity("", h5_file), self.dataset_prefix_pattern, self.tag_type_dict)
            for dataset_data in dataset_generator:
                temporal_value = temporal_value_getter(file_data, dataset_data)
                h5_object: 'typing.Union[h5py.Group, h5py.Dataset]' = dataset_data[0].Get()
                if isinstance(h5_object, h5py.Group):
                    for sub_dataset in GetAvailableDataSets(h5_file, [h5_object.name]):
                        yield temporal_value, sub_dataset
                else:
                    yield temporal_value, EntityData(h5_object)

