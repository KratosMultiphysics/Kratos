import abc
import typing
import h5py
import copy
from pathlib import Path

from .xdmf import EntityData
from .pattern import PatternEntity
from .pattern import PathPatternEntity
from .pattern import GetMachingEntities

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

    def Iterate(self) -> 'typing.Generator[HDF5PatternEntity]':
        for name, itr in self.__hdf5_object.items():
            yield HDF5PatternEntity(name, itr)

    def IsLeaf(self) -> bool:
        return isinstance(self.__hdf5_object, h5py.Dataset)

    def Name(self) -> str:
        return self.__name

    def Get(self) -> 'typing.Union[h5py.Dataset, h5py.Group]':
        return self.__hdf5_object

class DataSetGenerator(abc.ABC):
    @abc.abstractmethod
    def Iterate(self) -> 'typing.Generator[tuple[float, EntityData]]':
        pass

    @staticmethod
    def _GetDataSetPatterns(dataset_pattern: str) -> 'tuple[str, str]':
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

    @staticmethod
    def _HasTags(pattern: str, tag_type_dict: 'dict[str, typing.Any]') -> bool:
        for tag in tag_type_dict.keys():
            if pattern.find(tag) != -1:
                return True

        return False

class SingleMeshMultiFileSameDatasetsGenerator(DataSetGenerator):
    def __init__(self, dataset_pattern: str, temporal_value_tag_position: int = 0, tag_type_dict: 'dict[str, typing.Any]' = {"<time>": float, "<step>": int}) -> None:
        self.hdf5_file_name_pattern, self.dataset_prefix = self._GetDataSetPatterns(dataset_pattern)

        if self._HasTags(self.dataset_prefix, tag_type_dict):
            raise RuntimeError(f"Dataset prefix tags are not supported in SingleFileDatasetsGenerator [ dataset_pattern: {dataset_pattern} ].")

        self.temporal_value_tag_position = temporal_value_tag_position
        self.tag_type_dict = tag_type_dict

    def Iterate(self) -> 'typing.Generator[tuple[float, EntityData]]':
        # get the generator
        generator: 'typing.Generator[tuple[PathPatternEntity, typing.Any]]' = GetMachingEntities(PathPatternEntity(Path(".")), self.hdf5_file_name_pattern, self.tag_type_dict)

        # get the first file
        first_file_data = next(generator, None)
        if first_file_data is None:
            raise RuntimeError(f"No matching dataset files found at \"{Path('.').absolute()}\" for file_name_pattern = \"{self.hdf5_file_name_pattern}\".")

        with h5py.File(str(first_file_data[0].Get())) as h5_file:
            datasets = GetAvailableDataSets(h5_file, [self.dataset_prefix])

        # yield the first finding
        for dataset in datasets:
            yield first_file_data[self.temporal_value_tag_position + 1], dataset

        # now yield the rest of the files
        for file_data in generator:
            control_value = file_data[self.temporal_value_tag_position + 1]
            for dataset in datasets:
                file_name = str(file_data[0].Get().relative_to(Path(".")))
                copied_dataset = copy.deepcopy(dataset)
                copied_dataset.data.file_name = file_name
                yield control_value, copied_dataset

class SingleFileDatasetsGenerator(DataSetGenerator):
    def __init__(self, dataset_pattern: str, temporal_value_tag_position: int = 0, tag_type_dict: 'dict[str, typing.Any]' = {"<time>": float, "<step>": int}) -> None:
        self.hdf5_file_name, self.dataset_prefix_pattern = self._GetDataSetPatterns(dataset_pattern)

        if self._HasTags(self.hdf5_file_name, tag_type_dict):
            raise RuntimeError(f"File name tags are not supported in SingleFileDatasetsGenerator [ dataset_pattern: {dataset_pattern} ].")

        self.temporal_value_tag_position = temporal_value_tag_position
        self.tag_type_dict = tag_type_dict

        # remove the first "/" character for pattern mathing.
        self.dataset_prefix_pattern = self.dataset_prefix_pattern[1:]

    def Iterate(self) -> 'typing.Generator[tuple[float, EntityData]]':
        with h5py.File(self.hdf5_file_name) as h5_file:
            # get the generator
            generator: 'typing.Generator[tuple[HDF5PatternEntity, typing.Any]]' = GetMachingEntities(HDF5PatternEntity("", h5_file), self.dataset_prefix_pattern, self.tag_type_dict)

            for data in generator:
                yield data[self.temporal_value_tag_position + 1], EntityData(data[0].Get())

class MultiFileDatasetsGenerator(DataSetGenerator):
    def __init__(self, dataset_pattern: str, temporal_value_tag_position: int = 0, tag_type_dict: 'dict[str, typing.Any]' = {"<time>": float, "<step>": int}) -> None:
        self.hdf5_file_name_pattern, self.dataset_prefix = self._GetDataSetPatterns(dataset_pattern)

        if self._HasTags(self.dataset_prefix, tag_type_dict):
            raise RuntimeError(f"Dataset prefix tags are not supported in SingleFileDatasetsGenerator [ dataset_pattern: {dataset_pattern} ].")

        self.temporal_value_tag_position = temporal_value_tag_position
        self.tag_type_dict = tag_type_dict

    def Iterate(self) -> 'typing.Generator[tuple[float, EntityData]]':
        # get the generator
        generator: 'typing.Generator[tuple[PathPatternEntity, typing.Any]]' = GetMachingEntities(PathPatternEntity(Path(".")), self.hdf5_file_name_pattern, self.tag_type_dict)

        for data in generator:
            control_value = data[self.temporal_value_tag_position + 1]
            with h5py.File(data[0].Get()) as h5_file:
                datasets = GetAvailableDataSets(h5_file, [self.dataset_prefix])
                for dataset in datasets:
                    yield control_value, dataset