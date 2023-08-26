'''HDF5 core utils.

license: HDF5Application/license.txt
'''


__all__ = ['ParametersWrapper']

import re
from typing import Any
from operator import itemgetter
from pathlib import Path

from collections.abc import Mapping

import KratosMultiphysics

class Pattern:
    def __init__(self, tagged_pattern: str, tag_type_dict: 'dict[str, Any]') -> None:
        regex_special_chars_escape_map = {
            "+": r"\+",
            "^": r"\^",
            "$": r"\$",
            ".": r"\.",
            "|": r"\|",
            "?": r"\?",
            "*": r"\*",
            "(": r"\(",
            ")": r"\)",
            "[": r"\[",
            "]": r"\]",
            "{": r"\{",
            "}": r"\}"
        }

        # replace regex chars
        regex_pattern = tagged_pattern.replace("\\", "\\\\")
        for k, v in regex_special_chars_escape_map.items():
            regex_pattern = regex_pattern.replace(k, v)

        self.__converters: 'list[Any]' = []

        tags: 'list[str]' = re.findall(r"[\w+]?(<\w+>)", tagged_pattern)
        for tag in tags:
            if tag in tag_type_dict.keys():
                self.__converters.append(tag_type_dict[tag])

        for tag, value_type in tag_type_dict.items():
            if value_type == int:
                regex_pattern = regex_pattern.replace(tag, r"([-+]?[0-9]+)")
            elif value_type == float:
                regex_pattern = regex_pattern.replace(tag, r"([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)")
            else:
                raise RuntimeError(f"Unsupported tag type for tag = \"{tag}\" [ type = {value_type} ].")

        self.__pattern = re.compile(regex_pattern)
        if len(tags) == 1:
            self.__value_getter_method = lambda x: x
        else:
            self.__value_getter_method = lambda x: x[0]

    def GetData(self, input_str: str) -> 'tuple[bool, list[Any]]':
        values = self.__pattern.findall(input_str)
        if len(values) >= 1:
            converted_values = [converter(v) for v, converter in zip(self.__value_getter_method(values), self.__converters)]
            return True,  converted_values
        else:
            return False, []

    def GetNumberOfDataItems(self) -> int:
        return len(self.__converters)

def GetFilesList(path: Path, patterns: 'list[str]', tag_type_dict: 'dict[str, Any]', common_data: 'list[Any]') -> 'list[Path]':
    data = []
    if len(patterns) == 1:
        # this is the file pattern. hence we now try to find matching files
        pattern = Pattern(patterns[0], tag_type_dict)
        for itr_dir in path.iterdir():
            if itr_dir.is_file():
                is_valid, file_pattern_data = pattern.GetData(itr_dir.name)
                if is_valid:
                    data.append([itr_dir, *common_data, *file_pattern_data])
        return data
    else:
        # this is a dir pattern
        pattern = Pattern(patterns[0], tag_type_dict)
        for itr_dir in path.iterdir():
            if itr_dir.is_dir():
                is_valid, dir_pattern_data = pattern.GetData(itr_dir.name)
                if is_valid:
                    data.extend(GetFilesList(itr_dir, patterns[1:], tag_type_dict, [*common_data, *dir_pattern_data]))
        return data

def GetPatternSortedFilesList(path: Path, tagged_pattern: str, tag_type_dict: 'dict[str, Any]') -> 'list[str]':
    sub_patterns = tagged_pattern.split("/")
    file_paths_list = GetFilesList(path, sub_patterns, tag_type_dict, [])

    pattern = Pattern(tagged_pattern, tag_type_dict)

    file_paths_list = sorted(file_paths_list, key=itemgetter(*[v + 1 for v in range(pattern.GetNumberOfDataItems())]))
    return [v[0] for v in file_paths_list]

def EvaluatePattern(pattern: str, model_part: KratosMultiphysics.ModelPart, time_format='') -> str:
    time = model_part.ProcessInfo[KratosMultiphysics.TIME]
    pattern = pattern.replace("<time>", format(time, time_format))

    if KratosMultiphysics.STEP in model_part.ProcessInfo:
        step = model_part.ProcessInfo[KratosMultiphysics.STEP]
    else:
        step = 0
    pattern = pattern.replace("<step>", str(step))
    pattern = pattern.replace("<model_part_name>", model_part.Name)
    pattern = pattern.replace("<model_part_full_name>", model_part.FullName())
    return pattern

class ParametersWrapper(Mapping):
    '''A pythonic wrapper to KratosMultiphysics.Parameters.

    The idea is to reduce boilerplate and improve call-site readability by
    making Parameters more pythonic without breaking existing code.
    '''

    def __init__(self, params="{}"):
        if isinstance(params, self.__class__):
            self._parameters = params._parameters
        elif isinstance(params, str):
            self._parameters = KratosMultiphysics.Parameters(params)
        else:
            self._parameters = params

    def __getitem__(self, key):
        try:
            param = self._parameters[key]
        except:
            raise KeyError
        if param.IsString():
            value = param.GetString()
        elif param.IsDouble():
            value = param.GetDouble()
        elif param.IsBool():
            value = param.GetBool()
        elif param.IsInt():
            value = param.GetInt()
        else:
            value = self.__class__(param)
        return value

    def _convert_list_to_parameters(self, list_):
        dummy_params = KratosMultiphysics.Parameters()
        dummy_params.AddEmptyList('list')
        params = dummy_params['list']
        for v in list_:
            if isinstance(v, list):
                params.Append(self._convert_list_to_parameters(v))
            else:
                params.Append(self._convert_value_to_parameters(v))
        return params

    def _convert_value_to_parameters(self, value):
        params = KratosMultiphysics.Parameters()
        if isinstance(value, str):
            params.SetString(value)
        elif isinstance(value, float):
            params.SetDouble(value)
        elif isinstance(value, bool):
            params.SetBool(value)
        elif isinstance(value, int):
            params.SetInt(value)
        elif isinstance(value, KratosMultiphysics.Parameters):
            params = value
        elif isinstance(value, self.__class__):
            params = value._parameters
        else:
            raise TypeError()
        return params

    def __setitem__(self, key, value):
        if isinstance(value, list):
            params_object = self._convert_list_to_parameters(value)
        else:
            params_object = self._convert_value_to_parameters(value)
        if self._parameters.IsArray() or self._parameters.Has(key):
            self._parameters[key] = params_object
        else:
            self._parameters.AddValue(key, params_object)

    def __iter__(self):
        '''Return an iterator over keys or indices.

        This iterates over keys or indices, depending on if the Parameters
        instance is an array.
        '''
        if self._parameters.IsArray():
            yield from range(len(self))
        else:
            yield from self._parameters.keys()

    def __len__(self):
        return self._parameters.size()

    def __str__(self):
        return self._parameters.PrettyPrintJsonString()

    def __getattr__(self, name):
        '''Expose methods defined in KratosMultiphysics.Parameters.'''
        return getattr(self._parameters, name)

    def Get(self):
        return self._parameters

    def SetDefault(self, key, value=KratosMultiphysics.Parameters()):
        if key not in self:
            self[key] = value
