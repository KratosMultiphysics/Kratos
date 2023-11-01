import re
import abc
import typing
from pathlib import Path

import KratosMultiphysics as Kratos

class PatternEntity(abc.ABC):
    @abc.abstractmethod
    def Iterate(self) -> 'typing.Generator[PatternEntity, None, None]':
        pass

    @abc.abstractmethod
    def IsLeaf(self) -> bool:
        pass

    @abc.abstractmethod
    def Name(self) -> str:
        pass

    @abc.abstractmethod
    def Get(self) -> typing.Any:
        pass

class PathPatternEntity(PatternEntity):
    def __init__(self, path: Path) -> None:
        self.__path = path

    def Iterate(self) -> 'typing.Generator[PathPatternEntity, None, None]':
        for itr in self.__path.iterdir():
            yield PathPatternEntity(itr)

    def IsLeaf(self) -> bool:
        return self.__path.is_file()

    def Name(self) -> str:
        return self.__path.name

    def Get(self) -> Path:
        return self.__path

class Pattern:
    def __init__(self, tagged_pattern: str, tag_type_dict: 'dict[str, typing.Any]') -> None:
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

        self.__converters: 'list[typing.Any]' = []

        tags: 'list[str]' = re.findall(r"[\w+]?(<\w+>)", tagged_pattern)
        for tag in tags:
            if tag in tag_type_dict.keys():
                self.__converters.append(tag_type_dict[tag])

        for tag, value_type in tag_type_dict.items():
            if value_type == int:
                regex_pattern = regex_pattern.replace(tag, r"([0-9]+)")
            elif value_type == float:
                regex_pattern = regex_pattern.replace(tag, r"([0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)")
            else:
                raise RuntimeError(f"Unsupported tag type for tag = \"{tag}\" [ type = {value_type} ].")

        self.__pattern = re.compile(regex_pattern)
        if len(tags) == 1:
            self.__value_getter_method = lambda x: x
        else:
            self.__value_getter_method = lambda x: x[0]

    def GetData(self, input_str: str) -> 'tuple[bool, list[typing.Any]]':
        values = self.__pattern.findall(input_str)
        if len(values) >= 1:
            converted_values = [converter(v) for v, converter in zip(self.__value_getter_method(values), self.__converters)]
            return True,  converted_values
        else:
            return False, []

    def GetNumberOfDataItems(self) -> int:
        return len(self.__converters)

def __GetMachingEntities(starting_entity: PatternEntity, patterns: 'list[str]', tag_type_dict: 'dict[str, typing.Any]', common_data: 'list[typing.Any]') -> 'typing.Generator[tuple[PatternEntity, ...], None, None]':
    if len(patterns) == 1:
        # this is the file pattern. hence we now try to find matching files
        pattern = Pattern(patterns[0], tag_type_dict)
        for itr in starting_entity.Iterate():
            is_valid, file_pattern_data = pattern.GetData(itr.Name())
            if is_valid:
                yield itr, *common_data, *file_pattern_data
    else:
        # this is a dir pattern
        pattern = Pattern(patterns[0], tag_type_dict)
        for itr in starting_entity.Iterate():
            if not itr.IsLeaf():
                is_valid, dir_pattern_data = pattern.GetData(itr.Name())
                if is_valid:
                    for data in __GetMachingEntities(itr, patterns[1:], tag_type_dict, [*common_data, *dir_pattern_data]):
                        yield data

def GetMachingEntities(starting_entity: PatternEntity, tagged_pattern: str, tag_type_dict: 'dict[str, typing.Any]') -> 'typing.Generator[tuple[PatternEntity, ...], None, None]':
    sub_patterns = tagged_pattern.split("/")
    for data in __GetMachingEntities(starting_entity, sub_patterns, tag_type_dict, []):
        yield data

def EvaluatePattern(pattern: str, model_part: Kratos.ModelPart, time_format='') -> str:
    time = model_part.ProcessInfo[Kratos.TIME]
    pattern = pattern.replace("<time>", format(time, time_format))

    if Kratos.STEP in model_part.ProcessInfo:
        step = model_part.ProcessInfo[Kratos.STEP]
    else:
        step = 0
    pattern = pattern.replace("<step>", str(step))
    pattern = pattern.replace("<model_part_name>", model_part.Name)
    pattern = pattern.replace("<model_part_full_name>", model_part.FullName())
    return pattern

def IdentifyPattern(entity_name: str) -> 'tuple[str, dict[str, typing.Any]]':
    # all the tag types are assumed to be of float type
    # here the "-" sign is omitted.
    float_pattern = re.compile(r"([0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)")
    values = float_pattern.findall(entity_name)

    current_pos = 0
    current_value_index = 0
    initial_value_index = 0
    pattern = ""
    tags_dict: 'dict[str, typing.Any]' = {}
    while current_pos < len(entity_name):
        found_tag = False
        current_value_index = initial_value_index
        while current_value_index < len(values):
            current_value = values[current_value_index]
            if entity_name[current_pos:current_pos + len(current_value)] == current_value:
                tags_dict[f"<float_{current_value_index+1}>"] = float
                pattern += f"<float_{current_value_index+1}>"
                current_pos += len(current_value)
                initial_value_index += 1
                found_tag = True
                break
            current_value_index += 1

        if not found_tag:
            pattern += entity_name[current_pos]
            current_pos += 1

    return pattern, tags_dict
