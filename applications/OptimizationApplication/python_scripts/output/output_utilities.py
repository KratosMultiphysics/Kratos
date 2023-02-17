import re

import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.optimization_info import OptimizationInfo

def GetOptimizationInfoAvailableKeysForType(optimization_info: OptimizationInfo, required_type = None) -> 'list[str]':
    def __get_valid_keys(valid_keys: 'list[str]', current_prefix: str, d: dict):
        for k, v in d.items():
            if isinstance(v, dict):
                __get_valid_keys(valid_keys, f"{current_prefix}/{k}", d[k])
            elif required_type is None or isinstance(v, required_type):
                valid_keys.append(f"{current_prefix}/{k}")

    valid_keys = []
    __get_valid_keys(valid_keys, "", optimization_info.GetSolutionStepData(0))
    valid_keys = [valid_key[1:] for valid_key in valid_keys]
    return valid_keys

def GetPlaceHolders(pattern: str) -> 'list[str]':
    return ["<" + place_holder_start[:place_holder_start.find(">")+1] for place_holder_start in pattern.split("<")][1:]

def GetPatternMatchedOptimizationKeys(available_keys: 'list[str]', pattern: str, is_with_data_name: bool) -> 'list[str]':
    # first get the list of place holders
    # this matches any keyword starting with < and ending with > and keep them as place holders
    place_holders = GetPlaceHolders(pattern)

    # now build the regex pattern
    regex_pattern = pattern
    for place_holder in place_holders:
        regex_pattern = regex_pattern.replace(f"{place_holder}", r"(.*?[^\/]*)")

    # now add the last position pattern for data names
    if is_with_data_name:
        regex_pattern = regex_pattern + r"/(.*?$)"
        place_holders.append("<data_name>")

    # replace problamatic "/"
    regex_pattern = regex_pattern.replace("/", "\/")

    # get the matching list
    groups_list = [re.match(regex_pattern, key) for key in available_keys]
    matched_groups_list = [matched_group for matched_group in groups_list if not matched_group is None]

    results: 'list[str]' = []
    for matched_groups in matched_groups_list:
        current_result = ", ".join([f"{place_holder[1:-1]} = {value}" for place_holder, value in zip(place_holders, matched_groups.groups())])
        results.append(current_result)

    return sorted(list(set(results)))

def GetPrefixWithoutDataName(pattern: str, settings: Kratos.Parameters) -> str:
    place_holders = GetPlaceHolders(pattern)
    common_prefix = pattern
    for place_holder in place_holders:
        common_prefix = common_prefix.replace(place_holder, settings[place_holder[1:-1]].GetString())
    return common_prefix

def GetPlaceHolderWithValues(pattern: str, settings: Kratos.Parameters) -> str:
    place_holders = GetPlaceHolders(pattern)
    return ", ".join(f'"{place_holder[1:-1]}" = "{settings[place_holder[1:-1]].GetString()}"' for place_holder in place_holders)
