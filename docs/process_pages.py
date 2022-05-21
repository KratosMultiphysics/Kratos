from collections import OrderedDict
import json
from pathlib import Path

from entry_utilities import default_header_dict
from entry_utilities import GenerateEntryDataFromDir
from entry_utilities import GenerateEntryDataFromFile
from entry_utilities import GenerateEntryDataFromKratosExampleUrl
from entry_utilities import GetTaglessName
from entry_utilities import CreateNavigationBarEntry
from entry_utilities import GetDirMenuInfoFromJson
from entry_utilities import GenerateEntryDataFromExternalUrl

docs_absolute_path = Path(__file__)

spacing_information = {
    "dirs": [
        ["root", ""],
        ["folders", "  "],
        ["subfolders", "      "]
    ],
    "files": [
        ["unsupported", ""],
        ["folderitems", "    "],
        ["subfolderitems", "        "]
    ]
}

def GetEntryPathFromString(input_str: str, current_path: Path) -> Path:
    if input_str.startswith("/"):
        return docs_absolute_path / input_str[1:]
    else:
        return current_path / input_str

def CreateEntriesDicts(current_path: Path, navigation_level, max_navigation_level, default_headers_dict: dict) -> list:
    if navigation_level >= len(spacing_information["files"]):
        file_entry_type = "unsupported"
    else:
        file_entry_type = spacing_information["files"][navigation_level][0]

    if navigation_level >= len(spacing_information["dirs"]):
        dir_entry_type = "unsupported"
    else:
        dir_entry_type = spacing_information["dirs"][navigation_level][0]

    menu_data = GetDirMenuInfoFromJson(current_path)
    current_path_entry = GenerateEntryDataFromDir(current_path, dir_entry_type)

    if navigation_level == max_navigation_level:
        # if it is the leaf level, check and add landing pages for folders
        landing_url = ""
        if "landing_page" in menu_data.keys():
            landing_url = menu_data["landing_page"]
        else:
            raise RuntimeError("No landing page information found for leaf dir {0:s}. Please add it to {0:}/menu_info.json.".format(str(current_path)))
        landing_path = GetEntryPathFromString(landing_url, current_path)
        if landing_path.is_file():
            current_path_entry["path"] = landing_path
        else:
            raise RuntimeError("Landing page {:s} defined in {:s}/menu_info.json not found.".format(str(landing_path), str(current_path)))
        return [current_path_entry]
    else:
        list_of_entries = [current_path_entry]

    custom_entries_order = []
    if "custom_entries" in menu_data.keys():
        custom_entries_order = menu_data["custom_entries"]

    ignore_entries = []
    if "ignore_entries" in menu_data.keys():
        ignore_entries = menu_data["ignore_entries"]

    for i, custom_entry in enumerate(custom_entries_order):
        if custom_entry in ignore_entries:
            raise RuntimeError("Custom entry {:s} found in ignore entries list. These \"custom_entries\" and \"ignore_entries\" should be mutually exclusive.".format(custom_entry))

        if isinstance(custom_entry, dict):
            if "type" not in custom_entry.keys():
                raise RuntimeError("\"type\" key not found in {:s} at {:s}/menu_info.json".format(custom_entry, str(current_path)))
            if custom_entry["type"] == "kratos_example":
                current_dict = GenerateEntryDataFromKratosExampleUrl(current_path / custom_entry["file_name"], custom_entry["url"], file_entry_type, default_headers_dict)
                custom_entries_order[i] = GetTaglessName(current_dict["path"])
                list_of_entries.append(current_dict)
            elif custom_entry["type"] == "external":
                list_of_entries.append(GenerateEntryDataFromExternalUrl(custom_entry, file_entry_type))
            else:
                raise RuntimeError("Unsupported entry type \"{:s}\". \"kratos_example\" and \"external\" are the only supported entry types".format(custom_entry["type"]))
        elif isinstance(custom_entry, str):
            current_entry_path = GetEntryPathFromString(custom_entry, current_path)
            if current_entry_path.is_file():
                list_of_entries.append(GenerateEntryDataFromFile(current_entry_path, file_entry_type, default_headers_dict))
            elif current_entry_path.is_dir():
                list_of_entries.extend(CreateEntriesDicts(current_entry_path, navigation_level + 1, max_navigation_level, default_headers_dict))
            else:
                raise RuntimeError("Custom entry {:s} is not a file or a dir.".format(str(current_entry_path)))
        else:
            raise RuntimeError("Custom entry {:s} is not supported.".format(str(custom_entry)))

    temp_entries_list = []
    for iter_path in current_path.iterdir():
        if iter_path.is_file() and str(iter_path).endswith(".md"):
            file_name = GetTaglessName(iter_path)
            if file_name not in custom_entries_order and file_name not in ignore_entries:
                temp_entries_list.append(GenerateEntryDataFromFile(iter_path, file_entry_type, default_headers_dict))

    list_of_entries.extend(sorted(temp_entries_list, key=lambda x: x["title"]))

    temp_entries_list = []
    for iter_path in current_path.iterdir():
        if iter_path.is_dir():
            if not GetTaglessName(iter_path) in custom_entries_order:
                temp_entries_list.append(iter_path)

    temp_entries_list = sorted(temp_entries_list, key= lambda x: str(x))
    for iter_path in temp_entries_list:
        list_of_entries.extend(CreateEntriesDicts(iter_path, navigation_level + 1, max_navigation_level, default_headers_dict))

    return list_of_entries

def GetTypeInfo(item_type: str) -> str:
    for spacing_info in spacing_information["dirs"]:
        if item_type == spacing_info[0]:
            return spacing_info[1]
    for spacing_info in spacing_information["files"]:
        if item_type == spacing_info[0]:
            return spacing_info[1]

    raise Exception("Entry type {:s} is not defined.".format(item_type))

def GenerateStrings(list_of_dicts: list[dict]) -> list[str]:
    for itr_dict in list_of_dicts:
        itr_dict["str"] = CreateNavigationBarEntry(itr_dict)

    written_types = OrderedDict()
    written_types["root"]= False,
    written_types["folders"]= False,
    written_types["folderitems"]= False,
    written_types["subfolders"]= False,
    written_types["subfolderitems"]= False

    list_of_strings = []
    previous_dict_type = ""
    for dict_item in list_of_dicts:
        dict_type = dict_item["type"]
        if dict_type == "unsupported":
            print("Warning: Entry with unsupported level found at {:s}.".format(str(dict_item["path"])))
        else:
            if not written_types[dict_type]:
                written_types[dict_type] = True
                spacing_info = GetTypeInfo(dict_type)
                if spacing_info != "":
                    list_of_strings.append(spacing_info + dict_type + ":\n")
            elif previous_dict_type != dict_type:
                for k in reversed(written_types.keys()):
                    if k != dict_type:
                        written_types[k] = False
                    else:
                        break

        list_of_strings.append(dict_item["str"])
        previous_dict_type = dict_type

    return list_of_strings

if __name__ == "__main__":
    print("Creating top navigation bar...")
    # generate top navigation bar
    with open("_data/topnav.yml.orig", "r") as file_input:
        lines = file_input.readlines()

    list_of_entries = CreateEntriesDicts(
        Path("pages"),
        0,
        10,
        default_header_dict)

    list_of_entries = list_of_entries[1:]
    list_of_strings = GenerateStrings(list_of_entries)

    for entry in list_of_entries:
        print(str(entry["path"]))