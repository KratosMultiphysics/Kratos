import requests
from argparse import ArgumentParser
from collections import OrderedDict
from pathlib import Path
import subprocess

from entry_utilities import default_header_dict
from entry_utilities import file_navigation_levels
from entry_utilities import spacing_information
from entry_utilities import GenerateEntryDataFromDir
from entry_utilities import GenerateEntryDataFromFile
from entry_utilities import GenerateEntryDataFromKratosExampleUrl
from entry_utilities import GetTaglessName
from entry_utilities import AddMissingUrlEntry
from entry_utilities import CreateNavigationBarEntry
from entry_utilities import GetDirMenuInfoFromJson
from entry_utilities import GenerateEntryDataFromExternalUrl
from entry_utilities import IsLeafEntry
from entry_utilities import GetNavigationString

web_prefix = "Kratos/"

docs_absolute_path = Path(__file__)

def GetEntryPathFromString(input_str: str, current_path: Path) -> Path:
    if input_str.startswith("/"):
        return docs_absolute_path / input_str[1:]
    else:
        return current_path / input_str

def CreateEntriesDicts(current_path: Path, navigation_level: int, max_navigation_level: int, default_headers_dict: dict) -> list:
    menu_data = GetDirMenuInfoFromJson(current_path)
    current_path_entry = GenerateEntryDataFromDir(current_path, navigation_level)

    if navigation_level == max_navigation_level:
        # if it is the leaf level, check and add landing pages for folders
        landing_url = ""
        if "landing_page" in menu_data.keys():
            landing_url = menu_data["landing_page"]
        else:
            # check whether this dir has markdown files. If so throw an error stating it has reached maximum levels.
            for iter_dir in current_path.iterdir():
                if iter_dir.is_file() and str(iter_dir).endswith(".md"):
                    raise Exception("Found pages in {:s} which is above the maximum supported levels in the navigation bar. Please move them to one of the parent folders.".format(str(current_path)))
            return []
        landing_path = GetEntryPathFromString(landing_url, current_path)
        if landing_path.is_file():
            current_path_entry["path"] = landing_path
            current_path_entry["type"] = GetNavigationString(file_navigation_levels, navigation_level - 1)
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

    def check_and_add_sub_dirs(dir_path: Path):
        list_of_sub_entries = CreateEntriesDicts(dir_path, navigation_level + 1, max_navigation_level, default_headers_dict)

        # check whether the subdirectory has pages. If not do not add it to menu
        if len(list_of_sub_entries) == 1 and not IsLeafEntry(list_of_sub_entries[0]):
            return []

        found_sub_leaf_entry = False
        found_sub_entries = False

        for sub_entry in list_of_sub_entries:
            if sub_entry["navigation_level"] == navigation_level + 1:
                found_sub_entries = True
                if IsLeafEntry(sub_entry):
                    found_sub_leaf_entry = True
                    break
        if found_sub_entries and not found_sub_leaf_entry:
            # add a dummy entry
            dummy_dict = {
                "title": "",
                "output": "web",
                "path": Path("."),
                "url": "<dummy>",
                "type": GetNavigationString(file_navigation_levels, navigation_level + 1),
                "navigation_level": navigation_level + 1
            }
            list_of_sub_entries.insert(1, dummy_dict)

        return list_of_sub_entries

    for i, custom_entry in enumerate(custom_entries_order):
        if custom_entry in ignore_entries:
            raise RuntimeError("Custom entry {:s} found in ignore entries list. These \"custom_entries\" and \"ignore_entries\" should be mutually exclusive.".format(custom_entry))

        if isinstance(custom_entry, dict):
            if "type" not in custom_entry.keys():
                raise RuntimeError("\"type\" key not found in {:s} at {:s}/menu_info.json".format(custom_entry, str(current_path)))
            if custom_entry["type"] == "kratos_example":
                current_dict = GenerateEntryDataFromKratosExampleUrl(current_path / custom_entry["file_name"], custom_entry, navigation_level, default_headers_dict)
                custom_entries_order[i] = GetTaglessName(current_dict["path"])
                list_of_entries.append(current_dict)
            elif custom_entry["type"] == "external":
                list_of_entries.append(GenerateEntryDataFromExternalUrl(custom_entry, navigation_level))
            else:
                raise RuntimeError("Unsupported entry type \"{:s}\". \"kratos_example\" and \"external\" are the only supported entry types".format(custom_entry["type"]))
        elif isinstance(custom_entry, str):
            current_entry_path = GetEntryPathFromString(custom_entry, current_path)
            if current_entry_path.is_file():
                list_of_entries.append(GenerateEntryDataFromFile(current_entry_path, navigation_level, default_headers_dict))
            elif current_entry_path.is_dir():
                list_of_entries.extend(check_and_add_sub_dirs(current_entry_path))
            else:
                raise RuntimeError("Custom entry {:s} is not a file or a dir.".format(str(current_entry_path)))
        else:
            raise RuntimeError("Custom entry {:s} is not supported.".format(str(custom_entry)))

    temp_entries_list = []
    for iter_path in current_path.iterdir():
        if iter_path.is_file() and str(iter_path).endswith(".md"):
            file_name = GetTaglessName(iter_path)
            if file_name not in custom_entries_order and file_name not in ignore_entries:
                temp_entries_list.append(GenerateEntryDataFromFile(iter_path, navigation_level, default_headers_dict))

    list_of_entries.extend(sorted(temp_entries_list, key=lambda x: x["title"]))

    temp_entries_list = []
    for iter_path in current_path.iterdir():
        if iter_path.is_dir():
            if not GetTaglessName(iter_path) in custom_entries_order:
                temp_entries_list.append(iter_path)

    temp_entries_list = sorted(temp_entries_list, key= lambda x: str(x))
    for iter_path in temp_entries_list:
        list_of_entries.extend(check_and_add_sub_dirs(iter_path))

    return list_of_entries

def GenerateStrings(list_of_dicts: list[dict], is_locally_built) -> list[str]:
    for itr_dict in list_of_dicts:
        if is_locally_built:
            AddMissingUrlEntry(itr_dict, "")
        else:
            AddMissingUrlEntry(itr_dict, web_prefix)
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
                spacing_info = spacing_information[dict_type]
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

def UpdateRootEntry(list_of_entries: list) -> list:
    root_entry = list_of_entries[0]
    menu_info = GetDirMenuInfoFromJson(root_entry["path"])
    if "additional_menu_options" in menu_info.keys():
        for k, v in menu_info["additional_menu_options"].items():
            root_entry[k] = v

    list_of_entries[0] = root_entry
    return list_of_entries

def ClearDummies(list_of_entries: list) -> list:
    index = 0
    while index < len(list_of_entries):
        current_entry = list_of_entries[index]
        remove_dummy = False
        if "url" in current_entry.keys() and current_entry["url"] == "<dummy>":
            if index + 1 < len(list_of_entries):
                next_entry = list_of_entries[index + 1]
                if IsLeafEntry(next_entry) and next_entry["type"] == current_entry["type"]:
                    remove_dummy = True

        if remove_dummy:
            del list_of_entries[index]
        else:
            index += 1

    return list_of_entries

def CreateNavigatonBar(root_path: str, max_levels: int, default_header_dict: dict, is_locally_built) -> list:
    list_of_entries = CreateEntriesDicts(
        Path(root_path),
        0,
        max_levels,
        default_header_dict)
    list_of_entries = UpdateRootEntry(list_of_entries)
    list_of_entries = ClearDummies(list_of_entries)
    list_of_strings = GenerateStrings(list_of_entries, is_locally_built)
    return list_of_strings

def AddPythonSnippetOutputs(file_path: Path) -> None:
    with open(file_path, "r") as file_input:
        lines = file_input.readlines()

    found_python_snippet_block = False
    snippet_lines = []
    output_lines = []
    index = 0
    while index < len(lines):
        line = lines[index]
        index += 1
        if not found_python_snippet_block or line != "## POST_PROCESS_PAGES_PYTHON_OUTPUT_GENERATION\n":
            output_lines.append(line)

        if found_python_snippet_block and line.startswith("```"):
            found_python_snippet_block = False
            if "## POST_PROCESS_PAGES_PYTHON_OUTPUT_GENERATION\n" in snippet_lines:
                # check whether existing expected output is found
                temp_index = index
                is_existing_output_found = False
                while temp_index < len(lines):
                    if lines[temp_index].strip():
                        is_existing_output_found  = lines[temp_index] == "Expected output:\n"
                        is_existing_output_found &= lines[temp_index+1] == "```console\n"
                        break
                    temp_index += 1

                if not is_existing_output_found:
                    # create a temp file
                    temp_file_path = f"{file_path.name}.temp.py"
                    with open(temp_file_path, "w") as temp_file_output:
                        temp_file_output.writelines(snippet_lines)

                    subprocess_run = subprocess.run([GetPython3Command(), "-u", temp_file_path], stdout=subprocess.PIPE, universal_newlines=True, check=True)
                    output_lines.append("\n")
                    output_lines.append("Expected output:\n")
                    output_lines.append("```console\n")
                    output_lines.append(subprocess_run.stdout)
                    output_lines.append("```\n")

                    # remove the temp file
                    Path(temp_file_path).unlink()

                    if is_existing_output_found:
                        index += temp_index + lines[temp_index:].index("```\n") + 1

        if found_python_snippet_block:
            snippet_lines.append(line)

        if line == "```python\n":
            found_python_snippet_block = True
            snippet_lines = []

    with open(file_path, "w") as file_output:
        file_output.writelines(output_lines)

if __name__ == "__main__":
    parser = ArgumentParser(description="Process mark down files in pages folder to create navigation bars.")
    parser.add_argument("-t", "--build_type", dest="build_type", metavar="<build_type>",
                        choices=['local', 'web'], default="locally", help="type of the web page build")
    args = parser.parse_args()

    is_locally_built = args.build_type == "local"

    Path("_data").mkdir(parents=True, exist_ok=True)
    print("Creating top navigation bar...")
    # generate top navigation bar
    with open("pages/topnav.yml.orig", "r") as file_input:
        lines = file_input.readlines()
    list_of_strings = CreateNavigatonBar("pages", 2, default_header_dict, is_locally_built)
    lines.extend(list_of_strings)
    with open("_data/topnav.yml", "w") as file_output:
        file_output.writelines(lines)

    Path("_data/sidebars").mkdir(parents=True, exist_ok=True)
    # generate side bar
    for iter_dir in Path("pages").iterdir():
        if iter_dir.is_dir():
            for sub_itr_dir in iter_dir.iterdir():
                if sub_itr_dir.is_dir():
                    print("Creating side bar for {:s}...".format(str(sub_itr_dir)))
                    menu_info = GetDirMenuInfoFromJson(sub_itr_dir)
                    if "side_bar_name" not in menu_info.keys():
                        raise RuntimeError("No side bar name provied. Please add it to {:s}/menu_info.json.".format(str(sub_itr_dir)))
                    default_header_dict["sidebar"] = "<!>" + menu_info["side_bar_name"]
                    list_of_strings = CreateNavigatonBar(str(sub_itr_dir), 3, default_header_dict, True)
                    with open("_data/sidebars/{:s}.yml".format(menu_info["side_bar_name"]), "w") as file_output:
                        file_output.write("entries:\n")
                        file_output.writelines(list_of_strings)

    if is_locally_built:
        from KratosMultiphysics.testing.utilities import GetPython3Command
        for file_path in Path("pages").rglob("*.md"):
            AddPythonSnippetOutputs(file_path)

    # sub_itr_dir = Path("pages/1_Kratos/1_For_Users")
    # print("Creating side bar for {:s}...".format(str(sub_itr_dir)))
    # menu_info = GetDirMenuInfoFromJson(sub_itr_dir)
    # if "side_bar_name" not in menu_info.keys():
    #     raise RuntimeError("No side bar name provied. Please add it to {:s}/menu_info.json.".format(str(sub_itr_dir)))
    # default_header_dict["sidebar"] = "<!>" + menu_info["side_bar_name"]
    # list_of_strings = CreateNavigatonBar(str(sub_itr_dir), 3, default_header_dict)
    # with open("_data/sidebars/{:s}.yml".format(menu_info["side_bar_name"]), "w") as file_output:
    #     file_output.write("entries:\n")
    #     file_output.writelines(list_of_strings)