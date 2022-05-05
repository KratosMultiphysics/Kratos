from pathlib import Path
from collections import OrderedDict
import json

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

def GetAbsolutePath(input: str):
    if input[0] == "/":
        return input
    else:
        return "/" + input

def GetFileIndex(file_name: Path) -> int:
    relative_file_name = str(file_name.relative_to(file_name.parent))
    pos = relative_file_name.find("_")
    index = -1
    if pos != -1 and relative_file_name[:pos].isnumeric():
        index = int(relative_file_name[:pos])
    return index

def GetPrettyName(file_path: Path) -> str:
    file_name = str(file_path.relative_to(file_path.parent))
    if file_path.is_file():
        char_index = file_name.rfind(".")
        if char_index != -1:
            file_name = file_name[:char_index]

    char_index = file_name.find("_")
    if char_index != -1 and file_name[:char_index].isnumeric():
        file_name = file_name[char_index+1:]
    file_name = file_name.replace("_", " ")
    return file_name

def CreateNavigationBarEntry(entry_info: dict) -> str:
    if "title" not in entry_info.keys():
        raise Exception("title is not found in entry {:s}".format(entry_info))
    entry_string = "<TABBING>- title: {:s}\n".format(entry_info["title"])
    entry_order_list = ["product", "version", "url", "output", "folders", "subfolders", "folderitems", "subfolderitems"]
    for entry_order_item in entry_order_list:
        if entry_order_item in entry_info.keys():
            entry_string += "<TABBING>  {:s}: {:s}\n".format(entry_order_item, entry_info[entry_order_item])

    return entry_string

def AddTabbingToEntry(entry_string: str, spaces: str) -> str:
    return entry_string.replace("<TABBING>", spaces)

def GenerateListOfSortedRawEntries(current_dir: Path, is_valid_path) -> list[Path]:
    list_of_raw_entries_dict = {}
    max_index = -1
    for iter_dir in current_dir.iterdir():
        if is_valid_path(iter_dir):
            index = GetFileIndex(iter_dir)
            max_index = max(max_index, index)
            if index not in list_of_raw_entries_dict.keys():
                list_of_raw_entries_dict[index] = []
            list_of_raw_entries_dict[index].append(iter_dir)

    if -1 in list_of_raw_entries_dict.keys():
        list_of_raw_entries_dict[max_index + 1] = list(list_of_raw_entries_dict[-1])
        del list_of_raw_entries_dict[-1]

    list_of_raw_entries_list = []
    sorted_keys_list = sorted(list(list_of_raw_entries_dict.keys()))
    for sorted_key in sorted_keys_list:
        sorted_sub_list = sorted(list_of_raw_entries_dict[sorted_key], key=lambda x: str(x))
        list_of_raw_entries_list.extend(sorted_sub_list)

    return list_of_raw_entries_list

def GetEntryDict(current_path: Path) -> dict:
    entry_dict = {}
    if current_path.is_file():
        # get file information
        with open(str(current_path), "r") as file_input:
            lines = file_input.readlines()
        for line in lines:
            if line.startswith("title:"):
                entry_dict["title"] = line[6:-1].strip()
        entry_dict["url"] = GetAbsolutePath(str(current_path))[:-2] + "html"
    else:
        entry_dict["title"] = GetPrettyName(current_path)

    entry_dict["output"] = "web"
    entry_dict["path"] = current_path
    return entry_dict

def GetDirEntryDictFromJson(current_path: Path) -> dict:
    json_dict = {}
    if current_path.is_dir():
        if (current_path / "info.json").is_file():
            with open(str(current_path / "info.json"), "r") as file_input:
                json_dict = json.loads(file_input.read())
        return json_dict
    else:
        raise Exception("{:s} is not a valid directory.".format(str(current_path)))

def CreateNavigationBarStructure(
    current_dict: dict,
    navigation_level: int,
    max_navigation_level: int) -> list[dict]:

    list_of_dicts = []
    current_path = current_dict["path"]

    if navigation_level == max_navigation_level:
        file_tabbing_info = spacing_information["files"][navigation_level - 1]
        dir_json_info = GetDirEntryDictFromJson(current_path)
        if "url" in dir_json_info.keys():
            current_dict["url"] = GetAbsolutePath(str(current_path / dir_json_info["url"]))[:-2] + "html"
        current_dict["type"] = file_tabbing_info[0]
        current_dict["str"] = AddTabbingToEntry(CreateNavigationBarEntry(current_dict), file_tabbing_info[1])
        list_of_dicts.append(current_dict)
    else:
        file_tabbing_info = spacing_information["files"][navigation_level]
        dir_tabbing_info = spacing_information["dirs"][navigation_level]
        current_dict["type"] = dir_tabbing_info[0]
        current_dict["str"] = AddTabbingToEntry(CreateNavigationBarEntry(current_dict), dir_tabbing_info[1])
        list_of_dicts.append(current_dict)
        list_of_dir_entries = GenerateListOfSortedRawEntries(current_path, lambda x: x.is_dir())

        list_of_file_entries = GenerateListOfSortedRawEntries(current_path, lambda x: x.is_file() and str(x).endswith(".md"))
        for file_entry in list_of_file_entries:
            file_dict = GetEntryDict(file_entry)
            file_dict["type"] = file_tabbing_info[0]
            file_dict["str"] = AddTabbingToEntry(CreateNavigationBarEntry(file_dict), file_tabbing_info[1])
            list_of_dicts.append(file_dict)

        if len(list_of_file_entries) == 0 and navigation_level > 0:
            # add a dummy file entry
            dummy_dict = {
                "title": "",
                "path": Path("."),
                "output": "web",
                "type": file_tabbing_info[0]
            }
            dummy_dict["str"] = AddTabbingToEntry(CreateNavigationBarEntry(dummy_dict), file_tabbing_info[1])
            list_of_dicts.append(dummy_dict)

        for dir_path in list_of_dir_entries:
            dir_dict = GetEntryDict(dir_path)
            list_of_dicts.extend(CreateNavigationBarStructure(dir_dict, navigation_level + 1, max_navigation_level))

    # clean empty dummy dicts which are not required
    current_index = 0
    while (current_index < len(list_of_dicts) - 1):
        if list_of_dicts[current_index]["title"] == "" and list_of_dicts[current_index]["type"] == list_of_dicts[current_index+1]["type"]:
            del list_of_dicts[current_index]
        else:
            current_index += 1

    return list_of_dicts

def GetTypeInfo(item_type: str) -> str:
    for spacing_info in spacing_information["dirs"]:
        if item_type == spacing_info[0]:
            return spacing_info[1]
    for spacing_info in spacing_information["files"]:
        if item_type == spacing_info[0]:
            return spacing_info[1]

    raise Exception("Entry type {:s} is not defined.".format(item_type))

def GenerateStrings(list_of_dicts: list[dict]) -> list[str]:
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

def GetDataFromFileEntry(file_path: Path):
    if file_path.is_file() and str(file_path).endswith(".md"):
        with open(str(file_path), "r") as file_input:
            lines = file_input.readlines()
        meta_data_dict = {}
        if lines[0] == "---\n":
            closing_index = lines[1:].index("---\n") + 1
            for line in lines[1:closing_index]:
                line_data = line[:-1].strip().split(":")
                meta_data_dict[line_data[0].strip()] = line_data[1].strip()
            return meta_data_dict, lines[closing_index+1:]
        else:
            return meta_data_dict, lines
    else:
        return None, None

def FillDataFileEntries(file_path: Path, default_entries: OrderedDict):
    meta_data_dict, content = GetDataFromFileEntry(file_path)
    copied_default_entries = dict(default_entries)

    replacement_tag_dict = {
        "<PRETTY_FILE_NAME>": GetPrettyName(file_path),
        "<FILE_NAME>": str(file_path.relative_to(file_path.parent)),
        "<ABSOLUTE_FILE_NAME>": str(file_path.absolute())
    }

    for k1 in copied_default_entries.keys():
        for k2, v2 in replacement_tag_dict.items():
            copied_default_entries[k1] = copied_default_entries[k1].replace(k2, v2)

    # check if it is a file with content
    if content is not None:
        if meta_data_dict is not None:
            # now add missing entries from the defaults
            for k, v in copied_default_entries.items():
                if v.startswith("<!>"):
                    v = v[3:]
                    if v != meta_data_dict[k]:
                        print("Warning: {:s} of file at {:s} is \"{:s}\" which is different from forced defaults of \"{:s}\". Overwriting with forced default value.".format(k, str(file_path), meta_data_dict[k], v))
                    meta_data_dict[k] = v
                else:
                    if not k in meta_data_dict.keys():
                        meta_data_dict[k] = v
        else:
            meta_data_dict = copied_default_entries

        # now create ordered dict from the following order
        ordered_meta_data_dict = OrderedDict()
        for k in default_entries.keys():
            ordered_meta_data_dict[k] = meta_data_dict[k]

        with open(str(file_path), "w") as file_output:
            file_output.write("---\n")

            for k, v in ordered_meta_data_dict.items():
                file_output.write("{:s}: {:s}\n".format(k, v))

            file_output.write("---\n")
            file_output.writelines(content)

def GetFilesList(current_path: Path) -> list[Path]:
    list_of_files = []
    for iter_dir in current_path.iterdir():
        if iter_dir.is_file():
            list_of_files.append(iter_dir)
        else:
            list_of_files.extend(GetFilesList(iter_dir))

    return list_of_files

if __name__ == "__main__":
    list_of_files = GetFilesList(Path("pages"))
    defaults_dict = OrderedDict()
    defaults_dict["title"] = "<PRETTY_FILE_NAME>"
    defaults_dict["keywords"] = ""
    defaults_dict["tags"] = "[<FILE_NAME>]"
    defaults_dict["sidebar"] = "<DEFAULT_SIDE_BAR>"
    defaults_dict["summary"] = ""

    print("Processing pages for meta data completeness...")
    for file in list_of_files:
        FillDataFileEntries(file, defaults_dict)

    print("Creating top navigation bar...")
    # generate top navigation bar
    with open("_data/topnav.yml.orig", "r") as file_input:
        lines = file_input.readlines()
    list_of_entries = CreateNavigationBarStructure(
        {
            "title": "Topnav dropdowns",
            "path" : Path("pages")
        }, 0, 2)
    list_of_entries = GenerateStrings(list_of_entries)
    lines.extend(list_of_entries)

    with open("_data/topnav.yml", "w") as file_output:
        file_output.writelines(lines)

    for iter_dir in Path("pages").iterdir():
        if iter_dir.is_dir():
            for sub_itr_dir in iter_dir.iterdir():
                if sub_itr_dir.is_dir():
                    print("Creating side bar for {:s}...".format(str(sub_itr_dir)))
                    root_dict = GetEntryDict(sub_itr_dir)
                    json_settings = GetDirEntryDictFromJson(sub_itr_dir)
                    if "product" in json_settings:
                        root_dict["product"] = json_settings["product"]
                    else:
                        root_dict["product"] = root_dict["title"]
                    root_dict["title"] = "sidebar"

                    if "file_name" in json_settings:
                        file_name = json_settings["file_name"]
                    else:
                        file_name = str(root_dict["path"]).replace("/", "_")

                    list_of_files = GetFilesList(sub_itr_dir)
                    # forcefully update side bar with <!> prefix flag
                    defaults_dict["sidebar"] = "<!>{:s}".format(file_name)
                    for file in list_of_files:
                        FillDataFileEntries(file, defaults_dict)

                    list_of_entries = CreateNavigationBarStructure(
                            root_dict, 0, 3)

                    list_of_entries = GenerateStrings(list_of_entries)

                    lines = ["entries:\n"]
                    lines.extend(list_of_entries)

                    with open("_data/sidebars/{:s}.yml".format(file_name), "w") as file_output:
                        file_output.writelines(lines)
