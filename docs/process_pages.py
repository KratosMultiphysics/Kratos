import json
from pathlib import Path
from posixpath import split
import requests

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

def GetHeaderAndContent(current_file_path: Path) -> list:
    with open(str(current_file_path), "r") as file_input:
        lines = file_input.readlines()

    found_starting_line = False
    found_closing_line = False
    header = []
    content = []
    for i, line in enumerate(lines):
        if found_starting_line and line.startswith("---"):
            found_closing_line = True
            content = lines[i+1:]
            break
        if found_starting_line:
            header.append(line)
        if line.startswith("---"):
            found_starting_line = True

    if found_starting_line and found_closing_line:
        return header, content
    else:
        return [], lines

def GetPageHeader(header_lines: list) -> dict:
    header_dict = {}
    for header_line in header_lines:
        data = header_line[:-1].split(":")
        header_dict[data[0].strip()] = data[1].strip()

    return header_dict

# def WritePageHeader(header_lines: list, content_lines:str, header_dict: dict):
#     pass

def GetName(file_path: Path) -> str:
    return str(file_path.relative_to(file_path.parent))

def GetPrettyName(file_path: Path) -> str:
    file_name = GetName(file_path)
    file_name = file_name.replace("_", " ")
    return file_name

def GenerateEntryDataFromFile(current_file_path: Path, entry_type: str) -> dict:
    current_dict = {}
    with open(str(current_file_path), "r") as file_input:
        lines = file_input.readlines()
    for line in lines:
        if line.startswith("title:"):
            current_dict["title"] = line[6:-1].strip()
    # check for the header, if not present add the header

    if "title" not in list(current_dict.keys()):
        current_dict["title"] = GetPrettyName(current_file_path)

    current_dict["output"] = "web"
    current_dict["path"] = current_file_path
    current_dict["type"] = entry_type
    return current_dict

def GenerateEntryDataFromDir(current_dir_path: Path, entry_type: str) -> dict:
    return {
        "title": GetPrettyName(current_dir_path),
        "output": "web",
        "path": current_dir_path,
        "type": entry_type
    }

def GenerateEntryDataFromUrl(current_file_path: Path, url: str, entry_type: str) -> dict:
    raw_url = url
    original_folder_url = url[:url.rfind("/")]

    # get the raw url in case of this is from github
    if url.startswith("https://github.com"):
        tree_index = url.find("/tree/")
        raw_url = "https://raw.githubusercontent.com" + url[len("https://github.com"):tree_index] + url[tree_index+5:]

    folder_url = raw_url[:raw_url.rfind("/")]

    print("Downloading data from: " + url)
    r = requests.get(raw_url, allow_redirects=True)

    with open(str(current_file_path), "w") as file_output:
        if r.status_code == 200:
            data = r.text
            data = data.replace("<img src=\"", "<img src=\"{:s}/".format(folder_url))
            file_output.write(data)
        file_output.write("\n\n## Source: \n[{:s}]({:s})\n".format(original_folder_url, original_folder_url))
        print("Writing downloaded data to: " + str((current_file_path).absolute()))
    return GenerateEntryDataFromFile(current_file_path, entry_type)

def GetEntryPathFromString(input_str: str, current_path: Path) -> Path:
    if input_str.startswith("/"):
        return docs_absolute_path / input_str[1:]
    elif input_str.startswith("http"):
        return Path()
    else:
        return current_path / input_str

def CreateNavigationBarEntry(entry_info: dict) -> str:
    if "title" not in entry_info.keys():
        raise Exception("title is not found in entry {:s}".format(entry_info))
    entry_string = "<TABBING>- title: {:s}\n".format(entry_info["title"])
    entry_order_list = ["product", "version", "url", "output", "folders", "subfolders", "folderitems", "subfolderitems"]
    for entry_order_item in entry_order_list:
        if entry_order_item in entry_info.keys():
            entry_string += "<TABBING>  {:s}: {:s}\n".format(entry_order_item, entry_info[entry_order_item])

    return entry_string

def CreateEntriesDicts(current_path: Path, navigation_level, max_navigation_level, side_bar: str) -> list:
    if navigation_level == max_navigation_level:
        return []

    if navigation_level >= len(spacing_information["files"]):
        file_entry_type = "unsupported"
    else:
        file_entry_type = spacing_information["files"][navigation_level][0]

    if navigation_level >= len(spacing_information["dirs"]):
        dir_entry_type = "unsupported"
    else:
        dir_entry_type = spacing_information["dirs"][navigation_level][0]

    menu_data = {}
    if (current_path / "menu_info.json").is_file():
        with open(str(current_path / "menu_info.json"), "r") as file_input:
            menu_data = dict(json.loads(file_input.read()))

    custom_entries_order = []
    if "custom_entries" in menu_data.keys():
        custom_entries_order = menu_data["custom_entries"]

    ignore_entries = []
    if "ignore_entries" in menu_data.keys():
        ignore_entries = menu_data["ignore_entries"]

    list_of_entries = [GenerateEntryDataFromDir(current_path, dir_entry_type)]

    for i, custom_entry in enumerate(custom_entries_order):
        if custom_entry in ignore_entries:
            raise RuntimeError("Custom entry {:s} found in ignore entries list. These \"custom_entries\" and \"ignore_entries\" should be mutually exclusive.".format(custom_entry))

        if isinstance(custom_entry, dict):
            current_dict = GenerateEntryDataFromUrl(current_path / custom_entry["file_name"], custom_entry["url"], side_bar, file_entry_type)
            custom_entries_order[i] = custom_entry["file_name"]
            list_of_entries.append(current_dict)
        else:
            current_entry_path = GetEntryPathFromString(custom_entry, current_path)
            if current_entry_path.is_file():
                list_of_entries.append(GenerateEntryDataFromFile(current_entry_path, side_bar, file_entry_type))
            elif current_entry_path.is_dir():
                list_of_entries.extend(CreateEntriesDicts(current_entry_path, navigation_level + 1, max_navigation_level))
            else:
                raise RuntimeError("Custom entry {:s} is not a file or a dir.".format(str(current_entry_path)))

    temp_entries_list = []
    for iter_path in current_path.iterdir():
        if iter_path.is_file() and str(iter_path).endswith(".md"):
            file_name = GetName(iter_path)
            if file_name not in custom_entries_order and file_name not in ignore_entries:
                temp_entries_list.append(GenerateEntryDataFromFile(iter_path, side_bar, file_entry_type))

    list_of_entries.extend(sorted(temp_entries_list, key=lambda x: x["title"]))

    temp_entries_list = []
    for iter_path in current_path.iterdir():
        if iter_path.is_dir():
            if not GetName(iter_path) in custom_entries_order:
                temp_entries_list.append(iter_path)

    temp_entries_list = sorted(temp_entries_list, key= lambda x: str(x))
    for iter_path in temp_entries_list:
        list_of_entries.extend(CreateEntriesDicts(iter_path, navigation_level + 1, max_navigation_level, side_bar))

    return list_of_entries

if __name__ == "__main__":
    print(GetPageHeader(Path("pages/2_Applications/Shape_Optimization_Application/99_Examples/01_Strain_Energy_Minimization_3D_Hook.md")))
    # list_of_entries = CreateEntriesDicts(Path("pages/2_Applications/Shape_Optimization_Application"), 0, 10, "shape_optimization_application")
    # for entry in list_of_entries:
    #     print(str(entry["path"]))