from pathlib import Path
import requests
import json

__remote_tag = "remote_"

default_header_dict = {
    "title": "<PRETTY_FILE_NAME>",
    "keywords": "",
    "tags": "[<TAGLESS_NAME>]",
    "sidebar": "",
    "summary": ""
}

file_navigation_levels = [
    "unsupported",
    "folderitems",
    "subfolderitems"
]

dir_navigation_levels = [
    "root",
    "folders",
    "subfolders"
]

spacing_information = {
    "root" :"",
    "folders" :"  ",
    "subfolders" :"      ",
    "unsupported" :"",
    "folderitems" :"    ",
    "subfolderitems" :"        "
}

def GetNavigationString(navigations_list: list, navigation_level: int) -> str:
    if navigation_level >= len(navigations_list):
        return "unsupported"
    else:
        return navigations_list[navigation_level]

def GetName(file_path: Path) -> str:
    return str(file_path.relative_to(file_path.parent))

def GetTaglessName(file_path: Path) -> str:
    return GetName(file_path).replace(__remote_tag, "")

def GetPrettyName(file_path: Path) -> str:
    file_name = GetName(file_path)

    if file_name.startswith(__remote_tag):
        file_name = file_name[len(__remote_tag):]

    file_name = file_name.replace("_", " ")
    if file_name.rfind(".") != -1:
        file_name = file_name[:file_name.rfind(".")]

    if file_name.startswith(__remote_tag):
        file_name = file_name[len(__remote_tag):]

    return file_name

def ReplaceTagDataInDict(file_path: Path, input_dict: dict) -> dict:
    copied_default_entries = dict(input_dict)
    replacement_tag_dict = {
        "<PRETTY_FILE_NAME>": GetPrettyName(file_path),
        "<FILE_NAME>": str(file_path.relative_to(file_path.parent)),
        "<ABSOLUTE_FILE_NAME>": str(file_path.absolute()),
        "<TAGLESS_NAME>": GetTaglessName(file_path)
    }

    for k1 in copied_default_entries.keys():
        for k2, v2 in replacement_tag_dict.items():
            copied_default_entries[k1] = copied_default_entries[k1].replace(k2, v2)
    return copied_default_entries

def GetPageHeader(file_path: Path):
    with open(str(file_path), "r") as file_input:
        lines = file_input.readlines()

    found_starting_line = False
    found_closing_line = False
    header_lines = []
    content = []
    for i, line in enumerate(lines):
        if found_starting_line and line.startswith("---"):
            found_closing_line = True
            content = lines[i+1:]
            break
        if found_starting_line:
            header_lines.append(line)
        if line.startswith("---"):
            found_starting_line = True

    header_dict = {}
    if found_starting_line and found_closing_line:
        for header_line in header_lines:
            data = header_line[:-1].split(":")
            header_dict[data[0].strip()] = data[1].strip()
    else:
        content = lines

    return header_dict, content

def WritePageHeader(file_path: Path, default_header_dict: dict):
    header, content = GetPageHeader(file_path)
    modified_default_header_dict = ReplaceTagDataInDict(file_path, default_header_dict)

    for k, v in modified_default_header_dict.items():
        update_key = False
        tagless_v = v.replace("<!>", "")
        if k not in header.keys():
            update_key = True
            print("--- Adding missing entry \"{:s}\" with value \"{:s}\" at {:s}.".format(k, tagless_v, str(file_path)))

        if v.startswith("<!>"):
            if not update_key and tagless_v != header[k]:
                update_key = True
                print("--- Forcefully changing \"{:s}\" entry value \"{:s}\" to \"{:s}\" at {:s}.".format(k, header[k], tagless_v, str(file_path)))

        if update_key:
            header[k] = tagless_v

    with open(str(file_path), "w") as file_output:
        file_output.write("---\n")
        for k, v in header.items():
            file_output.write("{:s}: {:s}\n".format(k, v))
        file_output.write("---\n")
        file_output.writelines(content)

def GetDirMenuInfoFromJson(dir_path: Path) -> dict:
    menu_data = {}
    if (dir_path / "menu_info.json").is_file():
        with open(str(dir_path / "menu_info.json"), "r") as file_input:
            menu_data = dict(json.loads(file_input.read()))
    return menu_data

def GenerateEntryDataFromDir(dir_path: Path, navigation_level: int) -> dict:
    return {
        "title": GetPrettyName(dir_path),
        "output": "web",
        "path": dir_path,
        "type": GetNavigationString(dir_navigation_levels, navigation_level),
        "navigation_level": navigation_level
    }

def GenerateEntryDataFromFile(file_path: Path, navigation_level: int, default_header_dict: dict) -> dict:
    WritePageHeader(file_path, default_header_dict)
    header_dict, _ = GetPageHeader(file_path)

    entry_dict = {}
    if "title" not in list(header_dict.keys()):
        raise RuntimeError("title tag is not found in {:s}.".format(str(file_path)))

    entry_dict["title"] = header_dict["title"]
    entry_dict["output"] = "web"
    entry_dict["path"] = file_path
    entry_dict["type"] = GetNavigationString(file_navigation_levels, navigation_level)
    entry_dict["navigation_level"] = navigation_level
    return entry_dict

def GenerateEntryDataFromKratosExampleUrl(file_path: Path, custom_entry: dict, navigation_level: int, default_header_dict: dict) -> dict:
    raw_url = custom_entry["raw_url"]
    source_url = ""
    if "source_url" in custom_entry.keys():
        source_url = custom_entry["source_url"]
    original_folder_url = raw_url[:raw_url.rfind("/")]

    file_name = GetName(file_path)
    if not file_name.startswith(__remote_tag):
        file_path = file_path.parent / (__remote_tag + file_name)

    if not raw_url.startswith("https://raw.githubusercontent.com"):
        raise RuntimeError("Please provide the raw github raw_url. [ Provided raw_url = {:s} ].".format(raw_url))

    folder_url = raw_url[:raw_url.rfind("/")]

    print("Downloading data from: " + raw_url)
    r = requests.get(raw_url, allow_redirects=True)

    with open(str(file_path), "w") as file_output:
        if r.status_code == 200:
            data = r.text
            data = data.replace("<img src=\"", "<img src=\"{:s}/".format(folder_url))
            file_output.write(data)
        else:
            raise RuntimeError("Could not download {:s} [ Status = {:d}].".format(raw_url, r.status_code))
        if source_url != "":
            file_output.write("\n\n## Source: \n[{:s}]({:s})\n".format(source_url, source_url))
        print("--- Writing downloaded data to: " + str((file_path).absolute()))
    return GenerateEntryDataFromFile(file_path, navigation_level, default_header_dict)

def GenerateEntryDataFromExternalUrl(url_dict: dict, navigation_level: int) -> dict:
    if "title" not in url_dict.keys():
        raise RuntimeError("\"title\" key is not found in url dict {:s}.".format(str(url_dict)))
    if "url" not in url_dict.keys():
        raise RuntimeError("\"url\" key is not found in url dict {:s}.".format(str(url_dict)))
    entry_dict = {}
    entry_dict["title"] = url_dict["title"]
    entry_dict["output"] = "web"
    entry_dict["path"] = url_dict["url"]
    entry_dict["url"] = url_dict["url"]
    entry_dict["type"] = GetNavigationString(file_navigation_levels, navigation_level)
    entry_dict["navigation_level"] = navigation_level
    return entry_dict

def IsLeafEntry(entry_dict: dict) -> bool:
    return entry_dict["type"] in file_navigation_levels

def AddMissingUrlEntry(entry_info: dict, url_prefix: str):
    if "url" not in entry_info.keys():
        file_path = entry_info["path"]
        if file_path.is_file():
            if str(file_path).endswith(".md"):
                entry_info["url"] = "/" + url_prefix + str(file_path)[:-2] + "html"
            else:
                raise RuntimeError("Final entry {:s} is not a mark down file.".format(str(file_path)))

def CreateNavigationBarEntry(entry_info: dict) -> str:
    if "title" not in entry_info.keys():
        raise Exception("title is not found in entry {:s}".format(entry_info))
    # if "url" not in entry_info.keys():
    #     file_path = entry_info["path"]
    #     if file_path.is_file():
    #         if str(file_path).endswith(".md"):
    #             entry_info["url"] = "/" + str(file_path)[:-2] + "html"
    #         else:
    #             raise RuntimeError("Final entry {:s} is not a mark down file.".format(str(file_path)))
    if "url" in entry_info.keys() and entry_info["url"] == "<dummy>":
        del entry_info["url"]

    entry_string = "<TABBING>- title: {:s}\n".format(entry_info["title"])
    entry_order_list = ["product", "version", "url", "output", "folders", "subfolders", "folderitems", "subfolderitems"]
    for entry_order_item in entry_order_list:
        if entry_order_item in entry_info.keys():
            entry_string += "<TABBING>  {:s}: {:s}\n".format(entry_order_item, entry_info[entry_order_item])

    return entry_string.replace("<TABBING>", spacing_information[entry_info["type"]])

if __name__ == "__main__":
    print(GetPrettyName(Path("/test1/test2/hello_test.md")))
    print(GetPrettyName(Path("/test1/test2/!remote_hello_test.md")))
    WritePageHeader(
        Path("pages/2_Applications/Shape_Optimization_Application/99_Examples/02_Strain_Energy_Minimization_3D_Shell.md"),
        {
            "sidebar": "<!>Hellothere"
        })

    data = GenerateEntryDataFromKratosExampleUrl(
        Path("online_test.md"),
        "https://github.com/KratosMultiphysics/Examples/tree/master/shape_optimization/use_cases/11_Shape_Update_Optimization_Stanford_Bunny/README.md",
        "none",
        default_header_dict)
    print(data)