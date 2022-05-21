from pathlib import Path
import requests

from docs.old_process_pages import GetFileIndex

__remote_tag = "!remote_"

def GetName(file_path: Path) -> str:
    return str(file_path.relative_to(file_path.parent))

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
        "<ABSOLUTE_FILE_NAME>": str(file_path.absolute())
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
            print("Adding missing entry \"{:s}\" with value \"{:s}\" at {:s}.".format(k, tagless_v, str(file_path)))

        if v.startswith("<!>"):
            if not update_key and tagless_v != header[k]:
                update_key = True
                print("Forcefully changing \"{:s}\" entry value \"{:s}\" to \"{:s}\" at {:s}.".format(k, header[k], tagless_v, str(file_path)))

        if update_key:
            header[k] = tagless_v

    with open(str(file_path), "w") as file_output:
        file_output.write("---\n")
        for k, v in header.items():
            file_output.write("{:s}: {:s}\n".format(k, v))
        file_output.write("---\n")
        file_output.writelines(content)

def GenerateEntryDataFromDir(current_dir_path: Path, entry_type: str) -> dict:
    return {
        "title": GetPrettyName(current_dir_path),
        "output": "web",
        "path": current_dir_path,
        "type": entry_type
    }

def GenerateEntryDataFromFile(current_file_path: Path, entry_type: str) -> dict:
    header_dict, _ = GetPageHeader(current_file_path)

    current_dict = {}
    if "title" not in list(header_dict.keys()):
        current_dict["title"] = GetPrettyName(current_file_path)

    current_dict["output"] = "web"
    current_dict["path"] = current_file_path
    current_dict["type"] = entry_type
    return current_dict

def GenerateEntryDataFromUrl(current_file_path: Path, url: str, entry_type: str) -> dict:
    raw_url = url
    original_folder_url = url[:url.rfind("/")]

    file_name = GetName(current_file_path)
    if not file_name.startswith(__remote_tag):
        current_file_path = current_file_path.parent / (__remote_tag + file_name)

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

if __name__ == "__main__":
    print(GetPrettyName(Path("/test1/test2/hello_test.md")))
    print(GetPrettyName(Path("/test1/test2/!remote_hello_test.md")))
    WritePageHeader(
        Path("pages/2_Applications/Shape_Optimization_Application/99_Examples/02_Strain_Energy_Minimization_3D_Shell.md"),
        {
            "sidebar": "<!>Hellothere"
        })

    GenerateEntryDataFromUrl()