'''HDF5 file IO.

license: HDF5Application/license.txt
'''


import os
import typing
from pathlib import Path


import KratosMultiphysics
import KratosMultiphysics.HDF5Application as KratosHDF5
import KratosMultiphysics.kratos_utilities as kratos_utils
from .pattern import EvaluatePattern
from .pattern import GetMatchingEntities
from .pattern import PathPatternEntity

def GetFileName(model_part: KratosMultiphysics.ModelPart, file_name_pattern: str, time_format: str) -> str:
    if not file_name_pattern.endswith(".h5"):
        file_name_pattern += ".h5"
    return EvaluatePattern(file_name_pattern, model_part, time_format)

def KeepMostRecentFiles(file_name_pattern: str, max_files_to_keep: int, number_of_oldest_files_to_keep: int) -> None:
    if number_of_oldest_files_to_keep > max_files_to_keep:
        raise RuntimeError(f"Max files to keep should be higher than or equal to number_of_oldest_files_to_keep.")


    list_of_sorted_file_names = sorted(list(GetMatchingEntities(PathPatternEntity(Path(".")),
                                                   file_name_pattern,
                                                   {"<step>": int, "<time>": float})),
                                                   key=lambda args: tuple(args[1:]))

    if (len(list_of_sorted_file_names) >= max_files_to_keep):
        for file_name in list_of_sorted_file_names[number_of_oldest_files_to_keep:len(list_of_sorted_file_names) - max_files_to_keep + number_of_oldest_files_to_keep]:
            kratos_utils.DeleteFileIfExisting(str(file_name))

def CreateHDF5File(model_part: KratosMultiphysics.ModelPart, settings: KratosMultiphysics.Parameters) -> KratosHDF5.HDF5File:
    defaults = KratosMultiphysics.Parameters("""{
        "file_name"                     : "",
        "time_format"                   : "0.4f",
        "file_driver"                   : "sec2",
        "file_access_mode"              : "exclusive",
        "echo_level"                    : 0,
        "max_files_to_keep"             : "unlimited",
        "number_of_oldest_files_to_keep": 1
    }""")
    if model_part.IsDistributed():
        defaults["file_driver"].SetString("mpio")
    elif os.name == "nt":
        defaults["file_driver"].SetString("windows")
    settings.AddMissingParameters(defaults)

    current_file_settings = settings.Clone()

    # generate the new name from pattern
    file_name_pattern = current_file_settings["file_name"].GetString()
    file_name = GetFileName(model_part, file_name_pattern, current_file_settings["time_format"].GetString())
    current_file_settings["file_name"].SetString(file_name)

    # create directories
    dirname = Path(file_name).absolute().parent
    KratosMultiphysics.FilesystemExtensions.MPISafeCreateDirectories(str(dirname))

    if current_file_settings["file_access_mode"].GetString() in ["truncate", "exclusive", "read_write"]:
        # this is a writing file, hence check whether max files to keep is used.
        if current_file_settings.Has("max_files_to_keep") and current_file_settings["max_files_to_keep"].IsInt():
            KeepMostRecentFiles(file_name_pattern, current_file_settings["max_files_to_keep"].GetInt(), current_file_settings["number_of_oldest_files_to_keep"].GetInt())

    return KratosHDF5.HDF5File(model_part.GetCommunicator().GetDataCommunicator(), current_file_settings)

class OpenHDF5File:
    def __init__(self, parameters: KratosMultiphysics.Parameters, model_part: KratosMultiphysics.ModelPart) -> None:
        self.__model_part = model_part
        self.__parameters = parameters
        self.__file: 'typing.Optional[KratosHDF5.HDF5File]' = None

    def __enter__(self) -> KratosHDF5.HDF5File:
        self.__file = CreateHDF5File(self.__model_part, self.__parameters)
        return self.__file

    def __exit__(self, exit_type, exit_value, exit_traceback) -> None:
        self.__file.Close()
        self.__file = None

    @property
    def file(self) -> KratosHDF5.HDF5File:
        return self.__file

