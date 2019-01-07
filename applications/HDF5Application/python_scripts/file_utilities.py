import KratosMultiphysics as KM
import KratosMultiphysics.kratos_utilities as kratos_utils
from hdf5_io import IOObject
import os

class DeleteOldH5Files(IOObject):
    """Deleting h5-files from previous simulations"""
    def Execute(self, model_part, hdf5_file):
        file_path, file_name = os.path.split(hdf5_file.GetFileName())
        time_prefix = file_name[:len(file_name)-3]+"-"
        current_time = model_part.ProcessInfo[KM.TIME]
        if file_path == "": file_path = "." # os.listdir fails with empty path
        for name in os.listdir(file_path):
            if name.startswith(time_prefix):
                file_time = float(name.replace(".h5", "")[len(time_prefix):])
                if file_time > current_time:
                    kratos_utils.DeleteFileIfExisting(os.path.join(file_path, name))


def ExtendPathWithSaveFolder(file_path):
    """extending a file-path with an extra folder for the h5-files
    e.g. "foo/bar/my_file.h5"
      to "foo/bar/my_file__h5_files/my_file.h5"

    """
    raw_path, file_name = os.path.split(file_path)

    folder_name = file_name + "__h5_files"
    folder_path = os.path.join(raw_path, folder_name)

    return os.path.join(folder_path, file_name)

def CreateFolderIfNotExisting(folder_path, communicator):
    """Create a folder if it is not existing, also works in MPI"""
    if not os.path.isdir(folder_path) and communicator.MyPID() == 0:
        os.makedirs(folder_path)
    communicator.Barrier()

def CheckIfSaveH5FilesInFolder(settings, model_part):
    """If the h5 files are to be saved in a folder then
    do the necessary manipulations of the settings
    """
    file_settings = settings["file_settings"]
    output_time_settings = settings["output_time_settings"]
    if file_settings.Has("save_h5_files_in_folder"): #todo(msandre): discuss if this should be default
        if file_settings["save_h5_files_in_folder"].GetBool():
            # modify the file-path
            if output_time_settings.Has("file_name"):
                file_name = output_time_settings["file_name"].GetString()
            else:
                # use the name of the modelpart if nothing is specified
                file_name = model_part.Name
                output_time_settings.AddEmptyValue("file_name")

            full_path = ExtendPathWithSaveFolder(file_name)

            output_time_settings["file_name"].SetString(full_path)

            folder_path = os.path.split(full_path)[0]

            CreateFolderIfNotExisting(folder_path, model_part.GetCommunicator())

        # removing since it is not in the defaults
        file_settings.RemoveValue("save_h5_files_in_folder")
