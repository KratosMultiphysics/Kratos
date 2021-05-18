'''HDF5 operating system operations.

license: HDF5Application/license.txt

Main authors:
    Philipp Bucher
    Michael Andre
'''


import os


import KratosMultiphysics
from KratosMultiphysics.kratos_utilities import DeleteFileIfExisting


class DeleteOldH5Files:
    '''Delete h5-files from previous simulations.'''

    def __call__(self, model_part, hdf5_file):
        file_path, file_name = os.path.split(hdf5_file.GetFileName())
        time_prefix = file_name.replace(".h5", "") + "-"
        current_time = model_part.ProcessInfo[KratosMultiphysics.TIME]
        if file_path == "":
            file_path = "."  # os.listdir fails with empty path
        for name in os.listdir(file_path):
            if name.startswith(time_prefix):
                file_time = float(name.replace(".h5", "")[len(time_prefix):])
                if file_time > current_time:
                    DeleteFileIfExisting(
                        os.path.join(file_path, name))


def Create(settings):
    '''Return an operation specified by the setting's 'operation_type'.

    This method is normally not used directly, but rather it is imported
    in core.operations.model_part.Create using the 'module_name' setting.
    '''
    operation_type = settings['operation_type']
    if operation_type == 'delete_old_h5_files':
        return DeleteOldH5Files()
    else:
        raise ValueError(
            '"operation_type" has invalid value "' + operation_type + '"')
