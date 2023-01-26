'''HDF5 operating system operations.

license: HDF5Application/license.txt

Main authors:
    Philipp Bucher
    Michael Andre
'''

# --- Kratos Imports ---
import KratosMultiphysics
import KratosMultiphysics.HDF5Application as KratosHDF5
from KratosMultiphysics.kratos_utilities import DeleteFileIfExisting

# --- STD Imports ---
import pathlib


class DeleteOldH5Files(KratosMultiphysics.Operation):
    '''Delete h5-files from previous simulations.'''

    def __init__(self,
                 model_part: KratosMultiphysics.ModelPart,
                 _: KratosMultiphysics.Parameters,
                 file: KratosHDF5.HDF5File):
        super().__init__()
        self.__model_part = model_part
        self.__file = file

    def Execute(self) -> None:
        ## @todo this is a terrible idea (@matekelemen)
        file_path = pathlib.Path(self.__file.GetFileName())
        directory_path = file_path.absolute().parent
        file_name = file_path.name
        time_prefix = file_name.replace(".h5", "") + "-"
        current_time = self.__model_part.ProcessInfo[KratosMultiphysics.TIME]
        for path in directory_path.glob("*.h5"):
            name = path.name
            if name.startswith(time_prefix):
                file_time = float(name.replace(".h5", "")[len(time_prefix):])
                if file_time > current_time:
                    DeleteFileIfExisting(str(directory_path / name))


def Create(settings):
    '''Return an operation specified by the setting's 'operation_type'.

    This method is normally not used directly, but rather it is imported
    in core.operations.model_part.Create using the 'module_name' setting.
    '''
    operation_type = settings['operation_type'].GetString()
    if operation_type == 'delete_old_h5_files':
        return DeleteOldH5Files
    else:
        raise ValueError(f'"operation_type" has invalid value "{operation_type}"')
