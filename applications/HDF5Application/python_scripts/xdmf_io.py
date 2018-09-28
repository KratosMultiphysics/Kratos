import KratosMultiphysics
from hdf5_io import IOObject
from create_xdmf_file import WriteXdmfFile
import os



class XdmfOutput(IOObject):
    """Provides the interface for writing a model part to a file."""

    def __init__(self, settings):
        default_settings = KratosMultiphysics.Parameters("""
        {
            "file_name" : ""
        }
        """)

        self.settings = settings.Clone()
        self.settings.ValidateAndAssignDefaults(default_settings)

    def Execute(self, model_part, hdf5_file):
        path, file_name = os.path.split(hdf5_file.GetFileName())
        print(os.path.split(hdf5_file.GetFileName()))
        print(file_name[:7])
        print(os.getcwd())
        WriteXdmfFile(file_name[:7]+".h5", os.path.join(os.getcwd(),path))

