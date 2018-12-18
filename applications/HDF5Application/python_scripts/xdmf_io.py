from hdf5_io import IOObject
from create_xdmf_file import WriteXdmfFile
import os

class XdmfOutput(IOObject):
    """Output that creates the xdmf-file for the given h5-files"""

    def Execute(self, model_part, hdf5_file): pass

    def ExecutePostOutput(self, model_part, hdf5_file_name):
        model_part.GetCommunicator().Barrier()
        if model_part.GetCommunicator().MyPID() == 0:
            # write xdmf only on one rank!
            file_path, file_name = os.path.split(hdf5_file_name)
            base_file_name = "-".join(file_name.split("-")[:-1]) + ".h5"
            WriteXdmfFile(base_file_name, file_path)

