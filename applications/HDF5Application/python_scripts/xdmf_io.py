from hdf5_io import IOObject
from create_xdmf_file import WriteXdmfFile
import os

class XdmfOutput(IOObject):
    """Output that creates the xdmf-file for the given h5-files"""

    def Execute(self, model_part, hdf5_file): pass

    def ExecuteAfterClose(self, model_part, hdf5_file_name):
        model_part.GetCommunicator().Barrier()
        # write xdmf only on one rank!
        if model_part.GetCommunicator().MyPID() == 0:
            # removing the time-label in the file-name
            base_hdf5_file_name = "-".join(hdf5_file_name.split("-")[:-1]) + ".h5"
            WriteXdmfFile(base_hdf5_file_name)

