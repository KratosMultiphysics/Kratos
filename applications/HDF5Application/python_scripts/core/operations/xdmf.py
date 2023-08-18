'''HDF5 Xdmf operations.

license: HDF5Application/license.txt

Main authors:
    Philipp Bucher
    Michael Andre
'''

# --- Kratos Imports ---
import KratosMultiphysics
import KratosMultiphysics.HDF5Application as KratosHDF5


# in case the h5py-module is not installed (e.g. on clusters) we don't want it to crash the simulation!
# => in such a case the xdmf can be created manually afterwards locally
try:
    import KratosMultiphysics.HDF5Application.xdmf_utils
    from KratosMultiphysics.HDF5Application.xdmf_utils import WriteMultifileTemporalAnalysisToXdmf
except ImportError:
    # if we failed to import, then assign a dummy function that does nothing.
    WriteMultifileTemporalAnalysisToXdmf = lambda *args: None
    warn_msg = "XDMF-Writing is not available,\nOnly HDF5 files are written"
    KratosMultiphysics.Logger.PrintWarning(__name__, warn_msg)


class XdmfOutput(KratosMultiphysics.Operation):
    '''Output that creates the xdmf-file for the given h5-files.'''

    def __init__(self,
                 model_part: KratosMultiphysics.ModelPart,
                 _: KratosMultiphysics.Parameters,
                 file: KratosHDF5.HDF5File):
        super().__init__()
        self.__model_part = model_part
        self.__file = file

    def Execute(self) -> None:
        self.__model_part.GetCommunicator().GetDataCommunicator().Barrier()
        # write xdmf only on one rank!
        if self.__model_part.GetCommunicator().MyPID() == 0:
            WriteMultifileTemporalAnalysisToXdmf(self.__file.GetFileName(),
                                                 "/ModelData",
                                                 "/ResultsData")


def Create(settings):
    '''Return an operation specified by the setting's 'operation_type'.

    This method is normally not used directly, but rather it is imported
    in core.operations.model_part.Create using the 'module_name' setting.
    '''
    operation_type = settings['operation_type'].GetString()
    if operation_type == 'xdmf_output':
        return XdmfOutput
    else:
        raise ValueError(f'"operation_type" has invalid value "{operation_type}"')
