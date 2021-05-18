'''HDF5 Xdmf operations.

license: HDF5Application/license.txt

Main authors:
    Philipp Bucher
    Michael Andre
'''
import KratosMultiphysics

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


class XdmfOutput(object):
    '''Output that creates the xdmf-file for the given h5-files.'''

    def __call__(self, model_part, hdf5_file):
        model_part.GetCommunicator().GetDataCommunicator().Barrier()
        # write xdmf only on one rank!
        if model_part.GetCommunicator().MyPID() == 0:
            WriteMultifileTemporalAnalysisToXdmf(
                hdf5_file.GetFileName(), "/ModelData", "/ResultsData")


def Create(settings):
    '''Return an operation specified by the setting's 'operation_type'.

    This method is normally not used directly, but rather it is imported
    in core.operations.model_part.Create using the 'module_name' setting.
    '''
    operation_type = settings['operation_type']
    if operation_type == 'xdmf_output':
        return XdmfOutput()
    else:
        raise ValueError(
            '"operation_type" has invalid value "' + operation_type + '"')
