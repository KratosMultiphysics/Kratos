import KratosMultiphysics as KM
import KratosMultiphysics.HDF5Application.temporal_output_process_factory as output_factory
import KratosMultiphysics.HDF5Application.file_utilities as file_utils

def Factory(settings, Model):
    """Return a process for writing simulation results for a single mesh to HDF5.
    It also creates the xdmf-file on the fly
    """
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")

    params = settings["Parameters"]

    # setting default "file_settings"
    if not params.Has("file_settings"):
        file_params = KM.Parameters(r'''{
            "file_access_mode"      : "truncate",
            "save_h5_files_in_folder" : true
        }''')
        params.AddValue("file_settings", file_params)
    else:
        if not params["file_settings"].Has("file_access_mode"):
            params["file_settings"].AddEmptyValue("file_access_mode").SetString("truncate")
        if not params["file_settings"].Has("save_h5_files_in_folder"):
            params["file_settings"].AddEmptyValue("save_h5_files_in_folder").SetBool(True)

    model_part_name = params["model_part_name"].GetString() # name of modelpart must be specified!
    #todo(msandre): collapse older partitioned scripts to their serial counterparts like this
    if Model[model_part_name].GetCommunicator().TotalProcesses() > 1:
        factory_helper = output_factory.PartitionedTemporalOutputFactoryHelper()
    else:
        factory_helper = output_factory.TemporalOutputFactoryHelper()

    (temporal_output_process, _, list_of_results_output) = factory_helper.Execute(params, Model)
    for results_output in list_of_results_output:
        temporal_output_process.AddOutput(results_output)
    temporal_output_process._initial_output.AddOutput(file_utils.DeleteOldH5Files())

    # in case the h5py-module is not installed (e.g. on clusters) we don't want it to crash the simulation!
    # => in such a case the xdmf can be created manually afterwards locally
    try:
        from KratosMultiphysics.HDF5Application.xdmf_io import XdmfOutput
        temporal_output_process.AddOutput(XdmfOutput()) # xdmf should be the last in the list
    except ImportError:
        warn_msg  = "XDMF-Writing is not available,\nOnly HDF5 files are written"
        KM.Logger.PrintWarning("SingleMeshXdmfOutputProcess", warn_msg)

    return temporal_output_process
