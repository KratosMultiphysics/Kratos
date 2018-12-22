import KratosMultiphysics as KM
import temporal_output_process_factory
import file_utilities

def Factory(settings, Model):
    """Return a process for writing simulation results for a single mesh to HDF5.
    It also creates the xdmf-files on the fly
    """
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")

    params = settings["Parameters"]
    model_part_name = params["model_part_name"].GetString() # name of modelpart must be specified!
 #todo(msandre): collapse older partitioned scripts to their serial counterparts like this
    if Model[model_part_name].GetCommunicator().TotalProcesses() > 1:
        factory_helper = temporal_output_process_factory.PartitionedTemporalOutputFactoryHelper()
    else:
        factory_helper = temporal_output_process_factory.TemporalOutputFactoryHelper()

    (temporal_output_process, model_part_output, list_of_results_output) = factory_helper.Execute(params, Model)
    for results_output in list_of_results_output:
        temporal_output_process.AddOutput(results_output)

    temporal_output_process._initial_output.AddOutput(file_utilities.DeleteOldH5Files())

    # in case the h5py-module is not installed (e.g. on clusters) we don't want it to crash the simulation!
    # => in such a case the xdmf can be created manually afterwards locally
    try:
        from xdmf_io import XdmfOutput
        temporal_output_process.AddOutput(XdmfOutput()) # xdmf should be the last in the list
    except ImportError:
        KM.Logger.PrintWarning("XDMF-Writing is not available!")

    return temporal_output_process
