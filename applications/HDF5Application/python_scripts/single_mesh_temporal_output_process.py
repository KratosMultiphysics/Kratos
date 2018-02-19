import KratosMultiphysics
import KratosMultiphysics.HDF5Application as KratosHDF5
import hdf5_output

def FactoryHelper(settings, Model):
    """Return objects needed for constructing a temporal output process."""
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    model_part = Model[settings["model_part_name"].GetString()]
    hdf5_file_factory = hdf5_output.HDF5SerialFileFactory(settings["file_output_settings"])
    model_part_output = hdf5_output.ModelPartOutput(settings["model_part_output_settings"])
    results_output = hdf5_output.NodalResultsOutput(settings["results_settings"])
    temporal_output_process = hdf5_output.TemporalOutputProcess(
            model_part, hdf5_file_factory, settings["output_time_settings"], [model_part_output, results_output])
    return (temporal_output_process, model_part_output, [results_output])


def Factory(settings, Model):
    """Return a process for writing simulation results for a single mesh to HDF5."""
    (temporal_output_process, model_part_output, list_of_results_output) = FactoryHelper(settings, Model)
    for results_output in list_of_results_output:
        temporal_output_process.AddOutput(results_output)
    return temporal_output_process
