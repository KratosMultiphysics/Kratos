from single_mesh_temporal_output_process import FactoryHelper

def Factory(settings, Model):
    """Return a process for writing simulation results for multiple meshes to HDF5."""
    (temporal_output_process, model_part_output, list_of_results_output) = FactoryHelper(settings, Model)
    temporal_output_process.AddOutput(model_part_output)
    for results_output in list_of_results_output:
        temporal_output_process.AddOutput(results_output)
    return temporal_output_process
