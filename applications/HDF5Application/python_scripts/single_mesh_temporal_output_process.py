import temporal_output_process_factory


def Factory(settings, Model):
    """Return a process for writing simulation results for a single mesh to HDF5."""
    factory_helper = temporal_output_process_factory.TemporalOutputFactoryHelper()
    (temporal_output_process, model_part_output, list_of_results_output) = factory_helper.Execute(settings["Parameters"], Model)
    for results_output in list_of_results_output:
        temporal_output_process.AddOutput(results_output)
    return temporal_output_process
