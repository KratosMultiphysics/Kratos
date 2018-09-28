import temporal_output_process_factory
from xdmf_io import XdmfOutput

def GetFactoryHelper():
    """Return the FactoryHelper to be used for the TemporalOutput"""
    return temporal_output_process_factory.TemporalOutputFactoryHelper()

def IsRankZero():
    return True

def Factory(settings, Model):
    ###### """Return a process for writing simulation results for a single mesh to HDF5."""
    factory_helper = GetFactoryHelper()
    list_aux_initial_outputs = []
    if IsRankZero():
        list_aux_initial_outputs.append(XdmfOutput(settings["Parameters"]["xdmf_settings"]))
        settings["Parameters"].RemoveValue("xdmf_settings")

    (temporal_output_process, model_part_output, list_of_results_output) = factory_helper.Execute(settings["Parameters"], Model)
    list_of_results_output.extend(list_aux_initial_outputs)
    for results_output in list_of_results_output:
        temporal_output_process.AddOutput(results_output)
    # for aux_output in list_aux_initial_outputs:
    #     temporal_output_process._initial_output.AddOutput(aux_output)

    return temporal_output_process
