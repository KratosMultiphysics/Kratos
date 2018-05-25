import KratosMultiphysics
import KratosMultiphysics.HDF5Application as KratosHDF5
import hdf5_io

def Factory(settings, Model):
    """Return a process for writing a transient primal solution to HDF5."""
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    default_settings = KratosMultiphysics.Parameters("""
            {
                "model_part_name" : "MainModelPart",
                "file_settings" : {},
                "nodal_results_settings" : {},
                "element_results_settings" : {},
                "time_tag_precision" : 4,
                "file_name": "DEFAULT_NAME"
            }
            """)
    settings = settings["Parameters"].Clone()
    settings.ValidateAndAssignDefaults(default_settings)
    model_part = Model[settings["model_part_name"].GetString()]
    hdf5_file_factory = hdf5_io.HDF5ParallelFileFactory(settings["file_settings"])
    nodal_results_input = hdf5_io.PrimalBossakInput(settings["nodal_results_settings"])
    element_results_input = hdf5_io.ElementResultsInput(settings["element_results_settings"])
    input_time_settings = KratosMultiphysics.Parameters("""{}""")
    input_time_settings.AddEmptyValue("time_tag_precision")
    input_time_settings["time_tag_precision"].SetInt(settings["time_tag_precision"].GetInt())
    input_time_settings.AddEmptyValue("file_name")
    if settings["file_name"].GetString()=="DEFAULT_NAME":
        input_time_settings["file_name"].SetString(model_part.Name)
    else:
        input_time_settings["file_name"].SetString(settings["file_name"].GetString())
    temporal_input_process = hdf5_io.TemporalInputProcess(
        model_part, hdf5_file_factory, input_time_settings)
    temporal_input_process.AddInput(nodal_results_input)
    temporal_input_process.AddInput(element_results_input)
    return temporal_input_process
