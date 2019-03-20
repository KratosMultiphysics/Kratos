import KratosMultiphysics
import KratosMultiphysics.HDF5Application as KratosHDF5
import hdf5_io

def Factory(settings, Model):
    """Return a process to read a transient solution from HDF5."""
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    default_settings = KratosMultiphysics.Parameters("""
            {
                "model_part_name" : "MainModelPart",
                "file_settings" : {},
                "nodal_solution_step_data_settings" : {},
                "element_data_value_settings" : {},
                "nodal_data_value_settings" : {},
                "time_tag_precision" : 4,
                "file_name": "DEFAULT_NAME"
            }
            """)
    settings = settings["Parameters"].Clone()
    settings.ValidateAndAssignDefaults(default_settings)
    model_part = Model[settings["model_part_name"].GetString()]
    hdf5_file_factory = hdf5_io.HDF5ParallelFileFactory(settings["file_settings"])
    nodal_solution_step_input = hdf5_io.NodalSolutionStepDataInput(settings["nodal_solution_step_data_settings"])
    element_data_value_input = hdf5_io.ElementDataValueInput(settings["element_data_value_settings"])
    nodal_data_value_input = hdf5_io.NodalDataValueInput(settings["nodal_data_value_settings"])
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
    temporal_input_process.AddInput(nodal_solution_step_input)
    temporal_input_process.AddInput(element_data_value_input)
    temporal_input_process.AddInput(nodal_data_value_input)
    return temporal_input_process
