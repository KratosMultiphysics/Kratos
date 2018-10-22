import KratosMultiphysics
import KratosMultiphysics.HDF5Application as KratosHDF5
import hdf5_io

def Factory(settings, Model):
    """Return a process for initializing a model part from an existing HDF5 ouput file."""
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    default_settings = KratosMultiphysics.Parameters(r'''{
        "model_part_name" : "MainModelPart",
        "file_settings" : {
            "file_access_mode" : "read_only"
        },
        "nodal_solution_step_data_settings" : {},
        "element_data_value_settings" : {},
        "nodal_data_value_settings" : {},
        "process_info_results_settings" : {},
        "file_name": ""
    }''')

    settings = settings["Parameters"].Clone()
    settings.ValidateAndAssignDefaults(default_settings)
    model_part = Model[settings["model_part_name"].GetString()]
    hdf5_file_factory = hdf5_io.HDF5SerialFileFactory(settings["file_settings"])
    nodal_history_results_input = hdf5_io.NodalSolutionStepDataInput(settings["nodal_solution_step_data_settings"])
    element_data_results_input = hdf5_io.ElementDataValueInput(settings["element_data_value_settings"])
    nodal_data_results_input = hdf5_io.NodalDataValueInput(settings["nodal_data_value_settings"])

    initialization_settings = KratosMultiphysics.Parameters()
    initialization_settings.AddEmptyValue("file_name")
    initialization_settings["file_name"].SetString(settings["file_name"].GetString())

    initialization_process = hdf5_io.InitializationFromInputProcess(model_part,hdf5_file_factory,initialization_settings)
    initialization_process.AddInput(nodal_history_results_input)
    initialization_process.AddInput(element_data_results_input)
    initialization_process.AddInput(nodal_data_results_input)
    return initialization_process

