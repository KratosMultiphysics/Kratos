'''Write a model part once and its data at the end of each solution step.

This process provides the front end to the HDF5Application.core.

license: HDF5Application/license.txt
'''
import KratosMultiphysics
import KratosMultiphysics.HDF5Application.core as _core
import KratosMultiphysics.HDF5Application.utils as _utils


def Factory(process_settings, Model):
    """Return a process for writing simulation results for a single mesh to HDF5."""
    default_settings = KratosMultiphysics.Parameters("""
            {
                "model_part_name" : "MainModelPart",
                "file_settings" : {},
                "model_part_output_settings" : {},
                "nodal_solution_step_data_settings" : {},
                "nodal_data_value_settings": {},
                "element_data_value_settings" : {},
                "output_time_settings" : {}
            }
            """)
    settings = process_settings["Parameters"]
    settings.ValidateAndAssignDefaults(default_settings)
    new_settings = KratosMultiphysics.Parameters('''
            {
               "list_of_controllers": [{
                    "model_part_name" : "",
                    "process_step": "before_solution_loop",
                    "io_settings": {
                        "io_type": "serial_hdf5_file_io",
                        "file_name": "<identifier>.h5"
                    },
                    "list_of_operations": [{
                        "operation_type": "model_part_output"
                    },{
                        "operation_type": "nodal_solution_step_data_output"
                    },{
                        "operation_type": "nodal_data_value_output"
                    },{
                        "operation_type": "element_data_value_output"
                    }]
                },{
                    "model_part_name" : "",
                    "process_step": "finalize_solution_step",
                    "controller_settings": {
                        "controller_type": "temporal_controller"
                    },
                    "io_settings": {
                        "io_type": "serial_hdf5_file_io",
                        "file_name": "<identifier>-<time>.h5"
                    },
                    "list_of_operations": [{
                        "operation_type": "nodal_solution_step_data_output"
                    },{
                        "operation_type": "nodal_data_value_output"
                    },{
                        "operation_type": "element_data_value_output"
                    }]
                }]
            }
            ''')
    model_part_name = settings["model_part_name"].GetString()
    for current_settings in new_settings["list_of_controllers"]:
        current_settings["model_part_name"].SetString(model_part_name)
    model_part_settings = new_settings["list_of_controllers"][0]
    results_settings = new_settings["list_of_controllers"][1]
    for io_settings in [model_part_settings["io_settings"], results_settings["io_settings"]]:
        _utils.InsertSettings(settings["file_settings"], io_settings)
        if _utils.IsDistributed():
            io_settings["io_type"].SetString("parallel_hdf5_file_io")
    _utils.InsertArrayOfSettings(
        [settings["model_part_output_settings"], settings["nodal_solution_step_data_settings"], settings["nodal_data_value_settings"],
         settings["element_data_value_settings"]], model_part_settings["list_of_operations"])
    _utils.InsertArrayOfSettings([settings["nodal_solution_step_data_settings"], settings["nodal_data_value_settings"],
                                  settings["element_data_value_settings"]], results_settings["list_of_operations"])
    output_time_settings = settings["output_time_settings"]
    _utils.InsertSettings(output_time_settings,
                          results_settings["controller_settings"])
    _utils.CheckForDeprecatedFilename(
        output_time_settings, __name__, model_part_settings["io_settings"], results_settings["io_settings"])
    _utils.CheckForDeprecatedTemporalSettings(
        output_time_settings, __name__, results_settings["controller_settings"])

    return _core.Factory(new_settings["list_of_controllers"], Model)
