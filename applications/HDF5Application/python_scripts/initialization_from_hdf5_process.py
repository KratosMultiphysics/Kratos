'''Read model part data in the initialization step.

This process provides the front end to the HDF5Application.core.

license: HDF5Application/license.txt
'''
import KratosMultiphysics
import KratosMultiphysics.HDF5Application.core as _core
import KratosMultiphysics.HDF5Application.utils as _utils


def Factory(process_settings, Model):
    """Return a process for initializing a model part from an existing HDF5 output file."""
    default_settings = KratosMultiphysics.Parameters("""
            {
                "model_part_name" : "MainModelPart",
                "file_settings" : {},
                "nodal_solution_step_data_settings" : {},
                "nodal_data_value_settings": {},
                "element_data_value_settings" : {}
            }
            """)
    settings = process_settings["Parameters"]
    settings.ValidateAndAssignDefaults(default_settings)
    new_settings = KratosMultiphysics.Parameters('''
            {
               "list_of_controllers": [{
                    "model_part_name" : "",
                    "process_step": "initialize",
                    "io_settings": {
                        "io_type": "serial_hdf5_file_io",
                        "file_name": "<identifier>.h5",
                        "file_access_mode": "read_only"
                    },
                    "list_of_operations": [{
                        "operation_type": "nodal_solution_step_data_input"
                    },{
                        "operation_type": "nodal_data_value_input"
                    },{
                        "operation_type": "element_data_value_input"
                    }]
                }]
            }
            ''')
    model_part_name = settings["model_part_name"].GetString()
    new_settings["list_of_controllers"][0]["model_part_name"].SetString(
        model_part_name)
    results_settings = new_settings["list_of_controllers"][0]
    _utils.InsertSettings(
        settings["file_settings"], results_settings["io_settings"])
    if _utils.IsDistributed():
        results_settings["io_settings"]["io_type"].SetString(
            "parallel_hdf5_file_io")
    _utils.InsertArrayOfSettings([settings["nodal_solution_step_data_settings"], settings["nodal_data_value_settings"],
                                  settings["element_data_value_settings"]], results_settings["list_of_operations"])
    return _core.Factory(new_settings["list_of_controllers"], Model)
