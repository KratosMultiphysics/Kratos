'''Read model part data at the start of each solution step.

This process provides the front end to the HDF5Application.core.

license: HDF5Application/license.txt
'''
import KratosMultiphysics
import KratosMultiphysics.HDF5Application.core as _core
import KratosMultiphysics.HDF5Application.utils as _utils


def Factory(process_settings, Model):
    """Return a process to read a transient solution from HDF5."""
    default_settings = KratosMultiphysics.Parameters("""
            {
                "model_part_name" : "MainModelPart",
                "file_settings" : {},
                "nodal_solution_step_data_settings" : {},
                "nodal_data_value_settings": {},
                "element_data_value_settings" : {}
            }
            """)
    new_settings = KratosMultiphysics.Parameters('''
            {
               "list_of_controllers": [{
                    "model_part_name" : "",
                    "process_step": "initialize_solution_step",
                    "io_settings": {
                        "io_type": "serial_hdf5_file_io",
                        "file_name": "<identifier>-<time>.h5",
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
    settings = process_settings["Parameters"]
    if settings.Has('file_name'):
        _utils.CheckForDeprecatedFilename(
            settings, __name__, new_settings["list_of_controllers"][0]["io_settings"])
        settings.RemoveValue('file_name')
    if settings.Has('time_tag_precision'):
        depr_msg = '\nDEPRECATION-WARNING: "time_tag_precision" is ignored, please use "time_format" in "file_settings"!\n'
        KratosMultiphysics.Logger.PrintWarning(__name__, depr_msg)
        settings.RemoveValue('time_tag_precision')
    settings.ValidateAndAssignDefaults(default_settings)
    model_part_name = settings["model_part_name"].GetString()
    results_settings = new_settings["list_of_controllers"][0]
    results_settings["model_part_name"].SetString(model_part_name)
    _utils.InsertSettings(
        settings["file_settings"], results_settings["io_settings"])
    if _utils.IsDistributed():
        results_settings["io_settings"]["io_type"].SetString(
            "parallel_hdf5_file_io")
    _utils.InsertArrayOfSettings([settings["nodal_solution_step_data_settings"], settings["nodal_data_value_settings"],
                                  settings["element_data_value_settings"]], results_settings["list_of_operations"])
    return _core.Factory(new_settings["list_of_controllers"], Model)
