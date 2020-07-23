"""Load simulation results in the initialization step from HDF5.

license: HDF5Application/license.txt
"""


__all__ = ["Factory"]


import KratosMultiphysics
import KratosMultiphysics.HDF5Application.core as core
from KratosMultiphysics.HDF5Application.utils import ParametersWrapper
from KratosMultiphysics.HDF5Application.utils import IsDistributed
from KratosMultiphysics.HDF5Application.utils import CreateOperationSettings


def Factory(settings, Model):
    """Return the process for initialization from HDF5.

    The input settings are given in the following table:
    +-------------------------------------+------------+---------------------------------+
    | Setting                             | Type       | Default Value                   |
    +-------------------------------------+------------+---------------------------------+
    | "model_part_name"                   | String     | ""                              |
    +-------------------------------------+------------+---------------------------------+
    | "file_settings"                     | Parameters | "file_name": "<model_part_name>"|
    |                                     |            | "echo_level":  0                |
    +-------------------------------------+------------+---------------------------------+
    | "nodal_solution_step_data_settings" | Parameters | "prefix": "/ResultsData"        |
    |                                     |            | "list_of_variables": []         |
    +-------------------------------------+------------+---------------------------------+
    | "nodal_data_value_settings"         | Parameters | "prefix": "/ResultsData"        |
    |                                     |            | "list_of_variables": []         |
    +-------------------------------------+------------+---------------------------------+
    | "element_data_value_settings"       | Parameters | "prefix": "/ResultsData"        |
    |                                     |            | "list_of_variables": []         |
    +-------------------------------------+------------+---------------------------------+
    | "nodal_flag_value_settings"         | Parameters | "prefix": "/ResultsData"        |
    |                                     |            | "list_of_variables": []         |
    +-------------------------------------+------------+---------------------------------+
    | "element_flag_value_settings"       | Parameters | "prefix": "/ResultsData"        |
    |                                     |            | "list_of_variables": []         |
    +-------------------------------------+------------+---------------------------------+
    | "condition_flag_value_settings"     | Parameters | "prefix": "/ResultsData"        |
    |                                     |            | "list_of_variables": []         |
    +-------------------------------------+------------+---------------------------------+
    | "condition_data_value_settings"     | Parameters | "prefix": "/ResultsData"        |
    |                                     |            | "list_of_variables": []         |
    +-------------------------------------+------------+---------------------------------+
    """
    core_settings = CreateCoreSettings(settings["Parameters"])
    return InitializationFromHDF5ProcessFactory(core_settings, Model)


def InitializationFromHDF5ProcessFactory(core_settings, Model):
    return core.Factory(core_settings, Model)


def CreateCoreSettings(user_settings):
    """Return the core settings."""
    # Configure the defaults:
    core_settings = ParametersWrapper("""
        [{
            "model_part_name" : "",
            "process_step": "initialize",
            "io_settings": {
                "io_type": "serial_hdf5_file_io",
                "file_name": "<model_part_name>.h5",
                "file_access_mode": "read_only"
            },
            "list_of_operations": []
        }]
        """)
    # Apply the user settings:
    user_settings.ValidateAndAssignDefaults(
        KratosMultiphysics.Parameters("""
            {
                "model_part_name" : "MainModelPart",
                "file_settings" : {},
                "nodal_solution_step_data_settings" : {},
                "nodal_data_value_settings": {},
                "element_data_value_settings" : {},
                "nodal_flag_value_settings": {},
                "element_flag_value_settings" : {},
                "condition_data_value_settings" : {},
                "condition_flag_value_settings" : {}
            }
            """))
    user_settings = ParametersWrapper(user_settings)
    core_settings[0]["model_part_name"] = user_settings["model_part_name"]
    for key in user_settings["file_settings"]:
        core_settings[0]["io_settings"][key] = user_settings["file_settings"][key]
    core_settings[0]["file_access_mode"] = "read_only"
    if IsDistributed():
        core_settings[0]["io_settings"]["io_type"] = "parallel_hdf5_file_io"
    else:
        core_settings[0]["io_settings"]["io_type"] = "serial_hdf5_file_io"
    core_settings[0]["list_of_operations"] = [
        CreateOperationSettings("nodal_solution_step_data_input",
                                user_settings["nodal_solution_step_data_settings"]),
        CreateOperationSettings("nodal_data_value_input",
                                user_settings["nodal_data_value_settings"]),
        CreateOperationSettings("element_data_value_input",
                                user_settings["element_data_value_settings"]),
        CreateOperationSettings("nodal_flag_value_input",
                                user_settings["nodal_flag_value_settings"]),
        CreateOperationSettings("element_flag_value_input",
                                user_settings["element_flag_value_settings"]),
        CreateOperationSettings("condition_flag_value_input",
                                user_settings["condition_flag_value_settings"]),
        CreateOperationSettings("condition_data_value_input",
                                user_settings["condition_data_value_settings"])
    ]
    return core_settings
