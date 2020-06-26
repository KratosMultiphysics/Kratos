"""Store temporal simulation results for a single mesh with HDF5.

This process:
 - stores the initial model part in an .h5 file.
 - stores historical and non-historical results in one .h5 file per output step.

This process works with or without MPI.

license: HDF5Application/license.txt
"""


__all__ = ["Factory"]


import KratosMultiphysics
import KratosMultiphysics.HDF5Application.core as core
from KratosMultiphysics.HDF5Application.utils import ParametersWrapper
from KratosMultiphysics.HDF5Application.utils import IsDistributed
from KratosMultiphysics.HDF5Application.utils import CreateOperationSettings


def Factory(settings, Model):
    """Return the process for single mesh temporal output with HDF5.

    The input settings are given in the following table:
    +-------------------------------------+------------+------------------------------------+
    | Setting                             | Type       | Default Value                      |
    +-------------------------------------+------------+------------------------------------+
    | "model_part_name"                   | String     | "MainModelPart"                    |
    +-------------------------------------+------------+------------------------------------+
    | "file_settings"                     | Parameters | "file_name": "<model_part_name>.h5"|
    |                                     |            | "time_format": "0.4f"              |
    |                                     |            | "file_access_mode": "exclusive"    |
    |                                     |            | "echo_level":  0                   |
    +-------------------------------------+------------+------------------------------------+
    | "output_time_settings"              | Parameters | "time_frequency": 1.0              |
    |                                     |            | "step_frequency": 1                |
    +-------------------------------------+------------+------------------------------------+
    | "model_part_output_settings"        | Parameters | "prefix": "/ModelData"             |
    +-------------------------------------+------------+------------------------------------+
    | "nodal_solution_step_data_settings" | Parameters | "prefix": "/ResultsData"           |
    |                                     |            | "list_of_variables": []            |
    +-------------------------------------+------------+------------------------------------+
    | "nodal_data_value_settings"         | Parameters | "prefix": "/ResultsData"           |
    |                                     |            | "list_of_variables": []            |
    +-------------------------------------+------------+------------------------------------+
    | "element_data_value_settings"       | Parameters | "prefix": "/ResultsData"           |
    |                                     |            | "list_of_variables": []            |
    +-------------------------------------+------------+------------------------------------+
    | "nodal_flag_value_settings"         | Parameters | "prefix": "/ResultsData"           |
    |                                     |            | "list_of_variables": []            |
    +-------------------------------------+------------+------------------------------------+
    | "element_flag_value_settings"       | Parameters | "prefix": "/ResultsData"           |
    |                                     |            | "list_of_variables": []            |
    +-------------------------------------+------------+------------------------------------+
    | "condition_flag_value_settings"     | Parameters | "prefix": "/ResultsData"           |
    |                                     |            | "list_of_variables": []            |
    +-------------------------------------+------------+------------------------------------+
    | "condition_data_value_settings"     | Parameters | "prefix": "/ResultsData"           |
    |                                     |            | "list_of_variables": []            |
    +-------------------------------------+------------+------------------------------------+
    """
    core_settings = CreateCoreSettings(settings["Parameters"])
    return SingleMeshTemporalOutputProcessFactory(core_settings, Model)


def SingleMeshTemporalOutputProcessFactory(core_settings, Model):
    return core.Factory(core_settings, Model)


def CreateCoreSettings(user_settings):
    """Return the core settings.

    The core setting "io_type" cannot be overwritten by the user. It is
    automatically set depending on whether or not MPI is used.
    """
    # Configure the defaults:
    core_settings = ParametersWrapper("""
        [{
            "model_part_name": "",
            "process_step": "before_solution_loop",
            "io_settings": {
                "io_type": "serial_hdf5_file_io",
                "file_name": "<model_part_name>.h5"
            },
            "list_of_operations": []
        },{
            "model_part_name" : "",
            "process_step": "finalize_solution_step",
            "controller_settings": {
                "controller_type": "temporal_controller"
            },
            "io_settings": {
                "io_type": "serial_hdf5_file_io",
                "file_name": "<model_part_name>-<time>.h5"
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
        "output_time_settings" : {},
        "model_part_output_settings" : {},
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
    for i in core_settings:
        core_settings[i]["model_part_name"] = user_settings["model_part_name"]
        for key in user_settings["file_settings"]:
            core_settings[i]["io_settings"][key] = user_settings["file_settings"][key]
        if IsDistributed():
            model_part_output_type = "partitioned_model_part_output"
            core_settings[i]["io_settings"]["io_type"] = "parallel_hdf5_file_io"
        else:
            model_part_output_type = "model_part_output"
            core_settings[i]["io_settings"]["io_type"] = "serial_hdf5_file_io"
    core_settings[0]["list_of_operations"] = [
        CreateOperationSettings(model_part_output_type,
                                user_settings["model_part_output_settings"]),
        CreateOperationSettings("nodal_solution_step_data_output",
                                user_settings["nodal_solution_step_data_settings"]),
        CreateOperationSettings("nodal_data_value_output",
                                user_settings["nodal_data_value_settings"]),
        CreateOperationSettings("element_data_value_output",
                                user_settings["element_data_value_settings"]),
        CreateOperationSettings("nodal_flag_value_output",
                                user_settings["nodal_flag_value_settings"]),
        CreateOperationSettings("element_flag_value_output",
                                user_settings["element_flag_value_settings"]),
        CreateOperationSettings("condition_flag_value_output",
                                user_settings["condition_flag_value_settings"]),
        CreateOperationSettings("condition_data_value_output",
                                user_settings["condition_data_value_settings"])
    ]
    core_settings[1]["list_of_operations"] = [
        CreateOperationSettings("nodal_solution_step_data_output",
                                user_settings["nodal_solution_step_data_settings"]),
        CreateOperationSettings("nodal_data_value_output",
                                user_settings["nodal_data_value_settings"]),
        CreateOperationSettings("element_data_value_output",
                                user_settings["element_data_value_settings"]),
        CreateOperationSettings("nodal_flag_value_output",
                                user_settings["nodal_flag_value_settings"]),
        CreateOperationSettings("element_flag_value_output",
                                user_settings["element_flag_value_settings"]),
        CreateOperationSettings("condition_flag_value_output",
                                user_settings["condition_flag_value_settings"]),
        CreateOperationSettings("condition_data_value_output",
                                user_settings["condition_data_value_settings"])
    ]
    for key in user_settings["output_time_settings"]:
        core_settings[1]["controller_settings"][key] = user_settings["output_time_settings"][key]
    return core_settings
