"""Store temporal simulation results for a single mesh with HDF5 and Xdmf.

This process:
 - removes existing .h5 files at the start of the simulation.
 - stores the initial model part in an .h5 file.
 - stores historical and non-historical results in one .h5 file per output step.
 - stores Xdmf metadata for post-processing (e.g., Paraview or VisIt).

This process works with or without MPI.

license: HDF5Application/license.txt

Main authors:
    Philipp Bucher
    Michael Andre
"""


__all__ = ["Factory"]


import KratosMultiphysics
from KratosMultiphysics.HDF5Application import core
from KratosMultiphysics.HDF5Application.utils import ParametersWrapper
from KratosMultiphysics.HDF5Application.utils import IsDistributed
from KratosMultiphysics.HDF5Application.utils import CreateOperationSettings

def Factory(settings, Model):
    """Return a process for single mesh temporal results output with Xdmf and HDF5.

    The input settings are given in the following table:
    +-------------------------------------+------------+---------------------------------+
    | Setting                             | Type       | Default Value                   |
    +-------------------------------------+------------+---------------------------------+
    | "model_part_name"                   | String     | ""                              |
    +-------------------------------------+------------+---------------------------------+
    | "file_settings"                     | Parameters | "file_name": "<model_part_name>"|
    |                                     |            | "time_format": "0.4f"           |
    |                                     |            | "file_access_mode": "truncate"  |
    |                                     |            | "echo_level":  0                |
    +-------------------------------------+------------+---------------------------------+
    | "output_time_settings"              | Parameters | "time_frequency": 1.0           |
    |                                     |            | "step_frequency": 1             |
    +-------------------------------------+------------+---------------------------------+
    | "model_part_output_settings"        | Parameters | "prefix": "/ModelData"          |
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
    return SingleMeshXdmfOutputProcessFactory(core_settings, Model)


def SingleMeshXdmfOutputProcessFactory(core_settings, Model):
    return core.Factory(core_settings, Model)


def CreateCoreSettings(user_settings):
    """Return the core settings.

    The core setting "io_type" cannot be overwritten by the user. It is
    automatically set depending on whether or not MPI is used.
    """
    # Configure the defaults:
    core_settings = ParametersWrapper("""
        [{
            "model_part_name" : "",
            "process_step": "initialize",
            "io_settings": {
                "io_type": "mock_hdf5_file_io",
                "file_name": "<model_part_name>.h5"
            },
            "list_of_operations": [{
                "module_name": "operations.system",
                "operation_type": "delete_old_h5_files"
            }]
        },{
            "model_part_name": "",
            "process_step": "before_solution_loop",
            "io_settings": {
                "io_type": "serial_hdf5_file_io",
                "file_name": "<model_part_name>.h5",
                "file_access_mode": "truncate"
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
                "file_name": "<model_part_name>-<time>.h5",
                "file_access_mode": "truncate"
            },
            "list_of_operations": []
        },{
            "model_part_name" : "",
            "process_step": "finalize_solution_step",
            "controller_settings": {
                "controller_type": "temporal_controller"
            },
            "io_settings": {
                "io_type": "mock_hdf5_file_io",
                "file_name": "<model_part_name>.h5"
            },
            "list_of_operations": [{
                "module_name": "operations.xdmf",
                "operation_type": "xdmf_output"
            }]
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
    core_settings[0]["io_settings"]["io_type"] = "mock_hdf5_file_io"
    core_settings[3]["io_settings"]["io_type"] = "mock_hdf5_file_io"
    if IsDistributed():
        model_part_output_type = "partitioned_model_part_output"
        core_settings[1]["io_settings"]["io_type"] = "parallel_hdf5_file_io"
        core_settings[2]["io_settings"]["io_type"] = "parallel_hdf5_file_io"
    else:
        model_part_output_type = "model_part_output"
        core_settings[1]["io_settings"]["io_type"] = "serial_hdf5_file_io"
        core_settings[2]["io_settings"]["io_type"] = "serial_hdf5_file_io"
    core_settings[1]["list_of_operations"] = [
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
    core_settings[2]["list_of_operations"] = [
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
        core_settings[2]["controller_settings"][key] = user_settings["output_time_settings"][key]
        core_settings[3]["controller_settings"][key] = user_settings["output_time_settings"][key]
    return core_settings
