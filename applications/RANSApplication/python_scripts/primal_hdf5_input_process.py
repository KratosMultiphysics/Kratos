import KratosMultiphysics as Kratos
import KratosMultiphysics.RANSApplication as KratosRANS
import KratosMultiphysics.HDF5Application as KratosHDF5

from KratosMultiphysics import IsDistributedRun

import KratosMultiphysics.HDF5Application.core as core
from KratosMultiphysics.HDF5Application.utils import ParametersWrapper
from KratosMultiphysics.HDF5Application.utils import CreateOperationSettings

from KratosMultiphysics.RANSApplication import RansVariableUtilities

def Factory(parameters, model):
    settings = parameters["Parameters"]
    default_settings = Kratos.Parameters("""{
        "model_part_name": "",
        "file_settings"  : {
            "file_name"        : "./<model_part_name>-<time>.h5",
            "time_format"      : "0.6f",
            "echo_level"       : 0
        }
    }""")

    settings.ValidateAndAssignDefaults(default_settings)
    model_part = model[settings["model_part_name"].GetString()]

    core_settings = ParametersWrapper("""
        [
        {
            "model_part_name" : "",
            "process_step": "initialize",
            "io_settings": {
                "io_type": "serial_hdf5_file_io",
                "file_name": "<model_part_name>-<time>.h5",
                "file_access_mode": "read_only"
            },
            "list_of_operations": []
        },
        {
            "model_part_name" : "",
            "process_step": "initialize_solution_step",
            "io_settings": {
                "io_type": "serial_hdf5_file_io",
                "file_name": "<model_part_name>-<time>.h5",
                "file_access_mode": "read_only"
            },
            "list_of_operations": []
        }
        ]
        """)

    # set file settings
    file_settings = ParametersWrapper(settings["file_settings"])
    for i in core_settings:
        core_settings[i]["model_part_name"] = model_part.Name
        for key in file_settings:
            core_settings[i]["io_settings"][key] = file_settings[key]
        if IsDistributedRun():
            core_settings[i]["io_settings"]["io_type"] = "parallel_hdf5_file_io"
        else:
            core_settings[i]["io_settings"]["io_type"] = "serial_hdf5_file_io"

    list_of_solution_step_variables_parameters = ParametersWrapper("""{
        "list_of_variables" : ["ALL_VARIABLES"]
    }""")

    list_of_nodal_variables = ParametersWrapper("""{
        "list_of_variables" : ["ALL_VARIABLES"]
    }""")
    list_of_nodal_flags = ParametersWrapper("""{
        "list_of_variables" : ["SLIP", "INLET", "OUTLET"]
    }""")
    list_of_condition_data_variables = ParametersWrapper("""{
        "list_of_variables" : ["ALL_VARIABLES"]
    }""")

    core_settings[0]["list_of_operations"] = [
        CreateOperationSettings("nodal_flag_value_input", list_of_nodal_flags)
        # CreateOperationSettings("condition_data_value_input", list_of_condition_data_variables)
    ]

    core_settings[1]["list_of_operations"] = [
        CreateOperationSettings("nodal_flag_value_input", list_of_nodal_flags),
        CreateOperationSettings("nodal_solution_step_data_input", list_of_solution_step_variables_parameters),
        CreateOperationSettings("nodal_data_value_input", list_of_nodal_variables)
    ]

    return core.Factory(core_settings, model)