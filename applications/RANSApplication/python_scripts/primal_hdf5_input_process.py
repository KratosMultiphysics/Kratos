import os

import KratosMultiphysics as Kratos
import KratosMultiphysics
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
            "folder_name"      : ".",
            "time_format"      : "0.6f",
            "echo_level"       : 0
        }
    }""")

    settings.ValidateAndAssignDefaults(default_settings)
    model_part: Kratos.ModelPart = model[settings["model_part_name"].GetString()]

    is_steady = model_part.ProcessInfo[KratosRANS.RANS_IS_STEADY]
    file_settings = settings["file_settings"]
    file_settings.AddEmptyValue("file_name")
    if (is_steady):
        file_settings["file_name"].SetString(os.path.join(file_settings["folder_name"].GetString(), "<model_part_name>-final.h5"))
    else:
        file_settings["file_name"].SetString(os.path.join(file_settings["folder_name"].GetString(), "<model_part_name>-<time>.h5"))
    file_settings.RemoveValue("folder_name")

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
        "list_of_variables" : []
    }""")
    for variable_name in model_part.GetHistoricalVariablesNames():
        if variable_name.find("ADJOINT") == -1 and variable_name.find("SENSITIVITY") == -1:
            list_of_solution_step_variables_parameters["list_of_variables"].Append(variable_name)

    list_of_nodal_variables = ParametersWrapper("""{
        "list_of_variables" : ["RELAXED_ACCELERATION"]
    }""")
    list_of_nodal_flags = ParametersWrapper("""{
        "list_of_variables" : ["SLIP", "INLET", "OUTLET"]
    }""")
    list_of_condition_data_variables = ParametersWrapper("""{
        "list_of_variables" : ["RANS_IS_WALL_FUNCTION_ACTIVE"]
    }""")
    list_of_condition_flag_variables = ParametersWrapper("""{
        "list_of_variables" : ["SLIP"]
    }""")
    core_settings[0]["list_of_operations"] = [
        CreateOperationSettings("nodal_flag_value_input", list_of_nodal_flags),
        CreateOperationSettings("condition_data_value_input", list_of_condition_data_variables),
        CreateOperationSettings("condition_flag_value_input", list_of_condition_flag_variables)
    ]

    core_settings[1]["list_of_operations"] = [
        CreateOperationSettings("nodal_solution_step_data_input", list_of_solution_step_variables_parameters),
        CreateOperationSettings("nodal_data_value_input", list_of_nodal_variables)
    ]

    return core.Factory(core_settings, model, KratosMultiphysics.Process)