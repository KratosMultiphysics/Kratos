import os

import KratosMultiphysics as Kratos
import KratosMultiphysics.RANSApplication as KratosRANS
import KratosMultiphysics.HDF5Application as KratosHDF5

from KratosMultiphysics import IsDistributedRun

import KratosMultiphysics.HDF5Application.core as core
from KratosMultiphysics.HDF5Application.core import AssignOperationsToController
from KratosMultiphysics.HDF5Application.core import CreateControllerWithFileIO
from KratosMultiphysics.HDF5Application.core import AssignControllerToProcess
from KratosMultiphysics.HDF5Application.core.processes import OrderedOperationProcess
from KratosMultiphysics.HDF5Application.utils import ParametersWrapper
from KratosMultiphysics.HDF5Application.utils import CreateOperationSettings

from KratosMultiphysics.RANSApplication import RansVariableUtilities

def Factory(parameters, model):
    settings = parameters["Parameters"]
    default_settings = Kratos.Parameters("""{
        "model_part_name": "",
        "model_part_output_settings":{},
        "file_settings"  : {
            "folder_name"      : ".",
            "time_format"      : "0.6f",
            "file_access_mode" : "truncate",
            "max_files_to_keep": "unlimited",
            "use_steps"        : false,
            "echo_level"       : 0
        }
    }""")

    settings.ValidateAndAssignDefaults(default_settings)
    model_part = model[settings["model_part_name"].GetString()]

    is_steady = model_part.ProcessInfo[KratosRANS.RANS_IS_STEADY]
    file_settings = settings["file_settings"]
    file_settings.AddEmptyValue("file_name")
    if (is_steady):
        file_settings["file_name"].SetString(os.path.join(file_settings["folder_name"].GetString(), "<model_part_name>-final.h5"))
    else:
        if (file_settings.Has("use_steps") and file_settings["use_steps"].GetBool()):
            file_settings["file_name"].SetString(os.path.join(file_settings["folder_name"].GetString(), "<model_part_name>-<step>.h5"))
        else:
            file_settings["file_name"].SetString(os.path.join(file_settings["folder_name"].GetString(), "<model_part_name>-<time>.h5"))
    file_settings.RemoveValue("folder_name")

    core_settings = ParametersWrapper("""
        [{
            "model_part_name" : "",
            "process_step": "before_solution_loop",
            "io_settings": {
                "io_type": "serial_hdf5_file_io",
                "file_name": "<model_part_name>-<time>.h5"
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
        },{
            "model_part_name" : "",
            "process_step": "finalize",
            "io_settings": {
                "io_type": "serial_hdf5_file_io",
                "file_name": "<model_part_name>-<time>.h5"
            },
            "list_of_operations": []
        }        ]
        """)

    # set file settings
    file_settings = ParametersWrapper(settings["file_settings"])
    for i in core_settings:
        core_settings[i]["model_part_name"] = model_part.Name
        for key in file_settings:
            core_settings[i]["io_settings"][key] = file_settings[key]
        if IsDistributedRun():
            model_part_output_type = "partitioned_model_part_output"
            core_settings[i]["io_settings"]["io_type"] = "parallel_hdf5_file_io"
        else:
            model_part_output_type = "model_part_output"
            core_settings[i]["io_settings"]["io_type"] = "serial_hdf5_file_io"

    list_of_solution_step_variables = ParametersWrapper("""{
        "list_of_variables" : ["ALL_VARIABLES_FROM_VARIABLES_LIST"]
    }""")
    list_of_nodal_variables = ParametersWrapper("""{
        "list_of_variables" : ["RELAXED_ACCELERATION"]
    }""")
    list_of_nodal_flags = ParametersWrapper("""{
        "list_of_variables" : ["SLIP", "INLET", "OUTLET"]
    }""")
    list_of_condition_data_variables = ParametersWrapper("""{
        "list_of_variables" : ["RANS_IS_WALL_FUNCTION_ACTIVE"]
    }""")
    list_of_condition_flags = ParametersWrapper("""{
        "list_of_variables" : ["SLIP"]
    }""")

    list_of_element_variables = ParametersWrapper("""{
        "list_of_variables" : ["ADJOINT_STABILIZATION_COEFFICIENT", "ELEMENT_H", "ELEMENT_ERROR"]
    }""")

    initialize_list = []
    finalize_list = []
    if (is_steady):
        initialize_list = []
        finalize_list = [
            CreateOperationSettings(model_part_output_type, ParametersWrapper(settings["model_part_output_settings"]))
        ]
    else:
        initialize_list = [
            CreateOperationSettings(model_part_output_type, ParametersWrapper(settings["model_part_output_settings"]))
        ]
        finalize_list = []

    initialize_list.extend([
        CreateOperationSettings("nodal_solution_step_data_output", list_of_solution_step_variables),
        CreateOperationSettings("nodal_data_value_output", list_of_nodal_variables),
        CreateOperationSettings("nodal_flag_value_output", list_of_nodal_flags),
        CreateOperationSettings("condition_data_value_output", list_of_condition_data_variables),
        CreateOperationSettings("condition_flag_value_output", list_of_condition_flags)
    ])

    transient_list = [
        CreateOperationSettings(model_part_output_type, ParametersWrapper(settings["model_part_output_settings"])),
        CreateOperationSettings("nodal_solution_step_data_output", list_of_solution_step_variables),
        CreateOperationSettings("nodal_data_value_output", list_of_nodal_variables),
        CreateOperationSettings("nodal_flag_value_output", list_of_nodal_flags),
        CreateOperationSettings("condition_data_value_output", list_of_condition_data_variables),
        CreateOperationSettings("condition_flag_value_output", list_of_condition_flags)
    ]

    finalize_list.extend([
        CreateOperationSettings("nodal_solution_step_data_output", list_of_solution_step_variables),
        CreateOperationSettings("nodal_data_value_output", list_of_nodal_variables),
        CreateOperationSettings("nodal_flag_value_output", list_of_nodal_flags),
        CreateOperationSettings("condition_data_value_output", list_of_condition_data_variables),
        CreateOperationSettings("condition_flag_value_output", list_of_condition_flags)
    ])

    if (not is_steady):
        transient_list.extend([
            CreateOperationSettings(
                "element_data_value_output", list_of_element_variables)
        ])
        initialize_list.extend([
            CreateOperationSettings(
                "element_data_value_output", list_of_element_variables)
        ])
        finalize_list.extend([
            CreateOperationSettings(
                "element_data_value_output", list_of_element_variables)
        ])

    core_settings[0]["list_of_operations"] = initialize_list
    core_settings[1]["list_of_operations"] = transient_list
    core_settings[2]["list_of_operations"] = finalize_list

    core_settings[1]["controller_settings"]["step_frequency"] = 1

    process = PrimalHDF5OutputProcess()
    for i in core_settings:
        controller = CreateControllerWithFileIO(core_settings[i], model)
        AssignOperationsToController(core_settings[i]['list_of_operations'], controller)
        AssignControllerToProcess(core_settings[i], controller, process)

    Kratos.Logger.PrintInfo("PrimalHDF5OutputProcess", "PrimalHDF5OutputProcess created.")
    return process

class PrimalHDF5OutputProcess(OrderedOperationProcess):
    def IsOutputStep(self):
        return False
