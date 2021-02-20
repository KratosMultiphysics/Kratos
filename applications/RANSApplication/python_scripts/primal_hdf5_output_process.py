import KratosMultiphysics as Kratos
import KratosMultiphysics.RANSApplication as KratosRANS
import KratosMultiphysics.HDF5Application as KratosHDF5

from KratosMultiphysics import IsDistributedRun

import KratosMultiphysics.HDF5Application.core as core
from KratosMultiphysics.HDF5Application.core import AssignOperationsToController
from KratosMultiphysics.HDF5Application.core import CreateControllerWithFileIO
from KratosMultiphysics.HDF5Application.core import AssignControllerToProcess
from KratosMultiphysics.HDF5Application.core.processes import ControllerProcess
from KratosMultiphysics.HDF5Application.utils import ParametersWrapper
from KratosMultiphysics.HDF5Application.utils import CreateOperationSettings

from KratosMultiphysics.RANSApplication import RansVariableUtilities

def Factory(parameters, model):
    settings = parameters["Parameters"]
    default_settings = Kratos.Parameters("""{
        "model_part_name": "",
        "steady"         : true,
        "file_settings"  : {
            "file_name"        : "./<model_part_name>-<time>.h5",
            "time_format"      : "0.6f",
            "file_access_mode" : "truncate",
            "max_files_to_keep": "unlimited",
            "echo_level"       : 0
        }
    }""")

    settings.ValidateAndAssignDefaults(default_settings)
    model_part = model[settings["model_part_name"].GetString()]

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

        if (not settings["steady"].GetBool()):
            core_settings[i]["io_settings"]["max_files_to_keep"] = "unlimited"

    list_of_solution_step_variables = ParametersWrapper("""{
        "list_of_variables" : ["ALL_VARIABLES"]
    }""")
    list_of_nodal_variables = ParametersWrapper("""{
        "list_of_variables" : ["RELAXED_ACCELERATION"]
    }""")
    list_of_nodal_flags = ParametersWrapper("""{
        "list_of_variables" : ["SLIP", "INLET", "OUTLET"]
    }""")
    list_of_condition_data_variables = ParametersWrapper("""{
        "list_of_variables" : ["RANS_IS_WALL_FUNCTION_ACTIVE", "RANS_IS_INLET"]
    }""")

    core_settings[0]["list_of_operations"] = [
        CreateOperationSettings(model_part_output_type, file_settings),
        CreateOperationSettings("nodal_solution_step_data_output", list_of_solution_step_variables),
        CreateOperationSettings("nodal_data_value_output", list_of_nodal_variables),
        CreateOperationSettings("nodal_flag_value_output", list_of_nodal_flags)
        # CreateOperationSettings("condition_data_value_output", list_of_condition_data_variables)
    ]
    core_settings[1]["list_of_operations"] = [
        CreateOperationSettings("nodal_solution_step_data_output", list_of_solution_step_variables),
        CreateOperationSettings("nodal_data_value_output", list_of_nodal_variables),
        CreateOperationSettings("nodal_flag_value_output", list_of_nodal_flags)
        # CreateOperationSettings("condition_data_value_output", list_of_condition_data_variables)
    ]

    core_settings[1]["controller_settings"]["step_frequency"] = 1

    process = PrimalHDF5OutputProcess()
    for i in core_settings:
        controller = CreateControllerWithFileIO(core_settings[i], model)
        AssignOperationsToController(core_settings[i]['list_of_operations'], controller)
        AssignControllerToProcess(core_settings[i], controller, process)
    return process

class PrimalHDF5OutputProcess(ControllerProcess):
    def IsOutputStep(self):
        return False