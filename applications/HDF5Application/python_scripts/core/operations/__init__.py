'''HDF5 IO operations.

Exports operations from:
    - core.operations.model_part

license: HDF5Application/license.txt
'''

__all__ = ["CreateAggregatedOperation"]

import KratosMultiphysics as Kratos

from .model_part import *
from .xdmf import *
from .system import *
from .aggregated_operations import *

from KratosMultiphysics.HDF5Application.core.controllers import Factory as ControllerFactory

operation_types_map: 'dict[str, typing.Any]' = {
    "model_part_input": ModelPartInput,
    "model_part_output": ModelPartOutput,
    "partitioned_model_part_output": PartitionedModelPartOutput,
    "process_info_input": ProcessInfoInput,
    "process_info_output": ProcessInfoOutput,
    "nodal_solution_step_data_input" : NodalSolutionStepDataInput,
    "nodal_solution_step_data_output" : NodalSolutionStepDataOutput,
    "primal_bossak_input": PrimalBossakInput,
    "primal_bossak_output": PrimalBossakOutput,
    "nodal_data_value_input" : NodalDataValueInput,
    "nodal_data_value_output" : NodalDataValueOutput,
    "nodal_flag_value_input" : NodalFlagValueInput,
    "nodal_flag_value_output" : NodalFlagValueOutput,
    "element_data_value_input" : ElementDataValueInput,
    "element_data_value_output" : ElementDataValueOutput,
    "element_flag_value_input" : ElementFlagValueInput,
    "element_flag_value_output" : ElementFlagValueOutput,
    "element_gauss_point_value_output" : ElementGaussPointOutput,
    "condition_data_value_input" : ConditionDataValueInput,
    "condition_data_value_output" : ConditionDataValueOutput,
    "condition_flag_value_input" : ConditionFlagValueInput,
    "condition_flag_value_output" : ConditionFlagValueOutput,
    "condition_gauss_point_value_output" : ConditionGaussPointOutput,
    "xdmf_output": XdmfOutput,
    "delete_old_h5_file": DeleteOldH5Files
}

def CreateAggregatedOperation(model: Kratos.Model, parameters: Kratos.Parameters):
    default_settings = Kratos.Parameters("""{
        "model_part_name"    : "PLEASE_SPECIFY_MODEL_PART_NAME",
        "io_settings"        : {},
        "controller_settings": {
            "controller_type": "default_controller"
        },
        "list_of_operations": [
            {
                "operation_type": "PLEASE_SPECIFY_OPERATION_TYPE"
            }
        ]
    }""")
    parameters.AddMissingParameters(default_settings)

    model_part = model[parameters["model_part_name"].GetString()]

    controller = ControllerFactory(model, parameters["controller_settings"])

    aggregated_operation = AggregatedControlledOperations(model_part, parameters["io_settings"])

    operation_params: KratosMultiphysics.Parameters
    for operation_params in parameters["list_of_operations"].values():
        if not operation_params.Has("operation_type"):
            raise RuntimeError(f"\"operation_type\" not found in operation parameters. Parameters = {operation_params}")

        operation_type = operation_params["operation_type"].GetString()
        if not operation_type in operation_types_map.keys():
            raise RuntimeError(f"Unsupported operation_type = \"{operation_type}\". Followings are supported operation types:\n\t" + "\n\t".join(operation_types_map.keys()))

        operation_params.RemoveValue("operation_type")

        aggregated_operation.AddControlledOperation(ControlledOperation(operation_types_map[operation_type], operation_params, controller))

    return aggregated_operation