"""Construct a user-defined HDF5 IO process.

This process:
 - passes the json parameters directly to the core factory.

license: HDF5Application/license.txt
"""


__all__ = ["Factory"]

import typing

import KratosMultiphysics
from KratosMultiphysics.HDF5Application.core.processes import HDF5Process
from KratosMultiphysics.HDF5Application.core.processes import HDF5OutputProcess
from KratosMultiphysics.HDF5Application.core.operations.aggregated_operations import AggregatedControlledOperations
from KratosMultiphysics.HDF5Application.core.operations import CreateAggregatedOperation

def Factory(settings: KratosMultiphysics.Parameters, model: KratosMultiphysics.Model) -> 'typing.Union[HDF5Process, HDF5OutputProcess]':
    """Return a user-defined input/output process for HDF5.

    The input settings are a json array of parameters which maps to the
    structure of the HDF5 IO python core.

    For example:
        '''
        [{
            "model_part_name" : "MainModelPart",
            "process_step": "finalize_solution_step",
            "controller_settings": {
                "controller_type": "temporal_controller",
                "time_frequency": 0.5
            },
            "io_settings": {
                "file_name": "results/<model_part_name>-<time>.h5"
            },
            "list_of_operations": [{
                "operation_type": "model_part_output"
            },{
                "operation_type": "nodal_solution_step_data_output",
                "list_of_variables": ["DISPLACEMENT"]
            }]
        }]
        '''

    will store the model part and displacement every 0.5s in the following
    directory tree structure:

        ./results/MainModelPart-0.000.h5
        ./results/MainModelPart-0.500.h5
        ./results/MainModelPart-1.000.h5
        ...

    and internal file tree structure in each .h5 file:
        /ModelData/Conditions
        /ModelData/Elements
        ...
        /ResultsData/NodalSolutionStepData/DISPLACEMENT

    In the above example, the nonterminal symbols <model_part_name> and <time>
    are automatically replaced by the name of the model part and the current
    time.

    Alternatively, the simulations results can be stored in a single .h5 file
    containing the directory tree structure of the above example:

        '''
        [{
            "model_part_name" : "MainModelPart",
            "process_step": "finalize_solution_step",
            "controller_settings": {
                "controller_type": "temporal_controller",
                "time_frequency": 0.5
            },
            "io_settings": {
                "file_name": "results.h5",
                "file_access_mode": "read_write"
            },
            "list_of_operations": [{
                "prefix": "/<time>/<model_part_name>/ModelData",
                "operation_type": "model_part_output"
            },{
                "prefix": "/<time>/<model_part_name>/ResultsData",
                "operation_type": "nodal_solution_step_data_output",
                "list_of_variables": ["DISPLACEMENT"]
            }]
        }]
        '''

    will store the model part and displacement every 0.5s in the following
    directory tree structure:

        ./results.h5

    and internal file tree structure in each .h5 file:
        /0.000/MainModelPart/ModelData/...
        /0.000/MainModelPart/ResultsData/...
        ...
        /0.500/MainModelPart/ModelData/...
        /0.500/MainModelPart/ResultsData/...
        ...
        /1.000/MainModelPart/ModelData/...
        /1.000/MainModelPart/ResultsData/...

    Different groupings of model parts, files, locations within the solution
    algorithm, frequencies and IO operations can be configured by appending
    additional parameters to the json array.
    """

    defaults = KratosMultiphysics.Parameters("""{
        "model_part_name" : "MainModelPart",
        "process_step": "finalize_solution_step",
        "controller_settings": {
            "controller_type": "temporal_controller",
            "time_frequency": 0.5
        },
        "io_settings": {
            "file_name": "results.h5",
            "file_access_mode": "read_write"
        },
        "list_of_operations": [{
            "prefix": "/<time>/<model_part_name>/ModelData",
            "operation_type": "model_part_output"
        },{
            "prefix": "/<time>/<model_part_name>/ResultsData/",
            "operation_type": "nodal_solution_step_data_output",
            "list_of_variables": ["DISPLACEMENT"]
        }]
    }""")

    aggregated_operations_map: 'dict[str, list[AggregatedControlledOperations]]' = {}

    aggregated_operation_params: KratosMultiphysics.Parameters
    for aggregated_operation_params in settings.values():
        aggregated_operation_params.ValidateAndAssignDefaults(defaults)

        process_step = aggregated_operation_params["process_step"].GetString()
        if process_step not in aggregated_operations_map.keys():
            aggregated_operations_map[process_step]: 'list[AggregatedControlledOperations]' = []
        aggregated_operations_map[process_step].append(CreateAggregatedOperation(model, aggregated_operation_params))

    if "print_output" in aggregated_operations_map.keys():
        process: 'typing.Union[HDF5Process, HDF5OutputProcess]' = HDF5OutputProcess()
    else:
        process: 'typing.Union[HDF5Process, HDF5OutputProcess]' = HDF5Process()

    for process_step, aggregated_operations_list in aggregated_operations_map.items():
        if process_step == "print_output":
            list(map(lambda x: process.AddPrintOutput(x) , aggregated_operations_list))
        else:
            modified_process_step = KratosMultiphysics.StringUtilities.ConvertSnakeCaseToCamelCase(process_step)
            method_name = f"Add{modified_process_step}"
            list(map(lambda x: getattr(process, method_name)(x), aggregated_operations_list))

    return process