"""Construct a user-defined HDF5 IO process.

This process:
 - passes the json parameters directly to the core factory.

license: HDF5Application/license.txt
"""


__all__ = ["Factory"]


import KratosMultiphysics
from KratosMultiphysics.HDF5Application import core
from KratosMultiphysics.HDF5Application.utils import ParametersWrapper


def Factory(settings, Model):
    """Return a user-defined input/output process for HDF5.

    The input settings are a json array of parameters which maps to the
    structure of the HDF5 IO python core.

    The settings of each array item are given in the following table:
    +-----------------------+------------+-------------------------------------------+
    | Setting               | Type       | Default Value                             |
    +-----------------------+------------+-------------------------------------------+
    | "model_part_name"     | String     | "PLEASE_SPECIFY_MODEL_PART_NAME"          |
    +-----------------------+------------+-------------------------------------------+
    | "process_step"        | String     | "initialize"                              |
    +-----------------------+------------+-------------------------------------------+
    | "controller_settings" | Parameters | {                                         |
    |                       |            |   "controller_type": "default_controller" |
    |                       |            | }                                         |
    +-----------------------+------------+-------------------------------------------+
    | "io_settings"         | Parameters | "echo_level": 0                           |
    |                       |            | "file_access_mode": "exclusive"           |
    |                       |            | "file_driver": "sec2"                     |
    |                       |            | "file_name": "kratos"                     |
    |                       |            | "max_files_to_keep": "unlimited"          |
    |                       |            | "io_type": "serial_hdf5_file_io"          |
    +-----------------------+------------+-------------------------------------------+
    | "list_of_operations"  | Parameters | [{                                        |
    |                       | Array      |   "operation_type": "model_part_output"   |
    |                       |            |   "prefix": "/ModelData"                  |
    |                       |            | }]                                        |
    +-----------------------+------------+-------------------------------------------+


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
    # TODO: decide whether to pass a KratosMultiphysics.Process or a KratosMultiphysics.OutputProcess
    return core.Factory(ParametersWrapper(settings["Parameters"]), Model, KratosMultiphysics.Process)
