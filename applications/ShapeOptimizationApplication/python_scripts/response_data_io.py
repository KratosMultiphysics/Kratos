import json
from pathlib import Path

import KratosMultiphysics as Kratos
from KratosMultiphysics.response_functions.response_function_interface import ResponseFunctionInterface

def GetObjectiveValue(response: ResponseFunctionInterface, identifier: str, optimization_iteration: int, optimization_results_absolute_path: str, is_iteration_restart_files_written: bool) -> float:
    if is_iteration_restart_files_written:
        iteration_data_file_path = Path(optimization_results_absolute_path) / Path("Raw_Data/iteration_{:05d}.json".format(optimization_iteration))

        if not iteration_data_file_path.is_file():
            iteration_data_file_path.parent.mkdir(parents=True, exist_ok=True)
            # write dummy dict to create the file
            with open(str(iteration_data_file_path), "w") as file_output:
                file_output.write(json.dumps({}, indent=4, sort_keys=True))

        with open(str(iteration_data_file_path), "r") as file_input:
            data_dict = json.loads(file_input.read())

        if identifier not in data_dict.keys():
            data_dict[identifier] = {}

        if "value" not in data_dict[identifier].keys():
            response.CalculateValue()
            response_value = response.GetValue()
            data_dict[identifier]["value"] = response_value

            # now write the objective value
            with open(str(iteration_data_file_path), "w") as file_output:
                file_output.write(json.dumps(data_dict, indent=4, sort_keys=True))
        else:
            response_value = data_dict[identifier]["value"]
            Kratos.Logger.PrintInfo("ShapeOpt", "Found recorded objective value for iteration {:d} {:s} objective. Remove \"{:s}/value\" entry from {:s} to re-calculate objective value for this iteration.".format(optimization_iteration, identifier, identifier, str(iteration_data_file_path.absolute())))
    else:
        response.CalculateValue()
        response_value = response.GetValue()

    return response_value

def GetObjectiveSensitivityValues(response: ResponseFunctionInterface, identifier: str, optimization_iteration: int, optimization_results_absolute_path: str, is_iteration_restart_files_written: bool) -> dict:
    if is_iteration_restart_files_written:
        iteration_data_file_path = Path(optimization_results_absolute_path) / Path("Raw_Data/iteration_{:05d}.json".format(optimization_iteration))

        if not iteration_data_file_path.is_file():
            iteration_data_file_path.parent.mkdir(parents=True, exist_ok=True)
            # write dummy dict to create the file
            with open(str(iteration_data_file_path), "w") as file_output:
                file_output.write(json.dumps({}, indent=4, sort_keys=True))

        with open(str(iteration_data_file_path), "r") as file_input:
            data_dict = json.loads(file_input.read())

        if identifier not in data_dict.keys():
            data_dict[identifier] = {}

        if "sensitivities" not in data_dict[identifier].keys():
            response.CalculateGradient()
            sensitivities = response.GetNodalGradient(Kratos.SHAPE_SENSITIVITY)
            converted_sensitivities = {}
            for node_id, nodal_sensitivity in sensitivities.items():
                converted_sensitivities[node_id] = [nodal_sensitivity[0], nodal_sensitivity[1], nodal_sensitivity[2]]
            data_dict[identifier]["sensitivities"] = converted_sensitivities

            # now write the sensitivities
            with open(str(iteration_data_file_path), "w") as file_output:
                file_output.write(json.dumps(data_dict, indent=4, sort_keys=True))
        else:
            sensitivities = {}
            converted_sensitivities = data_dict[identifier]["sensitivities"]
            for node_id, nodal_sensitivity in converted_sensitivities.items():
                sensitivities[int(node_id)] = Kratos.Array3(nodal_sensitivity)
            Kratos.Logger.PrintInfo("ShapeOpt", "Found recorded sensitivity values for iteration {:d} {:s} objective. Remove \"{:s}/sensitivities\" entry from {:s} to re-calculate objective sensitivities for this iteration.".format(optimization_iteration, identifier, identifier, str(iteration_data_file_path.absolute())))
    else:
        response.CalculateGradient()
        sensitivities = response.GetNodalGradient(Kratos.SHAPE_SENSITIVITY)

    return sensitivities
