import KratosMultiphysics as Kratos

from KratosMultiphysics.RANSApplication.rans_analysis import RANSAnalysis
from KratosMultiphysics.RANSApplication.adjoint_rans_analysis import AdjointRANSAnalysis

from KratosMultiphysics.RANSApplication.formulations.utilities import AddFileLoggerOutput
from KratosMultiphysics.RANSApplication.formulations.utilities import RemoveFileLoggerOutput

from pathlib import Path
from shutil import copy
import math

def SolvePrimalProblem(kratos_primal_parameters_file_name, log_file_name = "primal_evaluation.log"):
    primal_filled_settings_file = Path(kratos_primal_parameters_file_name[:kratos_primal_parameters_file_name.rfind(".")] + "_final.json")

    if (not primal_filled_settings_file.is_file()):
        with open(kratos_primal_parameters_file_name, "r") as file_input:
            primal_parameters = Kratos.Parameters(file_input.read())

        # set the loggers
        default_severity, file_logger = AddFileLoggerOutput(log_file_name)

        # run the primal analysis
        primal_model = Kratos.Model()
        primal_simulation = RANSAnalysis(primal_model, primal_parameters)
        primal_simulation.Run()

        with open(str(primal_filled_settings_file), "w") as file_output:
            file_output.write(primal_parameters.PrettyPrintJsonString())

        del primal_simulation
        del primal_model

        # flush the primal output
        RemoveFileLoggerOutput(default_severity, file_logger)
        Kratos.Logger.PrintInfo("SolvePrimalProblem", "Solved primal evaluation at {}.".format(kratos_primal_parameters_file_name))
    else:
        # following is done to initialize default settings
        with open(str(primal_filled_settings_file), "r") as file_input:
            primal_parameters = Kratos.Parameters(file_input.read())
        Kratos.Logger.PrintInfo("SolvePrimalProblem", "Found existing completed primal evaluation at {}.".format(kratos_primal_parameters_file_name))

    return primal_parameters

def SolveAdjointProblem(model, adjoint_parameters_file_name, log_file_name, skip_if_already_run = False):
    if (skip_if_already_run):
        if Path(log_file_name).is_file():
            with open(log_file_name, "r") as log_file_input:
                last_line = log_file_input.readlines()[-1][:-1].strip()

            if (last_line == "AdjointRANSAnalysis: Analysis -END-"):
                Kratos.Logger.PrintInfo("SolveAdjointProblem", "Found existing completed adjoint evaluation at {}.".format(adjoint_parameters_file_name))
                return

    with open(adjoint_parameters_file_name, "r") as file_input:
        adjoint_parameters = Kratos.Parameters(file_input.read())

    # set the lift loggers
    default_severity, file_logger = AddFileLoggerOutput(log_file_name)

    adjoint_simulation = AdjointRANSAnalysis(model, adjoint_parameters)
    adjoint_simulation.Run()

    RemoveFileLoggerOutput(default_severity, file_logger)

    return adjoint_simulation

def RecursiveCopy(src, dest):
    src_path = Path(src)
    dest_path = Path(dest)
    for item in src_path.iterdir():
        if (item.is_file()):
            copy(str(item), dest)
        elif (item.is_dir()):
            new_path = dest_path / item.relative_to(src)
            new_path.mkdir(exist_ok=True)
            RecursiveCopy(str(item), str(new_path))

def GetDragValues(kratos_parameters, model_part_name):
    output_process = GetDragResponseFunctionOutputProcess(kratos_parameters, model_part_name)
    if (output_process is not None):
        output_file_name = output_process["Parameters"]["output_file_settings"]["file_name"].GetString()
        time_steps, reactions = ReadDrag(output_file_name)
        return time_steps, reactions
    else:
        raise RuntimeError("No \"compute_body_fitted_drag_process\" found in auxiliar_process_list.")

def CalculateTimeAveragedDrag(kratos_parameters, model_part_name, direction, start_time = 0.0):
    time_steps, reactions = GetDragValues(kratos_parameters, model_part_name)
    total_drag = 0.0
    for index, reaction in enumerate(reversed(reactions)):
        if (time_steps[len(time_steps) - index - 1] >= start_time):
            total_drag += reaction[0] * direction[0] + reaction[1] * direction[1] + reaction[2] * direction[2]
            if (kratos_parameters["solver_settings"]["time_scheme_settings"]["scheme_type"].GetString() == "steady"):
                break
    if len(time_steps) > 1:
        delta_time = time_steps[1] - time_steps[0]
        total_drag *= delta_time
    return total_drag

def ReadDrag(file_name):
    with open(file_name, "r") as file_input:
        lines = file_input.readlines()
    time_steps = []
    reaction = []
    for line in lines:
        line = line.strip()
        if len(line) == 0 or line[0] == '#':
            continue
        time_step_data = [float(v) for v in line.split()]
        time, fx, fy, fz = time_step_data
        time_steps.append(time)
        reaction.append([fx, fy, fz])
    return time_steps, reaction

def GetDragResponseFunctionOutputProcess(kratos_parameters, model_part_name):
    auxiliar_process_list = kratos_parameters["processes"]["auxiliar_process_list"]
    for process_settings in auxiliar_process_list:
        if (
            process_settings.Has("python_module") and process_settings["python_module"].GetString() == "compute_body_fitted_drag_process" and
            process_settings.Has("kratos_module") and process_settings["kratos_module"].GetString() == "KratosMultiphysics.FluidDynamicsApplication" and
            process_settings["Parameters"].Has("model_part_name") and process_settings["Parameters"]["model_part_name"].GetString() == model_part_name
            ):
            return process_settings

    return None

def GetResponseFunctionOutputProcess(kratos_parameters, model_part_name, response_function_parameters):
    output_processes_categories_list = kratos_parameters["output_processes"]
    for _, value in output_processes_categories_list.items():
        for process_settings in value:
            if (
                process_settings.Has("python_module") and process_settings["python_module"].GetString() == "response_function_output_process" and
                process_settings.Has("kratos_module") and process_settings["kratos_module"].GetString() == "KratosMultiphysics.FluidDynamicsApplication"):
                # found a reponse function output process

                process_parameters = process_settings["Parameters"]
                is_valid_respones_function = True

                is_valid_respones_function = is_valid_respones_function and process_parameters["response_type"].GetString() ==  response_function_parameters["response_type"].GetString()
                is_valid_respones_function = is_valid_respones_function and model_part_name.startswith(process_parameters["model_part_name"].GetString())
                is_valid_respones_function = is_valid_respones_function and (response_function_parameters["custom_settings"].IsKeysSubSetWithEquivalentValuesTo(process_parameters["response_settings"]))

                if is_valid_respones_function:
                    return process_settings

    return None

def ReadResponseValuesFile(file_name):
    with open(file_name, "r") as file_input:
        lines = file_input.readlines()
    time_steps = []
    values = []
    for line in lines:
        line = line.strip()
        if len(line) == 0 or line[0] == '#':
            continue
        time_step_data = [float(v) for v in line.split(sep=",")]
        time, value = time_step_data
        time_steps.append(time)
        values.append(value)
    return time_steps, values

def GetResponseValues(kratos_parameters, model_part_name, response_function_parameters):
    output_process = GetResponseFunctionOutputProcess(kratos_parameters, model_part_name, response_function_parameters)
    if (output_process is not None):
        file_path = Path(".")
        if output_process["Parameters"]["output_file_settings"].Has("output_path"):
            file_path = file_path / output_process["Parameters"]["output_file_settings"]["output_path"].GetString()
        output_file_name = str(file_path / str(output_process["Parameters"]["output_file_settings"]["file_name"].GetString() + ".dat"))
        time_steps, values = ReadResponseValuesFile(output_file_name)
        return time_steps, values
    else:
        raise RuntimeError("No \"response_function_output_process\" found in auxiliar_process_list matching {:s} model part".format(model_part_name))

def CalculateTimeAveragedResponseValue(kratos_parameters, model_part_name, response_function_parameters, start_time = 0.0):
    time_steps, values = GetResponseValues(kratos_parameters, model_part_name, response_function_parameters)
    total_value = 0.0
    for index, value in enumerate(reversed(values)):
        if (time_steps[len(time_steps) - index - 1] >= start_time):
            total_value += value
            if (kratos_parameters["solver_settings"]["time_scheme_settings"]["scheme_type"].GetString() == "steady"):
                break
    if len(time_steps) > 1:
        delta_time = time_steps[1] - time_steps[0]
        total_value *= delta_time
    return total_value

def WriteShapeSensitivities(model_part, file_name) -> dict:
    # now write shape sensitivities to results file
    shape_sensitivities_dict = {}
    with open(file_name, "w") as file_output:
        file_output.write("# node_id, shape_sensitivity_x, shape_sensitivity_y, shape_sensitivity_z\n")
        for node in model_part.Nodes:
            shape_sensitivity = node.GetSolutionStepValue(Kratos.SHAPE_SENSITIVITY)
            shape_sensitivities_dict[node.Id] = shape_sensitivity
            file_output.write("{:d},{:0.16e},{:0.16e},{:0.16e}\n".format(node.Id, shape_sensitivity[0], shape_sensitivity[1], shape_sensitivity[2]))

    return shape_sensitivities_dict

def CalculateShapeSensitivity(model, optimization_model_part_name, adjoint_parameters_file_name, log_file_name) -> dict:
    base_name = adjoint_parameters_file_name[:adjoint_parameters_file_name.rfind(".")]
    shape_sensitivity_file = Path(base_name + "_results.csv")

    if not shape_sensitivity_file.is_file():
        _ = SolveAdjointProblem(model, adjoint_parameters_file_name, log_file_name)
        return WriteShapeSensitivities(model[optimization_model_part_name], str(shape_sensitivity_file))
    else:
        Kratos.Logger.PrintInfo("SolveAdjointProblem", "Found existing completed adjoint evaluation at {}.".format(str(shape_sensitivity_file)))

        with open(str(shape_sensitivity_file), "r") as file_input:
            lines = file_input.readlines()

        lines = lines[1:-1]
        shape_sensitivities_dict = {}
        for line in lines:
            data_str = line[:-1].split(",")
            shape_sensitivities_dict[int(data_str[0])] = Kratos.Array3([float(data_str[1]), float(data_str[2]), float(data_str[3])])

        return shape_sensitivities_dict

def UpdateStringWithPlaceHolders(input: str, place_holder_dict: dict) -> str:
    for k, v in place_holder_dict.items():
        input = input.replace(k, v)
    return input

def UpdateFilesWithPlaceHolders(settings: Kratos.Parameters, place_holder_dict: dict):
    default_parameters = Kratos.Parameters("""{
        "original_file_name": "",
        "updated_file_name" : ""
    }""")

    settings.ValidateAndAssignDefaults(default_parameters)
    original_file_name = settings["original_file_name"].GetString()
    if original_file_name == "":
        raise Exception("Original file name cannot be empty.")

    updated_file_name = settings["updated_file_name"].GetString()
    if updated_file_name == "":
        raise Exception("Updated file name cannot be empty.")

    updated_file_name = UpdateStringWithPlaceHolders(updated_file_name, place_holder_dict)

    if original_file_name == updated_file_name:
        raise Exception("Original file name and updated file name cannot be the same. [ Original file name = {:s}, updated file name = {:s}".format(original_file_name, updated_file_name))

    with open(original_file_name, "r") as file_input:
        data = file_input.read()

    data = UpdateStringWithPlaceHolders(data, place_holder_dict)

    if Path(updated_file_name).is_file():
        Kratos.Logger.PrintWarning("UpdateFilesWithPlaceHolders", "Found an existing file at {:s}. This will overwrite it.".format(str(Path(updated_file_name).absolute())))

    with open(updated_file_name, "w") as file_output:
        file_output.write(data)

    Kratos.Logger.PrintInfo("UpdateFilesWithPlaceHolders", "{:s} is updated with place holder values at {:s}.".format(original_file_name, updated_file_name))


