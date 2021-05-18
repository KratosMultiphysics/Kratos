import KratosMultiphysics as Kratos
from KratosMultiphysics.kratos_utilities import DeleteFileIfExisting
from KratosMultiphysics.FluidDynamicsApplication.fluid_dynamics_analysis import FluidDynamicsAnalysis
from KratosMultiphysics.FluidDynamicsApplication.adjoint_fluid_analysis import AdjointFluidAnalysis

import os


def DeleteH5Files():
    for name in os.listdir():
        if name.endswith(".h5"):
            DeleteFileIfExisting(name)


def ReadParameters(parameters_file_name):
    with open(parameters_file_name, 'r') as file_input:
        return Kratos.Parameters(file_input.read())


def SolvePrimalProblem(kratos_parameters):
    test = FluidDynamicsAnalysis(Kratos.Model(), kratos_parameters)
    test.Run()

    return test


def SolveAdjointProblem(kratos_parameters):
    test = AdjointFluidAnalysis(Kratos.Model(), kratos_parameters)
    test.Run()

    return test


def ComputeAdjointSensitivity(node_ids, kratos_parameters, adjoint_problem_solving_method):
    test = adjoint_problem_solving_method(kratos_parameters)

    sensitivities = Kratos.Matrix(len(node_ids), 2)
    for i, node_id in enumerate(node_ids):
        node_sensitivity = test._GetSolver().main_model_part.GetNode(
            node_id).GetSolutionStepValue(Kratos.SHAPE_SENSITIVITY)
        sensitivities[i, 0] = node_sensitivity[0]
        sensitivities[i, 1] = node_sensitivity[1]

    return sensitivities


def ComputeFiniteDifferenceSensitivity(node_ids, step_size,
                                       kratos_parameters,
                                       objective_value_evaluation_method,
                                       primal_output_process_inclusion_method):
    model_import_settings = kratos_parameters["solver_settings"]["model_import_settings"]
    input_type = model_import_settings["input_type"].GetString()
    if (input_type != "mdpa"):
        raise Exception(
            "Unsupported \"input_type\". Only mdpa is supported. [ \"input_type\" = \"{:s}\" ].".format(input_type))

    model_part_file_name = model_import_settings["input_filename"].GetString(
    )
    with open(model_part_file_name + '.mdpa', 'r') as file_input:
        lines = file_input.readlines()

    node_block_start = lines.index('Begin Nodes\n')
    node_block_end = lines.index('End Nodes\n')
    node_lines = lines[node_block_start:node_block_end]

    sensitivity = Kratos.Matrix(len(node_ids), 2)
    for i, node_id in enumerate(node_ids):
        coord = [
            float(node_lines[node_id].split()[i]) for i in range(1, 4)
        ]

        for k in range(2):
            coord[k] += step_size
            _writeNodalCoordinates(
                list(lines), node_block_start, node_lines, node_id, coord,
                model_part_file_name)
            sensitivity[i, k] = objective_value_evaluation_method(
                kratos_parameters.Clone())
            coord[k] -= step_size

    # now compute the unperturbed
    with open(model_part_file_name + '.mdpa', 'w') as model_part_file:
        model_part_file.writelines(lines)

    unperturbed_settings = kratos_parameters.Clone()
    primal_output_process_inclusion_method(unperturbed_settings)
    ref_value = objective_value_evaluation_method(unperturbed_settings)

    # update finite difference sensitivities to correct values
    for i in range(sensitivity.Size1()):
        for j in range(sensitivity.Size2()):
            sensitivity[i, j] = (sensitivity[i, j] - ref_value) / step_size

    return sensitivity


def _writeNodalCoordinates(lines, node_block_start, node_lines, node_id,
                           coords, model_part_file_name):
    old_line = node_lines[node_id]
    components = old_line.split()
    if int(components[0]) != node_id:
        raise RuntimeError('Error parsing file ' + model_part_file_name)
    new_line = '{:5d}'.format(node_id) + ' ' \
        + '{:19.10f}'.format(coords[0]) + ' ' \
        + '{:19.10f}'.format(coords[1]) + ' ' \
        + '{:19.10f}'.format(coords[2]) + '\n'
    lines[node_block_start + node_id] = new_line
    with open(model_part_file_name + '.mdpa', 'w') as model_part_file:
        model_part_file.writelines(lines)


class FiniteDifferenceBodyFittedDragShapeSensitivityAnalysis:
    @staticmethod
    def ComputeSensitivity(node_ids, step_size,
                           kratos_parameters,
                           drag_direction,
                           drag_model_part_name,
                           primal_problem_solving_method,
                           primal_output_process_inclusion_method,
                           is_steady = False):

        FiniteDifferenceBodyFittedDragShapeSensitivityAnalysis._AddResponseFunctionOutput(
            kratos_parameters, drag_model_part_name)

        output_file_name = drag_model_part_name + "_drag.dat"

        def compute_drag(kratos_parameters):
            primal_problem_solving_method(kratos_parameters)
            return FiniteDifferenceBodyFittedDragShapeSensitivityAnalysis._GetTimeAveragedDrag(drag_direction, output_file_name, is_steady)

        sensitivities = ComputeFiniteDifferenceSensitivity(node_ids, step_size,
                                                           kratos_parameters, compute_drag, primal_output_process_inclusion_method)

        DeleteFileIfExisting(output_file_name)
        return sensitivities

    @staticmethod
    def _GetTimeAveragedDrag(direction, drag_file_name, is_steady):
        time_steps, reactions = FiniteDifferenceBodyFittedDragShapeSensitivityAnalysis._ReadDrag(
            drag_file_name)
        total_drag = 0.0
        for reaction in reversed(reactions):
            total_drag += reaction[0]*direction[0] + \
                reaction[1]*direction[1]+reaction[2]*direction[2]
            if (is_steady):
                break
        if len(time_steps) > 1:
            delta_time = time_steps[1] - time_steps[0]
            total_drag *= delta_time
        return total_drag

    @staticmethod
    def _ReadDrag(file_name):
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

    @staticmethod
    def _AddResponseFunctionOutput(kratos_parameters, drag_model_part_name):
        response_parameters = Kratos.Parameters(R'''
        {
            "python_module": "compute_body_fitted_drag_process",
            "kratos_module": "KratosMultiphysics.FluidDynamicsApplication",
            "process_name": "ComputeBodyFittedDragProcess",
            "Parameters": {
                "model_part_name": "<DRAG_MODEL_PART_NAME>",
                "interval": [
                    0.0,
                    1e30
                ],
                "write_drag_output_file": true,
                "print_drag_to_screen": false,
                "print_format": "22.15e"
            }
        }
        ''')

        response_parameters["Parameters"]["model_part_name"].SetString(
            drag_model_part_name)

        kratos_parameters["processes"]["auxiliar_process_list"].Append(
            response_parameters)


class FiniteDifferenceVelocityPressureNormSquareShapeSensitivityAnalysis:
    @staticmethod
    def ComputeSensitivity(node_ids, step_size,
                           kratos_parameters,
                           norm_model_part_name,
                           primal_problem_solving_method,
                           primal_output_process_inclusion_method,
                           is_steady = False):

        FiniteDifferenceVelocityPressureNormSquareShapeSensitivityAnalysis._AddResponseFunctionOutput(
            kratos_parameters, norm_model_part_name)

        def compute_drag(kratos_parameters):
            primal_problem_solving_method(kratos_parameters)
            return FiniteDifferenceVelocityPressureNormSquareShapeSensitivityAnalysis._ObjectiveValueEvaluation(is_steady)

        sensitivities = ComputeFiniteDifferenceSensitivity(node_ids, step_size,
                                                           kratos_parameters, compute_drag, primal_output_process_inclusion_method)

        DeleteFileIfExisting(
            "velocity_pressure_norm_square_response_output.dat")
        return sensitivities

    @staticmethod
    def _ObjectiveValueEvaluation(is_steady):
        with open("velocity_pressure_norm_square_response_output.dat", "r") as file_input:
            lines = file_input.readlines()

        value = 0.0
        total_time = 0.0
        delta_time = 0.0
        for line in reversed(lines[2:]):
            data = line.strip().split(",")
            delta_time = float(data[0]) - total_time
            total_time = float(data[0])
            value += float(data[1])
            if (is_steady):
                break

        return value * delta_time

    @staticmethod
    def _AddResponseFunctionOutput(kratos_parameters, norm_model_part_name):
        response_parameters = Kratos.Parameters(R'''
        {
            "python_module": "response_function_output_process",
            "kratos_module": "KratosMultiphysics.FluidDynamicsApplication",
            "process_name": "ResponseFunctionOutputProcess",
            "Parameters": {
                "response_type": "norm_square",
                "model_part_name": "<PARENT_MODEL_PART>",
                "response_settings": {
                    "main_model_part_name": "<PARENT_MODEL_PART>",
                    "norm_model_part_name": "<SUBMODEL_PART>",
                    "entities": [
                        "conditions"
                    ]
                },
                "output_file_settings": {
                    "file_name": "velocity_pressure_norm_square_response_output"
                }
            }
        }
        ''')
        seperator_index = norm_model_part_name.rfind('.')
        response_parameters["Parameters"]["model_part_name"].SetString(
            norm_model_part_name[:seperator_index])
        response_parameters["Parameters"]["response_settings"]["main_model_part_name"].SetString(
            norm_model_part_name[:seperator_index])
        response_parameters["Parameters"]["response_settings"]["norm_model_part_name"].SetString(
            norm_model_part_name)

        kratos_parameters["output_processes"].AddEmptyList("response_function_outputs")
        kratos_parameters["output_processes"]["response_function_outputs"].Append(response_parameters)
