import numpy as np
import math

from pathlib import Path

import KratosMultiphysics as Kratos
import KratosMultiphysics.StatisticsApplication as KratosStats

from KratosMultiphysics.RANSApplication.formulations.utilities import ExecutionScope
from KratosMultiphysics.RANSApplication.formulations.utilities import SolveProblem

def __Bisect(data):
    if len(data) != 3:
        raise Exception("Bisection requires 3 points of evaluations")

    data = sorted(data, key=lambda x:x[0])

    lower_bound = -1
    for data_item in data:
        if not data_item[1]:
            lower_bound = data_item[0]
            lower_value = data_item[1]

    upper_bound = -1
    for data_item in reversed(data):
        if data_item[1]:
            upper_bound = data_item[0]
            upper_value = data_item[1]

    if lower_bound == -1:
        lower_bound = data[0][0] / 2.0
        lower_value = None

    if upper_bound == -1:
        upper_bound = data[2][0] * 2.0
        upper_value = None

    data = [
        [lower_bound, lower_value],
        [(lower_bound + upper_bound) * 0.5, None],
        [upper_bound, upper_value]
    ]

    data = sorted(data, key=lambda x: x[0])
    return data

def __IsPlateau(slope, max_slope):
    return math.fabs(slope) <= max_slope

def __CalculateTimeSeriesSlope(time_series_values, time_range):
    time_values = time_series_values[:, 0]
    max_time = np.max(time_values)
    indices_range = np.where(np.logical_and(time_values >= time_range[0], time_values <= time_range[1]))[0]
    ranged_time_values = max_time - np.take(time_values, indices_range)
    ranged_values = np.take(time_series_values[:, 1], indices_range)[:-1] / ranged_time_values[1:]
    return np.polyfit(ranged_time_values[1:], ranged_values, 1)[0]


def __ExecuteAnalysis(analysis_class_type, adjoint_parameters, stabilization_coefficient, solve_id):
    Kratos.Logger.PrintInfo("Stabilization Analysis", "Solving adjoint problem with {:e} stabilization coefficient.".format(stabilization_coefficient))
    scheme_settings = adjoint_parameters["solver_settings"]["scheme_settings"]
    if not scheme_settings.Has("stabilization_coefficient"):
        scheme_settings.AddEmptyValue("stabilization_coefficient")
    scheme_settings["stabilization_coefficient"].SetDouble(stabilization_coefficient)

    with ExecutionScope("adjoint_stabilization_analysis/{:d}".format(solve_id)):
        execute_analysis = not Path("adjoint_stabilization_data.dat").is_file()
        if not execute_analysis:
            if Path("stabilization_coefficient.json").is_file():
                with open("stabilization_coefficient.json", "r") as file_input:
                    kratos_parameters = Kratos.Parameters(file_input.read())
                execute_analysis = kratos_parameters["solver_settings"]["scheme_settings"]["stabilization_coefficient"].GetDouble() != stabilization_coefficient
            else:
                execute_analysis = True

        if execute_analysis:
            _, simulation = SolveProblem(analysis_class_type, adjoint_parameters, "stabilization_coefficient")
            return simulation.time_series_data
        else:
            Kratos.Logger.PrintInfo("Stabilization Analysis", "Found existing solution at {:s}".format(str(Path(".").absolute())))
            return np.loadtxt("adjoint_stabilization_data.dat")

def ComputeStabilizationCoefficient(analysis_class_type, stabilization_settings, execution_method=__ExecuteAnalysis):
    default_parameters = Kratos.Parameters("""{
        "initial_coefficient_bounds"  : [0.0, 1.0],
        "tolerance"                   : 1e-4,
        "plateau_time_range"          : [20.0, 40.0],
        "plateau_max_slope"           : 0.1,
        "max_iterations"              : 20,
        "adjoint_parameters_file_name": "PLEASE_SPECIFY_ADJOINT_KRATOS_PARAMETERS_JSON_FILE_NAME"
    }""")

    stabilization_settings.ValidateAndAssignDefaults(default_parameters)

    class StabilizationAnalysisClass(analysis_class_type):
        def __init__(self, model, parameters):
            super().__init__(model, parameters)
            self.time_series_data = []

        def OutputSolutionStep(self):
            super().OutputSolutionStep()

            rms = KratosStats.SpatialMethods.Historical.NormMethods.RootMeanSquare(
                        self._GetSolver().main_model_part,
                        Kratos.SHAPE_SENSITIVITY,
                        "magnitude")

            self.time_series_data.append([self.time, rms])

        def Finalize(self):
            super().Finalize()

            self.time_series_data = np.array(self.time_series_data)

            # now write time series data for restarting
            np.savetxt("adjoint_stabilization_data.dat", self.time_series_data, header="#time, L2_norm_shape_sensitivity")

    plateau_max_slope = stabilization_settings["plateau_max_slope"].GetDouble()
    plateau_time_range = stabilization_settings["plateau_time_range"].GetVector()
    coefficient_bounds = stabilization_settings["initial_coefficient_bounds"].GetVector()
    tolerance = stabilization_settings["tolerance"].GetDouble()

    with open(stabilization_settings["adjoint_parameters_file_name"].GetString(), "r") as file_input:
        adjoint_parameters = Kratos.Parameters(file_input.read())

    coefficient_data_list = [
        [coefficient_bounds[0], None],
        [(coefficient_bounds[0] + coefficient_bounds[1]) * 0.5, None],
        [coefficient_bounds[1], None]
    ]

    solve_id = 0
    iteration = 0
    max_iterations = stabilization_settings["max_iterations"].GetInt()
    while (coefficient_data_list[2][0] - coefficient_data_list[0][0] > tolerance and iteration <= max_iterations):
        iteration += 1
        msg = "Summary of bisection results in iteration {:d}/{:d}:".format(iteration, max_iterations)

        for index, coefficient_data in enumerate(coefficient_data_list):
            if coefficient_data[1] is None:
                solve_id += 1
                coefficient_data[1] = __IsPlateau(
                                        __CalculateTimeSeriesSlope(
                                            execution_method(
                                                StabilizationAnalysisClass,
                                                adjoint_parameters.Clone(),
                                                coefficient_data[0],
                                                solve_id
                                            ),
                                            plateau_time_range
                                        ),
                                    plateau_max_slope
                                    )

            msg += "\n       point {:d}: stabilization coefficient = {:e}, is plateau = {:d}".format(index + 1, coefficient_data[0], int(coefficient_data[1]))

        Kratos.Logger.PrintInfo("Stabilization Analysis", msg)

        coefficient_data_list = __Bisect(coefficient_data_list)

    stabilization_coefficient = (coefficient_data_list[2][0] + coefficient_data_list[0][0]) / 2.0
    Kratos.Logger.PrintInfo("Stabilization Analysis", "Computed adjoint stabilization coefficient = {:e}".format(stabilization_coefficient))
    return stabilization_coefficient
