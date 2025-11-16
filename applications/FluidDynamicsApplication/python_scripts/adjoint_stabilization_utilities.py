import numpy as np
import math
import os
import sys

from pathlib import Path

import KratosMultiphysics as Kratos
import KratosMultiphysics.StatisticsApplication as KratosStats

def AddFileLoggerOutput(log_file_name):
    file_logger = Kratos.FileLoggerOutput(log_file_name)
    default_severity = Kratos.Logger.GetDefaultOutput().GetSeverity()
    Kratos.Logger.GetDefaultOutput().SetSeverity(Kratos.Logger.Severity.WARNING)
    Kratos.Logger.AddOutput(file_logger)

    return default_severity, file_logger

def RemoveFileLoggerOutput(default_severity, file_logger):
    Kratos.Logger.Flush()
    Kratos.Logger.RemoveOutput(file_logger)
    Kratos.Logger.GetDefaultOutput().SetSeverity(default_severity)

class ExecutionScope:
    def __init__(self, execution_path):
        self.currentPath = Path.cwd()
        self.scope = Path(execution_path)

    def __enter__(self):
        self.scope.mkdir(parents=True, exist_ok=True)
        sys.path.append(str(self.scope.absolute()))
        os.chdir(str(self.scope))

    def __exit__(self, exc_type, exc_value, traceback):
        os.chdir(self.currentPath)
        sys.path.remove(str(self.scope.absolute()))

def SolveProblem(analysis_class_type, kratos_parameters, execution_prefix):
    # set the loggers
    default_severity, file_logger = AddFileLoggerOutput(execution_prefix + ".log")

    # run the primal analysis
    model = Kratos.Model()
    primal_simulation = analysis_class_type(model, kratos_parameters)
    primal_simulation.Run()

    with open(execution_prefix + ".json", "w") as file_output:
        file_output.write(kratos_parameters.PrettyPrintJsonString())

    # flush the primal output
    RemoveFileLoggerOutput(default_severity, file_logger)
    Kratos.Logger.PrintInfo("SolvePrimalProblem", "Solved primal evaluation at {}.".format(execution_prefix + ".json"))

    return model, primal_simulation

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

    for i, data_item in enumerate(data[1:]):
        if data_item[0] == data[i-1][0]:
            data[i][1] = data[i-1][1]

    data = sorted(data, key=lambda x: x[0])
    return data

def __IsPlateau(slope, max_slope):
    return slope <= max_slope

def _CalculateTimeSeriesSlope(time_series_values, time_range):
    if np.isfinite(time_series_values[-1, 1]):
        time_values = time_series_values[:, 0]
        max_time = np.max(time_values)
        indices_range = np.where(np.logical_and(time_values >= time_range[0], time_values <= time_range[1]))[0]
        ranged_time_values = max_time - np.take(time_values, indices_range)
        ranged_values = np.take(time_series_values[:, 1], indices_range)[:-1] / ranged_time_values[1:]
        if len(ranged_time_values[1:]) == 0:
            raise Exception("Stabilization analysis plateau time range is not within the adjoint evaluation time range.")
        return np.polyfit(ranged_time_values[1:], ranged_values, 1)[0]
    else:
        return 1e+300


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
    while (coefficient_data_list[2][0] - coefficient_data_list[0][0] > tolerance and iteration < max_iterations):
        iteration += 1
        msg = "Summary of bisection results in iteration {:d}/{:d}:".format(iteration, max_iterations)

        for index, coefficient_data in enumerate(coefficient_data_list):
            if coefficient_data[1] is None:
                if (index > 0 and coefficient_data_list[index-1][1] == 0) or (index == 0):
                    solve_id += 1
                    coefficient_data[1] = __IsPlateau(
                                            _CalculateTimeSeriesSlope(
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
                    Kratos.Logger.PrintInfo("Stabilization Analysis", "Stabilization coefficient {:f} resulted in a plateau : {:d}".format(coefficient_data[0], coefficient_data[1]))
                else:
                    Kratos.Logger.PrintInfo("Stabilization Analysis", "Found plateau in previous calculations, therefore assuming stabilization coefficient {:f} as plateau.".format(coefficient_data[0]))
                    coefficient_data[1] = 1

            msg += "\n       point {:d}: stabilization coefficient = {:e}, is plateau = {:d}".format(index + 1, coefficient_data[0], int(coefficient_data[1]))

        Kratos.Logger.PrintInfo("Stabilization Analysis", msg)

        coefficient_data_list = __Bisect(coefficient_data_list)

    if iteration == max_iterations:
        Kratos.Logger.PrintWarning("Stabilization Analysis", "Adjoint stabilization computation reached maximum iterations.")

    stabilization_coefficient = (coefficient_data_list[2][0] + coefficient_data_list[0][0]) / 2.0
    Kratos.Logger.PrintInfo("Stabilization Analysis", "Computed adjoint stabilization coefficient = {:e}".format(stabilization_coefficient))
    return stabilization_coefficient

