#    |  /           |
#    ' /   __| _` | __|  _ \   __|
#    . \  |   (   | |   (   |\__ `
#   _|\_\_|  \__,_|\__|\___/ ____/
#                   Multi-Physics
#
#  License:		 BSD License
#					 license: OptimizationApplication/license.txt
#
#  Main authors:    Suneth Warnakulasuriya
#

import time as timer
import datetime
import KratosMultiphysics as Kratos
from functools import wraps
import typing


def AddFileLoggerOutput(logger_file_name):
    logger_file = Kratos.FileLoggerOutput(logger_file_name)
    default_severity = Kratos.Logger.GetDefaultOutput().GetSeverity()
    Kratos.Logger.GetDefaultOutput().SetSeverity(Kratos.Logger.Severity.WARNING)
    Kratos.Logger.AddOutput(logger_file)

    return default_severity, logger_file

def RemoveFileLoggerOutput(default_severity, logger_file):
    Kratos.Logger.Flush()
    Kratos.Logger.RemoveOutput(logger_file)
    Kratos.Logger.GetDefaultOutput().SetSeverity(default_severity)

class FileLogger:
    """@brief A context responsible for managing the lifetime of logger files."""
    def __init__(self, logger_file_name: str):
        self.__logger_file_name = logger_file_name

    def __enter__(self):
        self.__default_severity, self.__logger_file = AddFileLoggerOutput(self.__logger_file_name)

    def __exit__(self, exit_type, exit_value, exit_traceback):
        RemoveFileLoggerOutput(self.__default_severity, self.__logger_file)

class TimeLogger:
    """@brief A context responsible for outputting execution times."""
    def __init__(self, topic: str, entry_msg: str , exit_msg: str, end_new_line: bool=True):
        self.topic = topic
        self.entry_msg = entry_msg
        self.exit_msg = exit_msg
        self.start_time = None
        self.end_new_line = end_new_line

    def __enter__(self):
        if self.entry_msg is not None:
            Kratos.Logger.PrintInfo(self.topic, self.entry_msg)
        self.start_time = timer.time()

    def __exit__(self, exit_type, exit_value, exit_traceback):
        if self.exit_msg is not None:
            elapsed_time = timer.time() - self.start_time
            if self.end_new_line:
                Kratos.Logger.PrintInfo(self.topic, "{:s} - [ Elapsed time: {:s} ]".format(self.exit_msg, str(datetime.timedelta(seconds=round(elapsed_time))))+"\n")
            else:
                Kratos.Logger.PrintInfo(self.topic, "{:s} - [ Elapsed time: {:s} ]".format(self.exit_msg, str(datetime.timedelta(seconds=round(elapsed_time)))))

class OptimizationAlgorithmTimeLogger:
    """@brief A context responsible for outputting execution times."""
    def __init__(self, optimizer_name: str, optimization_itr: int):
        self.start_time = None
        self.optimization_itr = optimization_itr
        self.optimizer_name = optimizer_name

    def __enter__(self):
        self.start_time = timer.time()

    def __exit__(self, exit_type, exit_value, exit_traceback):

        current_time = timer.time()
        elapsed_time = current_time - self.start_time
        elapsed_time_string = str(datetime.timedelta(seconds=round(elapsed_time)))

        iteration_text = f"{self.optimizer_name} EoF Iteration {self.optimization_itr}"
        iteration_output = f"{'#'}  {iteration_text} [Elapsed Time: {elapsed_time_string}]  {'#'}"

        divided_line = len(iteration_output) * '#'

        to_print = f"{divided_line}\n{iteration_output}\n{divided_line}\n"

        Kratos.Logger.PrintInfo(to_print)

class OptimizationAnalysisTimeLogger:

    def __enter__(self):
        start_text = "Optimization Start"

        self.start_time = datetime.datetime.now()
        start_time_string = self.start_time.strftime("%Y-%m-%d %H:%M:%S")

        time_string = f"**   {start_time_string}   **"
        separator_string = len(time_string) * '*'
        center_string = f"**{start_text.center(len(str(time_string))-4)}**"
        final_string = f"{separator_string}\n{center_string}\n{time_string}\n{separator_string}"

        Kratos.Logger.PrintInfo(final_string)

    def __exit__(self, exit_type, exit_value, exit_traceback):

        end_text = "Optimization End"

        end_time = datetime.datetime.now()
        end_time_string = end_time.strftime("%Y-%m-%d %H:%M:%S")

        elapsed_time = end_time - self.start_time
        elapsed_time_string = str(datetime.timedelta(seconds=round(elapsed_time.total_seconds())))

        time_string = f"**   {end_time_string}  [Elapsed Time: {elapsed_time_string}]   **"
        separator_string = len(time_string) * '*'
        center_string = f"**{end_text.center(len(str(time_string))-4)}**"
        final_string = f"{separator_string}\n{center_string}\n{time_string}\n{separator_string}"

        Kratos.Logger.PrintInfo(final_string)

def DictLogger(title: str, data: dict):

    # First do the formatting and converting to strings
    for key, value in data.items():
        if isinstance(value, (float, int)):
            if abs(value) >= 1e6 or abs(value) < 1e-6:
                data[key] = "{:.6e}".format(value)
            else:
                data[key] = "{:.6f}".format(value)

    # Determine the maximum length of labels
    max_label_len = max(len(str(label)) for label in data.keys())
    max_value_len = max(len(str(value)) for value in data.values())

    title_len = len(str(title)) + 8

    # Calculate the row width
    if title_len > max_label_len + max_value_len + 7:
        row_width = title_len
        max_value_len = row_width - max_label_len - 7
    else:
        row_width = max_label_len + max_value_len + 7

    # Create format strings for the labels and values
    label_format = f"| {{:<{max_label_len}}} |"

    # Build the table
    table = '-' * row_width + "\n"
    table += f"|{str(title).center(row_width - 2)}|\n"
    table += '-' * row_width + "\n"

    # Add the data to the table
    for label, value in data.items():
        table += f"{label_format.format(label)} {value:>{max_value_len}} |\n"

    table += '-' * row_width

    # now print
    Kratos.Logger.PrintInfo(table)

def ListLogger(title: str, data: 'list[tuple[str, typing.Union[int, float, str]]]'):
    # First do the formatting and converting to strings
    sting_data: 'list[tuple[str, typing.Union[int, float, str]]]' = []
    for key, value in data:
        if isinstance(value, (float, int)):
            if abs(value) >= 1e6 or abs(value) < 1e-6:
                sting_data.append((key, "{:.6e}".format(value)))
            else:
                sting_data.append((key, "{:.6f}".format(value)))

    # Determine the maximum length of labels
    max_label_len = max(len(str(label)) for label, _  in data)
    max_value_len = max(len(str(value)) for _, value in data)

    title_len = len(str(title)) + 8

    # Calculate the row width
    if title_len > max_label_len + max_value_len + 7:
        row_width = title_len
        max_value_len = row_width - max_label_len - 7
    else:
        row_width = max_label_len + max_value_len + 7

    # Create format strings for the labels and values
    label_format = f"| {{:<{max_label_len}}} |"

    # Build the table
    table = '-' * row_width + "\n"
    table += f"|{str(title).center(row_width - 2)}|\n"
    table += '-' * row_width + "\n"

    # Add the data to the table
    for label, value in data:
        table += f"{label_format.format(label)} {value:>{max_value_len}} |\n"

    table += '-' * row_width

    # now print
    Kratos.Logger.PrintInfo(table)

def time_decorator(arg1=None, arg2=None, methodName=None):
    def inner_func(func):
        @wraps(func)
        def wrapper(self,*args, **kwargs):
            start_time = timer.perf_counter()
            if arg1:
                Kratos.Logger.Print(f"{func.__qualname__.split('.')[0]}::{func.__name__}: {arg1}")
            result = func(self,*args, **kwargs)
            end_time = timer.perf_counter()
            elapsed_time = end_time - start_time
            if arg2:
                Kratos.Logger.Print(f"{func.__qualname__.split('.')[0]}::{func.__name__}: {arg2} - [ Elapsed time: {datetime.timedelta(seconds=round(elapsed_time))}] \n")
            elif methodName:
                method_to_call = getattr(self, methodName)
                Kratos.Logger.Print(f"{func.__qualname__.split('.')[0]}::{func.__name__}: {method_to_call()} finished - [ Elapsed time: {datetime.timedelta(seconds=round(elapsed_time))}] \n")
            else:
                Kratos.Logger.Print(f"{func.__qualname__.split('.')[0]}::{func.__name__}: Finished - [ Elapsed time: {datetime.timedelta(seconds=round(elapsed_time))}] \n")
            return result
        return wrapper
    return inner_func