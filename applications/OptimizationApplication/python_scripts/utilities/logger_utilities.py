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
import os
import sys
import datetime
import KratosMultiphysics as Kratos

def GetTerminalWidth():
    try:
        if sys.stdout.isatty():
            return os.get_terminal_size().columns
    except OSError:
        pass
    return 80

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
    def __init__(self, topic: str, entry_msg: str , exit_msg: str):
        self.topic = topic
        self.entry_msg = entry_msg
        self.exit_msg = exit_msg
        self.start_time = None

    def __enter__(self):
        if self.entry_msg is not None:
            Kratos.Logger.PrintInfo(self.topic, self.entry_msg)
        self.start_time = timer.time()

    def __exit__(self, exit_type, exit_value, exit_traceback):
        if self.exit_msg is not None:
            elapsed_time = timer.time() - self.start_time
            Kratos.Logger.PrintInfo(self.topic, "{:s} - [ Elapsed time: {:s} ]".format(self.exit_msg, str(datetime.timedelta(seconds=round(elapsed_time))))+"\n")

class OptimizationAlgorithmTimeLogger:
    """@brief A context responsible for outputting execution times."""
    def __init__(self, optimizer_name: str, optimization_itr: int):
        self.start_time = None
        self.optimization_itr = optimization_itr
        self.optimizer_name = optimizer_name

    def __enter__(self):
        self.start_time = timer.time()

    def __exit__(self, exit_type, exit_value, exit_traceback):

        terminal_width = GetTerminalWidth()
        border_symbol = '@'
        empty_space_percentage = 0.1

        # Calculate the width of the content (excluding borders and empty spaces)
        content_width = terminal_width - 2 - int(terminal_width * empty_space_percentage) * 2

        current_time = timer.time()
        elapsed_time = current_time - self.start_time
        elapsed_time_string = str(datetime.timedelta(seconds=round(elapsed_time)))

        iteration_text = f"{self.optimizer_name} EOF Iteration {self.optimization_itr}"
        iteration_output = f"{border_symbol}{' ' * int(content_width * empty_space_percentage)}{iteration_text} [Elapsed Time: {elapsed_time_string}]{' ' * int(content_width * empty_space_percentage)}{border_symbol}"

        total_width = terminal_width

        left_width = int(total_width * 0.1)
        center_width = int(total_width * 0.8)
        right_width = total_width - left_width - center_width

        left_part = ' ' * left_width
        center_part = '#' * center_width
        right_part = ' ' * right_width

        divided_line = f"{left_part}{center_part}{right_part}"

        to_print = f"{divided_line}\n{iteration_output.center(terminal_width)}\n{divided_line}"

        Kratos.Logger.Print(to_print)

class OptimizationAnalysisTimeLogger:

    def __enter__(self):
        terminal_width = GetTerminalWidth()
        border_symbol = '#'
        start_text = "Optimization Start"

        border = border_symbol * terminal_width
        start_centered_text = start_text.center(terminal_width - 2, ' ')
        separator = border_symbol * terminal_width

        self.start_time = datetime.datetime.now()
        start_time_string = self.start_time.strftime("%Y-%m-%d %H:%M:%S")

        start_output = f"{border}\n{border_symbol}{start_centered_text}{border_symbol}\n{start_time_string.center(terminal_width - 2)}\n{separator}"

        Kratos.Logger.Print(start_output)

    def __exit__(self, exit_type, exit_value, exit_traceback):

        terminal_width = GetTerminalWidth()
        border_symbol = '#'
        end_text = "Optimization End"

        border = border_symbol * terminal_width
        end_centered_text = end_text.center(terminal_width - 2, ' ')
        separator = border_symbol * terminal_width

        end_time = datetime.datetime.now()
        end_time_string = end_time.strftime("%Y-%m-%d %H:%M:%S")

        elapsed_time = end_time - self.start_time
        elapsed_time_string = str(datetime.timedelta(seconds=round(elapsed_time.total_seconds())))

        test_string = end_time_string + "  [Elapsed Time: " + elapsed_time_string + "]"

        end_output = f"{border}\n{border_symbol}{end_centered_text}{border_symbol}\n{test_string.center(terminal_width)}\n{separator}"

        Kratos.Logger.Print(end_output)
