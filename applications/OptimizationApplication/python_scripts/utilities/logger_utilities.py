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
    """@brief A context responsible for outputing execution times."""
    def __init__(self, topic: str, entry_msg: str, exit_msg: str):
        self.topic = topic
        self.entry_msg = entry_msg
        self.exit_msg = exit_msg

    def __enter__(self):
        Kratos.Logger.PrintInfo(self.topic, self.entry_msg)
        self.start_time = timer.time()

    def __exit__(self, exit_type, exit_value, exit_traceback):
        Kratos.Logger.PrintInfo(self.topic, "{:s} - [ Elapsed time: {:s} ]".format(self.exit_msg, str(datetime.timedelta(seconds=round(timer.time() - self.start_time)))))