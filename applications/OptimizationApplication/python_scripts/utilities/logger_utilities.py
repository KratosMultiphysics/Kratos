# ==============================================================================
#  KratosOptimizationApplication
#
#  License:         BSD License
#                   license: OptimizationApplication/license.txt
#
#  Main authors:    Suneth Warnakulasuriya
#
# ==============================================================================

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