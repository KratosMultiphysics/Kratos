import KratosMultiphysics
from KratosMultiphysics import Logger
Logger.GetDefaultOutput().SetSeverity(Logger.Severity.INFO)
import KratosMultiphysics.DEMApplication

import main_script as Main

solution = Main.Solution()
solution.Run()
