import os
import KratosMultiphysics as Kratos
from KratosMultiphysics import Logger
from KratosMultiphysics.DEMApplication import *
import KratosMultiphysics.KratosUnittest as KratosUnittest
import main_script

import KratosMultiphysics.kratos_utilities as kratos_utils

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)


class GhostsTestSolution(main_script.Solution):

    @classmethod
    def GetMainPath(self):
        return os.path.join(os.path.dirname(os.path.realpath(__file__)), "ghost_tests_files")

    def GetProblemNameWithPath(self):
        return os.path.join(self.main_path, self.DEM_parameters["problem_name"].GetString())

    def Finalize(self):
        super(GhostsTestSolution, self).Finalize()
        #self.procedures.RemoveFoldersWithResults(self.main_path, self.problem_name)
        pass

class GhostAnalytics(KratosUnittest.TestCase):

    def setUp(self):
        pass

    @classmethod
    def ghost_Analytics_1(self):
        GhostsTestSolution().Run()

    @classmethod
    def ghost_Analytics_2(self):
        pass

if __name__ == "__main__":
    #Kratos.Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
    KratosUnittest.main()
