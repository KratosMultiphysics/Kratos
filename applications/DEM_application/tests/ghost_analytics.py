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

    def MakeAnalyticsMeasurements(self):
        self.face_watcher.MakeMeasurements(self.GetAnalyticFacesModelParts())
        #self.particle_watcher.MakeMeasurements()
        
        times = []
        neighbour_ids = []
        normal_relative_vel = []
        tangential_relative_vel = []
        my_id = []
        masses = []
        n_particles = []

        '''print("PRINTING IN PYTHON!!!!!!!!!!!!!!!!!!!!!!!!")
        print(my_id)
        print(neighbour_ids)
        print(masses)
        print(normal_relative_vel)
        print(tangential_relative_vel)
        print(n_particles)
        print(times)'''
        self.face_watcher.GetTotalFlux(times, n_particles, masses)
        #self.face_watcher.GetAllFacesData(self.GetAnalyticFacesModelParts(), times, neighbour_ids, masses, normal_relative_vel, tangential_relative_vel)
        #self.face_watcher.GetTimeStepsData(my_id, neighbour_ids, masses, normal_relative_vel, tangential_relative_vel)

        #print(self.GetAnalyticFacesModelParts())
        print("IN PYTHON, TIME====")
        print(times)
        print("IN PYTHON, N_PARTICLES====")
        print(n_particles)
        print("IN PYTHON, MASS====")
        print(masses)
        
    def Finalize(self):
        super(GhostsTestSolution, self).Finalize()
        self.procedures.RemoveFoldersWithResults(self.main_path, self.problem_name)

class GhostAnalytics(KratosUnittest.TestCase):

    def setUp(self):
        pass

    @classmethod
    def ghost_Analytics_1(self):
        GhostsTestSolution().Run()

    @classmethod
    def ghost_Analytics_2(self):
        pass

    def tearDown(self):
        file_to_remove = os.path.join("ghost_tests_files", "TimesPartialRelease")
        kratos_utils.DeleteFileIfExisting(GetFilePath(file_to_remove))


if __name__ == "__main__":
    Kratos.Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
    KratosUnittest.main()
