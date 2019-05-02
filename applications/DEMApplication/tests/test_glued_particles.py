import os
import KratosMultiphysics as Kratos
from KratosMultiphysics import Logger
Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
from KratosMultiphysics.DEMApplication import *
import KratosMultiphysics.KratosUnittest as KratosUnittest
import main_script

import KratosMultiphysics.kratos_utilities as kratos_utils

this_working_dir_backup = os.getcwd()

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

def CreateAndRunObjectInOneOpenMPThread(my_obj):
    omp_utils = Kratos.OpenMPUtils()
    if "OMP_NUM_THREADS" in os.environ:
        initial_number_of_threads = os.environ['OMP_NUM_THREADS']
        omp_utils.SetNumThreads(1)

    model = Kratos.Model()
    my_obj(model).Run()

    if "OMP_NUM_THREADS" in os.environ:
        omp_utils.SetNumThreads(int(initial_number_of_threads))

class GluedParticlesTestSolution(main_script.Solution):

    def GetParametersFileName(self):
        return os.path.join(self.GetMainPath(), "ProjectParametersDEM.json")

    @classmethod
    def GetMainPath(self):
        return os.path.join(os.path.dirname(os.path.realpath(__file__)), "glued_particles_tests_files")

    def GetProblemNameWithPath(self):
        return os.path.join(self.main_path, self.DEM_parameters["problem_name"].GetString())

    def FinalizeTimeStep(self, time):
        tolerance = 1e-4
        for node in self.spheres_model_part.Nodes:
            angular_velocity = node.GetSolutionStepValue(Kratos.ANGULAR_VELOCITY)
            if node.Id == 1:
                if time > 0.01:
                    self.CheckValue("Angular Velocity at time "+ str(time), angular_velocity[0], 2.0, tolerance)

                if time > 0.499999 and time < 0.5000001:
                    self.CheckValue("X Coordinate at time 0.5", node.X, -1.0,     tolerance)
                    self.CheckValue("Y Coordinate at time 0.5", node.Y,  0.6638445208099379, tolerance)
                    self.CheckValue("Z Coordinate at time 0.5", node.Z,  0.21679366644461678, tolerance)

                if time > 0.999999 and time < 1.0000001:
                    self.CheckValue("X Coordinate at time 1.0", node.X, -1.0, tolerance)
                    self.CheckValue("Y Coordinate at time 1.0", node.Y,  0.6359488, tolerance)
                    self.CheckValue("Z Coordinate at time 1.0", node.Z, -0.1657309642449, tolerance)

    @classmethod
    def CheckValue(self, explaining_string, value, expected_value, tolerance):
        if value > expected_value + tolerance or value < expected_value - tolerance:
            raise ValueError('Incorrect value for ' + explaining_string +': expected value was '+ str(expected_value) + ' but received ' + str(value))

    def Finalize(self):
        super(GluedParticlesTestSolution, self).Finalize()
        #self.procedures.RemoveFoldersWithResults(self.main_path, self.problem_name)


class TestGluedParticles(KratosUnittest.TestCase):

    def setUp(self):
        pass

    @classmethod
    def test_Glued_Particles_1(self):
        CreateAndRunObjectInOneOpenMPThread(GluedParticlesTestSolution)

    def tearDown(self):
        file_to_remove = os.path.join("glued_particles_tests_files", "TimesPartialRelease")
        kratos_utils.DeleteFileIfExisting(GetFilePath(file_to_remove))

        os.chdir(this_working_dir_backup)


if __name__ == "__main__":
    Kratos.Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
    KratosUnittest.main()
