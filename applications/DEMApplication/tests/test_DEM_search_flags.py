import os
import KratosMultiphysics
from KratosMultiphysics import Logger
Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
import KratosMultiphysics.DEMApplication as DEM
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.DEMApplication.DEM_analysis_stage

import auxiliary_functions_for_tests

this_working_dir_backup = os.getcwd()

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class DEM3D_SearchFlagTestSolution(KratosMultiphysics.DEMApplication.DEM_analysis_stage.DEMAnalysisStage, KratosUnittest.TestCase):
    '''
    Test definition:
    Two particles are created in two different locations and assigned velocities in order to make them collide. A FEM surface is also created between the two particles. With boths search flags set to FALSE, the test expect no collision between particles or against the middle wall.
    '''

    @classmethod
    def GetMainPath(self):
        return os.path.join(os.path.dirname(os.path.realpath(__file__)), "DEM3D_search_flags_tests_files")

    def GetProblemNameWithPath(self):
        return os.path.join(self.main_path, self.DEM_parameters["problem_name"].GetString())

    def Finalize(self):
        tolerance = 1.0e-5
        for node in self.spheres_model_part.Nodes:
            if node.Id == 1:
                final_vel = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Z)
                self.assertAlmostEqual(-3.9, final_vel, delta=tolerance)
            if node.Id == 2:
                final_vel = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Z)
                self.assertAlmostEqual(1.0, final_vel, delta=tolerance)
        self.procedures.RemoveFoldersWithResults(str(self.main_path), str(self.problem_name), '')
        super().Finalize()


    def ReadModelParts(self, max_node_Id=0, max_elem_Id=0, max_cond_Id=0):
        properties = KratosMultiphysics.Properties(0)
        properties_walls = KratosMultiphysics.Properties(0)
        self.spheres_model_part.AddProperties(properties)
        self.rigid_face_model_part.AddProperties(properties_walls)

        translational_scheme = DEM.ForwardEulerScheme()
        translational_scheme.SetTranslationalIntegrationSchemeInProperties(properties, True)
        rotational_scheme = DEM.ForwardEulerScheme()
        rotational_scheme.SetRotationalIntegrationSchemeInProperties(properties, True)

        element_name = "SphericParticle3D"
        DEM.PropertiesProxiesManager().CreatePropertiesProxies(self.spheres_model_part)

        coordinates = KratosMultiphysics.Array3()
        coordinates[0] = 0.0
        coordinates[1] = 0.0
        coordinates[2] = 0.0025002
        radius = 0.0025
        self.creator_destructor.CreateSphericParticle(self.spheres_model_part, coordinates, properties, radius, element_name)

        # Second particle to check search flag against particles
        coordinates[2] = -0.0026
        self.creator_destructor.CreateSphericParticle(self.spheres_model_part, coordinates, properties, radius, element_name)

        for node in self.spheres_model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Z, -3.9)
            if node.Id == 2:
                node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Z, 1)

        self.rigid_face_model_part.CreateNewNode(3, -0.01, 0.01, 0.0)
        self.rigid_face_model_part.CreateNewNode(4, 0.01, 0.01, 0.0)

        self.rigid_face_model_part.CreateNewNode(5, -0.01, -0.01, 0.0)
        self.rigid_face_model_part.CreateNewNode(6, 0.01, -0.01, 0.0)

        condition_name = "RigidFace3D3N"
        self.rigid_face_model_part.CreateNewCondition(condition_name, 7, [5, 6, 3], self.rigid_face_model_part.GetProperties()[0])
        self.rigid_face_model_part.CreateNewCondition(condition_name, 8, [3, 6, 4], self.rigid_face_model_part.GetProperties()[0])


class TestDEM3DSearchFlag(KratosUnittest.TestCase):

    def setUp(self):
        pass

    @classmethod
    def test_DEM3D_search(self):
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "DEM3D_search_flags_tests_files")
        parameters_file_name = os.path.join(path, "ProjectParametersDEM.json")
        model = KratosMultiphysics.Model()
        auxiliary_functions_for_tests.CreateAndRunStageInSelectedNumberOfOpenMPThreads(DEM3D_SearchFlagTestSolution, model, parameters_file_name, auxiliary_functions_for_tests.GetHardcodedNumberOfThreads())



if __name__ == "__main__":
    Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
    KratosUnittest.main()
