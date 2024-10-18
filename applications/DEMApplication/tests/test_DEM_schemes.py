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

def SetHardcodedProperties(properties, properties_walls):

    properties[DEM.PARTICLE_DENSITY] = 4000.0
    properties[KratosMultiphysics.YOUNG_MODULUS] = 1.0e9
    properties[KratosMultiphysics.POISSON_RATIO] = 0.20

    properties_walls[DEM.COMPUTE_WEAR] = 0
    properties_walls[KratosMultiphysics.YOUNG_MODULUS] = 1.0e20
    properties_walls[KratosMultiphysics.POISSON_RATIO] = 0.23

class DEM3D_ForwardEulerTestSolution(KratosMultiphysics.DEMApplication.DEM_analysis_stage.DEMAnalysisStage, KratosUnittest.TestCase):

    def Initialize(self):
        super().Initialize()

    def GetMainPath(self):
        return os.path.join(os.path.dirname(os.path.realpath(__file__)), "test_schemes")

    def GetProblemNameWithPath(self):
        return os.path.join(self.main_path, self.DEM_parameters["problem_name"].GetString())

    def CheckValues(self, x_vel, dem_pressure):
        tol = 1.0e-10
        x_vel_ref = 0.02677117811730761
        self.assertAlmostEqual(x_vel, x_vel_ref, delta=tol)

        dem_pressure_ref = 21615.368565335117
        self.assertAlmostEqual(dem_pressure, dem_pressure_ref, delta=tol)

    def Finalize(self):
        node = self.spheres_model_part.GetNode(1)
        x_vel = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_X)

        node = self.rigid_face_model_part.GetNode(5)
        dem_pressure = node.GetSolutionStepValue(DEM.DEM_PRESSURE)

        self.CheckValues(x_vel, dem_pressure)
        self.procedures.RemoveFoldersWithResults(str(self.main_path), str(self.problem_name), '')
        super().Finalize()

    def ReadModelParts(self, max_node_Id=0, max_elem_Id=0, max_cond_Id=0):
        properties = KratosMultiphysics.Properties(0)
        properties_walls = KratosMultiphysics.Properties(0)
        SetHardcodedProperties(properties, properties_walls)
        self.spheres_model_part.AddProperties(properties)
        self.rigid_face_model_part.AddProperties(properties_walls)

        translational_scheme = DEM.ForwardEulerScheme()
        translational_scheme.SetTranslationalIntegrationSchemeInProperties(properties, True)
        rotational_scheme = DEM.ForwardEulerScheme()
        rotational_scheme.SetRotationalIntegrationSchemeInProperties(properties, True)

        element_name = "SphericContinuumParticle3D"
        DEM.PropertiesProxiesManager().CreatePropertiesProxies(self.spheres_model_part)

        coordinates = KratosMultiphysics.Array3()
        coordinates[0] = -1
        coordinates[1] = 0.0
        coordinates[2] = 0.0
        radius = 1
        self.creator_destructor.CreateSphericParticle(self.spheres_model_part, coordinates, properties, radius, element_name)

        coordinates = KratosMultiphysics.Array3()
        coordinates[0] = 0.95
        coordinates[1] = 0.0
        coordinates[2] = 0.0
        radius = 1
        self.creator_destructor.CreateSphericParticle(self.spheres_model_part, coordinates, properties, radius, element_name)

        for node in self.spheres_model_part.Nodes:
            node.SetSolutionStepValue(DEM.COHESIVE_GROUP, 1)

        node = self.spheres_model_part.GetNode(2)
        node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_X, 0.0)

        node = self.spheres_model_part.GetNode(1)
        node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_X, 0.1)

        self.rigid_face_model_part.CreateNewNode(3, -5, 5, -1.008)
        self.rigid_face_model_part.CreateNewNode(4, 5, 5, -1.008)

        self.rigid_face_model_part.CreateNewNode(5, -5, -5, -1.008)
        self.rigid_face_model_part.CreateNewNode(6, 5, -5, -1.008)

        condition_name = "RigidFace3D3N"
        self.rigid_face_model_part.CreateNewCondition(condition_name, 7, [5, 6, 3], self.rigid_face_model_part.GetProperties()[0])
        self.rigid_face_model_part.CreateNewCondition(condition_name, 8, [3, 6, 4], self.rigid_face_model_part.GetProperties()[0])

class DEM3D_TaylorTestSolution(DEM3D_ForwardEulerTestSolution):

    def CheckValues(self, x_vel, dem_pressure):
        tol = 1.0e-10
        x_vel_ref = 0.026581603224742834
        self.assertAlmostEqual(x_vel, x_vel_ref, delta=tol)

        dem_pressure_ref = 21598.828562997718
        self.assertAlmostEqual(dem_pressure, dem_pressure_ref, delta=tol)

class DEM3D_SymplecticTestSolution(DEM3D_ForwardEulerTestSolution):

    def CheckValues(self, x_vel, dem_pressure):
        tol = 1.0e-10
        x_vel_ref = 0.026396099600870276
        self.assertAlmostEqual(x_vel, x_vel_ref, delta=tol)

        dem_pressure_ref = 21582.363307819433
        self.assertAlmostEqual(dem_pressure, dem_pressure_ref, delta=tol)
class DEM3D_VerletTestSolution(DEM3D_ForwardEulerTestSolution):

    def CheckValues(self,x_vel, dem_pressure):
        tol = 1.0e-10
        x_vel_ref = 0.026480762680319958
        self.assertAlmostEqual(x_vel, x_vel_ref, delta=tol)

        dem_pressure_ref = 21609.179544197963
        self.assertAlmostEqual(dem_pressure, dem_pressure_ref, delta=tol)
class TestDEMSchemes(KratosUnittest.TestCase):

    def test_ForwardEuler(self):
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "test_schemes")
        parameters_file_name = os.path.join(path, "ProjectParametersDEM_ForwardEuler.json")
        model = KratosMultiphysics.Model()
        auxiliary_functions_for_tests.CreateAndRunStageInSelectedNumberOfOpenMPThreads(DEM3D_ForwardEulerTestSolution, model, parameters_file_name, auxiliary_functions_for_tests.GetHardcodedNumberOfThreads())

    def test_Taylor(self):
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "test_schemes")
        parameters_file_name = os.path.join(path, "ProjectParametersDEM_Taylor.json")
        model = KratosMultiphysics.Model()
        auxiliary_functions_for_tests.CreateAndRunStageInSelectedNumberOfOpenMPThreads(DEM3D_TaylorTestSolution, model, parameters_file_name, auxiliary_functions_for_tests.GetHardcodedNumberOfThreads())

    def test_Symplectic(self):
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "test_schemes")
        parameters_file_name = os.path.join(path, "ProjectParametersDEM_Symplectic.json")
        model = KratosMultiphysics.Model()
        auxiliary_functions_for_tests.CreateAndRunStageInSelectedNumberOfOpenMPThreads(DEM3D_SymplecticTestSolution, model, parameters_file_name, auxiliary_functions_for_tests.GetHardcodedNumberOfThreads())

    def test_Verlet(self):
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "test_schemes")
        parameters_file_name = os.path.join(path, "ProjectParametersDEM_Verlet.json")
        model = KratosMultiphysics.Model()
        auxiliary_functions_for_tests.CreateAndRunStageInSelectedNumberOfOpenMPThreads(DEM3D_VerletTestSolution, model, parameters_file_name, auxiliary_functions_for_tests.GetHardcodedNumberOfThreads())


if __name__ == "__main__":
    Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
    KratosUnittest.main()
