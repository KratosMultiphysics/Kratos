import os
import KratosMultiphysics
from KratosMultiphysics import Logger
import KratosMultiphysics.DEMApplication as DEM
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.DEMApplication.DEM_analysis_stage
import auxiliary_functions_for_tests

Logger.GetDefaultOutput().SetSeverity(Logger.Severity.INFO)

this_working_dir_backup = os.getcwd()

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class DEM2D_RestitutionTestSolution(KratosMultiphysics.DEMApplication.DEM_analysis_stage.DEMAnalysisStage, KratosUnittest.TestCase):

    def Initialize(self):
        super().Initialize()
        for node in self.spheres_model_part.Nodes:
            self.initial_normal_vel = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y)

    @classmethod
    def GetMainPath(self):
        return os.path.join(os.path.dirname(os.path.realpath(__file__)), "DEM2D_restitution_tests_files")

    def GetProblemNameWithPath(self):
        return os.path.join(self.main_path, self.DEM_parameters["problem_name"].GetString())

    def Finalize(self):
        tolerance = 0.05
        for node in self.spheres_model_part.Nodes:
            final_vel = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y)
            Logger.PrintInfo("initial velocity:", self.initial_normal_vel)
            Logger.PrintInfo("final velocity:", final_vel)
            restitution_coefficient = -final_vel / self.initial_normal_vel
            Logger.PrintInfo("restitution_coefficient", restitution_coefficient)
            Logger.PrintInfo("ref:", self.coeff)
            Logger.PrintInfo("upper bound:", restitution_coefficient*tolerance)
            Logger.PrintInfo("lower bound:", restitution_coefficient/tolerance)
            self.assertAlmostEqual(self.coeff, restitution_coefficient, delta=tolerance)

        self.procedures.RemoveFoldersWithResults(str(self.main_path), str(self.problem_name), '')
        super().Finalize()

    def ReadModelParts(self, max_node_Id=0, max_elem_Id=0, max_cond_Id=0):
        properties = KratosMultiphysics.Properties(0)
        properties_walls = KratosMultiphysics.Properties(0)
        self.SetHardcodedProperties(properties, properties_walls)
        self.spheres_model_part.AddProperties(properties)
        self.rigid_face_model_part.AddProperties(properties_walls)

        translational_scheme = DEM.ForwardEulerScheme()
        translational_scheme.SetTranslationalIntegrationSchemeInProperties(properties, True)
        rotational_scheme = DEM.ForwardEulerScheme()
        rotational_scheme.SetRotationalIntegrationSchemeInProperties(properties, True)

        element_name = "CylinderParticle2D"
        DEM.PropertiesProxiesManager().CreatePropertiesProxies(self.spheres_model_part)

        coordinates = KratosMultiphysics.Array3()
        coordinates[0] = 0.0
        coordinates[1] = 0.00227
        coordinates[2] = 0.0
        radius = 0.0025
        self.creator_destructor.CreateSphericParticle(self.spheres_model_part, coordinates, properties, radius, element_name)

        for node in self.spheres_model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Y, -6.9)

        self.rigid_face_model_part.CreateNewNode(3, -0.01, 0.0, 0.0)
        self.rigid_face_model_part.CreateNewNode(4, 0.01, 0.0, 0.0)

        condition_name = "RigidEdge2D2N"
        self.rigid_face_model_part.CreateNewCondition(condition_name, 7, [3, 4], self.rigid_face_model_part.GetProperties()[0])

    @classmethod
    def SetHardcodedProperties(self, properties, properties_walls):
        self.coeff = 1.0

class DEM2D_RestitutionTestSolution_2(DEM2D_RestitutionTestSolution):

    def ReadMaterialsFile(self):
        materials_file_abs_path = os.path.join(self.GetMainPath(), "MaterialsDEM2.json")
        with open(materials_file_abs_path, 'r') as materials_file:
            self.DEM_material_parameters = KratosMultiphysics.Parameters(materials_file.read())

    @classmethod
    def SetHardcodedProperties(self, properties, properties_walls):
        self.coeff = 0.5

class TestDEM2DRestitution(KratosUnittest.TestCase):

    def setUp(self):
        pass

    @classmethod
    def test_DEM2D_restitution_1(self):
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "DEM2D_restitution_tests_files")
        parameters_file_name = os.path.join(path, "ProjectParametersDEM.json")
        model = KratosMultiphysics.Model()
        auxiliary_functions_for_tests.CreateAndRunStageInSelectedNumberOfOpenMPThreads(DEM2D_RestitutionTestSolution, model, parameters_file_name, auxiliary_functions_for_tests.GetHardcodedNumberOfThreads())

    @classmethod
    def test_DEM2D_restitution_2(self):
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "DEM2D_restitution_tests_files")
        parameters_file_name = os.path.join(path, "ProjectParametersDEM.json")
        materials_file_name = os.path.join(path, "MaterialsDEM2.json")
        model = KratosMultiphysics.Model()
        auxiliary_functions_for_tests.CreateAndRunStageInSelectedNumberOfOpenMPThreads(DEM2D_RestitutionTestSolution_2, model, parameters_file_name, auxiliary_functions_for_tests.GetHardcodedNumberOfThreads())

if __name__ == "__main__":
    Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
    KratosUnittest.main()
