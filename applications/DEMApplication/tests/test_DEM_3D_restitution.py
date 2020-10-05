import os
import KratosMultiphysics
from KratosMultiphysics import Logger
Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
import KratosMultiphysics.DEMApplication as DEM
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.DEMApplication.DEM_analysis_stage

import KratosMultiphysics.kratos_utilities as kratos_utils

import auxiliary_functions_for_tests

this_working_dir_backup = os.getcwd()

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class DEM3D_RestitutionTestSolution(KratosMultiphysics.DEMApplication.DEM_analysis_stage.DEMAnalysisStage, KratosUnittest.TestCase):

    def Initialize(self):
        super().Initialize()
        for node in self.spheres_model_part.Nodes:
            self.initial_normal_vel = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Z)

    @classmethod
    def GetMainPath(self):
        return os.path.join(os.path.dirname(os.path.realpath(__file__)), "DEM3D_restitution_tests_files")

    def GetProblemNameWithPath(self):
        return os.path.join(self.main_path, self.DEM_parameters["problem_name"].GetString())

    def Finalize(self):
        tolerance = 1.0001
        for node in self.spheres_model_part.Nodes:
            final_vel = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Z)
            Logger.PrintInfo("initial velocity:",self.initial_normal_vel)
            Logger.PrintInfo("final velocity:",final_vel)
            restitution_coefficient = -final_vel / self.initial_normal_vel
            Logger.PrintInfo("restitution_coefficient",restitution_coefficient)
            Logger.PrintInfo("ref:",self.coeff)
            Logger.PrintInfo("upper bound:",restitution_coefficient*tolerance)
            Logger.PrintInfo("lower bound:",restitution_coefficient/tolerance)
            self.assertAlmostEqual(self.coeff, restitution_coefficient, delta=tolerance)
        super().Finalize()


    def ReadModelParts(self, max_node_Id=0, max_elem_Id=0, max_cond_Id=0):
        properties = KratosMultiphysics.Properties(0)
        properties_walls = KratosMultiphysics.Properties(0)
        self.SetHardcodedProperties(properties, properties_walls)
        self.spheres_model_part.AddProperties(properties)
        self.rigid_face_model_part.AddProperties(properties_walls)

        DiscontinuumConstitutiveLaw = getattr(DEM, properties[DEM.DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME])()
        DiscontinuumConstitutiveLaw.SetConstitutiveLawInProperties(properties, False)

        translational_scheme = DEM.ForwardEulerScheme()
        translational_scheme.SetTranslationalIntegrationSchemeInProperties(properties, True)
        rotational_scheme = DEM.ForwardEulerScheme()
        rotational_scheme.SetRotationalIntegrationSchemeInProperties(properties, True)

        element_name = "SphericParticle3D"
        DEM.PropertiesProxiesManager().CreatePropertiesProxies(self.spheres_model_part)

        coordinates = KratosMultiphysics.Array3()
        coordinates[0] = 0.0
        coordinates[1] = 0.0
        coordinates[2] = 0.00255
        radius = 0.0025
        self.creator_destructor.CreateSphericParticle(self.spheres_model_part, coordinates, properties, radius, element_name)

        for node in self.spheres_model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Z, -3.9)

        self.rigid_face_model_part.CreateNewNode(3, -0.01, 0.01, 0.0)
        self.rigid_face_model_part.CreateNewNode(4, 0.01, 0.01, 0.0)

        self.rigid_face_model_part.CreateNewNode(5, -0.01, -0.01, 0.0)
        self.rigid_face_model_part.CreateNewNode(6, 0.01, -0.01, 0.0)

        condition_name = "RigidFace3D3N"
        self.rigid_face_model_part.CreateNewCondition(condition_name, 7, [5, 6, 3], self.rigid_face_model_part.GetProperties()[0])
        self.rigid_face_model_part.CreateNewCondition(condition_name, 8, [3, 6, 4], self.rigid_face_model_part.GetProperties()[0])


    @classmethod
    def SetHardcodedProperties(self, properties, properties_walls):
        properties[DEM.PARTICLE_DENSITY] = 4000.0
        properties[KratosMultiphysics.YOUNG_MODULUS] = 3.8e11
        properties[KratosMultiphysics.POISSON_RATIO] = 0.23
        properties[DEM.STATIC_FRICTION] = 0.0
        properties[DEM.DYNAMIC_FRICTION] = 0.0
        properties[DEM.PARTICLE_COHESION] = 0.0
        properties[DEM.COEFFICIENT_OF_RESTITUTION] = 1.0
        self.coeff = 1.0963606640305437
        properties[KratosMultiphysics.PARTICLE_MATERIAL] = 1
        properties[DEM.ROLLING_FRICTION] = 0.0
        properties[DEM.DEM_CONTINUUM_CONSTITUTIVE_LAW_NAME] = "DEMContinuumConstitutiveLaw"
        properties[DEM.DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME] = "DEM_D_Hertz_viscous_Coulomb"

        properties_walls[DEM.STATIC_FRICTION] = 0.0
        properties_walls[DEM.DYNAMIC_FRICTION] = 0.0
        properties_walls[DEM.WALL_COHESION] = 0.0
        properties_walls[DEM.COMPUTE_WEAR] = 0
        properties_walls[DEM.SEVERITY_OF_WEAR] = 0.001
        properties_walls[DEM.IMPACT_WEAR_SEVERITY] = 0.001
        properties_walls[DEM.BRINELL_HARDNESS] = 200.0
        properties_walls[KratosMultiphysics.YOUNG_MODULUS] = 1.0e20
        properties_walls[KratosMultiphysics.POISSON_RATIO] = 0.23



class DEM3D_RestitutionTestSolution_2(DEM3D_RestitutionTestSolution):

    @classmethod
    def SetHardcodedProperties(self, properties, properties_walls):
        properties[DEM.PARTICLE_DENSITY] = 4000.0
        properties[KratosMultiphysics.YOUNG_MODULUS] = 3.8e11
        properties[KratosMultiphysics.POISSON_RATIO] = 0.23
        properties[DEM.STATIC_FRICTION] = 0.0
        properties[DEM.DYNAMIC_FRICTION] = 0.0
        properties[DEM.PARTICLE_COHESION] = 0.0
        properties[DEM.COEFFICIENT_OF_RESTITUTION] = 0.5
        self.coeff = 0.5472000544114178
        properties[KratosMultiphysics.PARTICLE_MATERIAL] = 1
        properties[DEM.ROLLING_FRICTION] = 0.0
        properties[DEM.DEM_CONTINUUM_CONSTITUTIVE_LAW_NAME] = "DEMContinuumConstitutiveLaw"
        properties[DEM.DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME] = "DEM_D_Hertz_viscous_Coulomb"

        properties_walls[DEM.STATIC_FRICTION] = 0.0
        properties_walls[DEM.DYNAMIC_FRICTION] = 0.0
        properties_walls[DEM.WALL_COHESION] = 0.0
        properties_walls[DEM.COMPUTE_WEAR] = 0
        properties_walls[DEM.SEVERITY_OF_WEAR] = 0.001
        properties_walls[DEM.IMPACT_WEAR_SEVERITY] = 0.001
        properties_walls[DEM.BRINELL_HARDNESS] = 200.0
        properties_walls[KratosMultiphysics.YOUNG_MODULUS] = 1.0e20
        properties_walls[KratosMultiphysics.POISSON_RATIO] = 0.23


class TestDEM3DRestitution(KratosUnittest.TestCase):

    def setUp(self):
        pass

    @classmethod
    def test_DEM3D_restitution_1(self):
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "DEM3D_restitution_tests_files")
        parameters_file_name = os.path.join(path, "ProjectParametersDEM.json")
        model = KratosMultiphysics.Model()
        auxiliary_functions_for_tests.CreateAndRunStageInSelectedNumberOfOpenMPThreads(DEM3D_RestitutionTestSolution, model, parameters_file_name, 1)

    @classmethod
    def test_DEM3D_restitution_2(self):
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "DEM3D_restitution_tests_files")
        parameters_file_name = os.path.join(path, "ProjectParametersDEM.json")
        model = KratosMultiphysics.Model()
        auxiliary_functions_for_tests.CreateAndRunStageInSelectedNumberOfOpenMPThreads(DEM3D_RestitutionTestSolution_2, model, parameters_file_name, 1)


    def tearDown(self):
        file_to_remove = os.path.join("DEM3D_restitution_tests_files", "TimesPartialRelease")
        kratos_utils.DeleteFileIfExisting(GetFilePath(file_to_remove))
        os.chdir(this_working_dir_backup)


if __name__ == "__main__":
    Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
    KratosUnittest.main()
