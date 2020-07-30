import os
import KratosMultiphysics
from KratosMultiphysics import Logger
Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
import KratosMultiphysics.DEMApplication as DEM
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.DEMApplication.DEM_analysis_stage

import KratosMultiphysics.kratos_utilities as kratos_utils

this_working_dir_backup = os.getcwd()

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

def CreateAndRunStageInOneOpenMPThread(my_obj, model, parameters_file_name):
    omp_utils = KratosMultiphysics.OpenMPUtils()
    if "OMP_NUM_THREADS" in os.environ:
        initial_number_of_threads = os.environ['OMP_NUM_THREADS']
        omp_utils.SetNumThreads(1)

    with open(parameters_file_name,'r') as parameter_file:
        project_parameters = KratosMultiphysics.Parameters(parameter_file.read())

    my_obj(model, project_parameters).Run()

    if "OMP_NUM_THREADS" in os.environ:
        omp_utils.SetNumThreads(int(initial_number_of_threads))

class DEM2D_RestitutionTestSolution(KratosMultiphysics.DEMApplication.DEM_analysis_stage.DEMAnalysisStage):

    def Initialize(self):
        super(DEM2D_RestitutionTestSolution, self).Initialize()
        for node in self.spheres_model_part.Nodes:
            self.initial_normal_vel = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y)

    @classmethod
    def GetMainPath(self):
        return os.path.join(os.path.dirname(os.path.realpath(__file__)), "DEM2D_restitution_tests_files")

    def GetProblemNameWithPath(self):
        return os.path.join(self.main_path, self.DEM_parameters["problem_name"].GetString())

    @classmethod
    def CheckRestitution(self, reference, restitution_coefficient, tolerance):
        if not (reference < restitution_coefficient*tolerance and reference > restitution_coefficient/tolerance):
            raise ValueError('Incorrect value for COEFFICIENT_OF_RESTITUTION: expected value was '+ str(reference) + ' but received ' + str(restitution_coefficient))

    def Finalize(self):
        tolerance = 1.0+1.0e-4
        for node in self.spheres_model_part.Nodes:
            final_vel = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y)
            Logger.PrintInfo("initial velocity:",self.initial_normal_vel)
            Logger.PrintInfo("final velocity:",final_vel)
            restitution_coefficient = -final_vel / self.initial_normal_vel
            Logger.PrintInfo("restitution_coefficient",restitution_coefficient)
            Logger.PrintInfo("ref:",self.coeff)
            Logger.PrintInfo("upper bound:",restitution_coefficient*tolerance)
            Logger.PrintInfo("lower bound:",restitution_coefficient/tolerance)
            self.CheckRestitution(self.coeff, restitution_coefficient, tolerance)
        super(DEM2D_RestitutionTestSolution, self).Finalize()


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

        element_name = "CylinderParticle2D"
        DEM.PropertiesProxiesManager().CreatePropertiesProxies(self.spheres_model_part)

        coordinates = KratosMultiphysics.Array3()
        coordinates[0] = 0.0
        coordinates[1] = 0.00255
        radius = 0.0025
        self.creator_destructor.CreateSphericParticle(self.spheres_model_part, coordinates, properties, radius, element_name)

        for node in self.spheres_model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Y, -3.9)

        self.rigid_face_model_part.CreateNewNode(3, -0.01, 0.0, 0.0)
        self.rigid_face_model_part.CreateNewNode(4, 0.01, 0.0, 0.0)

        condition_name = "RigidEdge2D2N"
        self.rigid_face_model_part.CreateNewCondition(condition_name, 7, [3, 4], self.rigid_face_model_part.GetProperties()[0])

    @classmethod
    def SetHardcodedProperties(self, properties, properties_walls):
        properties[DEM.PARTICLE_DENSITY] = 4000.0
        properties[KratosMultiphysics.YOUNG_MODULUS] = 3.8e11
        properties[KratosMultiphysics.POISSON_RATIO] = 0.23
        properties[DEM.FRICTION] = 0.0
        properties[DEM.PARTICLE_COHESION] = 0.0
        properties[DEM.COEFFICIENT_OF_RESTITUTION] = 1.0
        self.coeff = 1.0985200123566359      # reference value at dt=1.0e-6
        properties[KratosMultiphysics.PARTICLE_MATERIAL] = 1
        properties[DEM.ROLLING_FRICTION] = 0.0
        properties[DEM.DEM_CONTINUUM_CONSTITUTIVE_LAW_NAME] = "DEMContinuumConstitutiveLaw"
        properties[DEM.DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME] = "DEM_D_Hertz_viscous_Coulomb"

        properties_walls[DEM.FRICTION] = 0.0
        properties_walls[DEM.WALL_COHESION] = 0.0
        properties_walls[DEM.COMPUTE_WEAR] = 0
        properties_walls[DEM.SEVERITY_OF_WEAR] = 0.001
        properties_walls[DEM.IMPACT_WEAR_SEVERITY] = 0.001
        properties_walls[DEM.BRINELL_HARDNESS] = 200.0
        properties_walls[KratosMultiphysics.YOUNG_MODULUS] = 1.0e20
        properties_walls[KratosMultiphysics.POISSON_RATIO] = 0.23



class DEM2D_RestitutionTestSolution_2(DEM2D_RestitutionTestSolution):

    @classmethod
    def SetHardcodedProperties(self, properties, properties_walls):
        properties[DEM.PARTICLE_DENSITY] = 4000.0
        properties[KratosMultiphysics.YOUNG_MODULUS] = 3.8e11
        properties[KratosMultiphysics.POISSON_RATIO] = 0.23
        properties[DEM.FRICTION] = 0.0
        properties[DEM.PARTICLE_COHESION] = 0.0
        properties[DEM.COEFFICIENT_OF_RESTITUTION] = 0.5
        self.coeff = 0.5455811757241611      # reference value at dt=1.0e-6
        properties[KratosMultiphysics.PARTICLE_MATERIAL] = 1
        properties[DEM.ROLLING_FRICTION] = 0.0
        properties[DEM.DEM_CONTINUUM_CONSTITUTIVE_LAW_NAME] = "DEMContinuumConstitutiveLaw"
        properties[DEM.DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME] = "DEM_D_Hertz_viscous_Coulomb"

        properties_walls[DEM.FRICTION] = 0.0
        properties_walls[DEM.WALL_COHESION] = 0.0
        properties_walls[DEM.COMPUTE_WEAR] = 0
        properties_walls[DEM.SEVERITY_OF_WEAR] = 0.001
        properties_walls[DEM.IMPACT_WEAR_SEVERITY] = 0.001
        properties_walls[DEM.BRINELL_HARDNESS] = 200.0
        properties_walls[KratosMultiphysics.YOUNG_MODULUS] = 1.0e20
        properties_walls[KratosMultiphysics.POISSON_RATIO] = 0.23


class TestDEM2DRestitution(KratosUnittest.TestCase):

    def setUp(self):
        pass

    @classmethod
    def test_DEM2D_restitution_1(self):
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "DEM2D_restitution_tests_files")
        parameters_file_name = os.path.join(path, "ProjectParametersDEM.json")
        model = KratosMultiphysics.Model()
        CreateAndRunStageInOneOpenMPThread(DEM2D_RestitutionTestSolution, model, parameters_file_name)

    @classmethod
    def test_DEM2D_restitution_2(self):
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "DEM2D_restitution_tests_files")
        parameters_file_name = os.path.join(path, "ProjectParametersDEM.json")
        model = KratosMultiphysics.Model()
        CreateAndRunStageInOneOpenMPThread(DEM2D_RestitutionTestSolution_2, model, parameters_file_name)


    def tearDown(self):
        file_to_remove = os.path.join("DEM2D_restitution_tests_files", "TimesPartialRelease")
        kratos_utils.DeleteFileIfExisting(GetFilePath(file_to_remove))
        os.chdir(this_working_dir_backup)


if __name__ == "__main__":
    Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
    KratosUnittest.main()
