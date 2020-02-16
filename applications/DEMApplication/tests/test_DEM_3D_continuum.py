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


class DEM3D_ContinuumTestSolution(KratosMultiphysics.DEMApplication.DEM_analysis_stage.DEMAnalysisStage):

    def Initialize(self):
        super(DEM3D_ContinuumTestSolution, self).Initialize()
        for node in self.spheres_model_part.Nodes:
            self.initial_normal_vel = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Z)

    @classmethod
    def GetMainPath(self):
        return os.path.join(os.path.dirname(os.path.realpath(__file__)), "test_DEM_3D_continuum")

    def GetProblemNameWithPath(self):
        return os.path.join(self.main_path, self.DEM_parameters["problem_name"].GetString())

    def CheckValues(self, x_vel, z_vel, x_force, z_force, dem_pressure, z_elastic, shear, x_tangential):
        tol = 1.0000+1.0e-3
        #DEM reference values
        x_vel_ref           =0.02043790
        z_vel_ref           =-0.8771226
        x_force_ref         =-25240.795
        z_force_ref         =1963237.70

        #FEM reference values
        dem_pressure_ref    =21558.5
        z_elastic_ref       =-273320
        shear_ref           =271.471
        x_tangential_ref    =4524.52


        if not (abs(x_vel_ref) < abs(x_vel*tol) and abs(x_vel_ref) > abs(x_vel/tol)):
            raise ValueError('Incorrect value for VELOCITY_X: expected value was '+ str(x_vel_ref) + ' but received ' + str(x_vel))

        if not (abs(z_vel_ref) < abs(z_vel*tol) and abs(z_vel_ref) > abs(z_vel/tol)):
            raise ValueError('Incorrect value for VELOCITY_Z: expected value was '+ str(z_vel_ref) + ' but received ' + str(z_vel))

        if not (abs(x_force_ref) < abs(x_force*tol) and abs(x_force_ref) > abs(x_force/tol)):
            raise ValueError('Incorrect value for FORCE_X: expected value was '+ str(x_force_ref) + ' but received ' + str(x_force))

        if not (abs(z_force_ref) < abs(z_force*tol) and abs(z_force_ref) > abs(z_force/tol)):
            raise ValueError('Incorrect value for FORCE_Z: expected value was '+ str(z_force_ref) + ' but received ' + str(z_force))

        if not (abs(dem_pressure_ref) < abs(dem_pressure*tol) and abs(dem_pressure_ref) > abs(dem_pressure/tol)):
            raise ValueError('Incorrect value for DEMPRESSURE: expected value was '+ str(dem_pressure_ref) + ' but received ' + str(dem_pressure))

        if not (abs(z_elastic_ref) < abs(z_elastic*tol) and abs(z_elastic_ref) > abs(z_elastic/tol)):
            raise ValueError('Incorrect value for ELASTIC_FORCE_Z: expected value was '+ str(z_elastic_ref) + ' but received ' + str(z_elastic))

        if not (abs(shear_ref) < abs(shear*tol) and abs(shear_ref) > abs(shear/tol)):
            raise ValueError('Incorrect value for SHEAR: expected value was '+ str(shear_ref) + ' but received ' + str(shear))

        if not (abs(x_tangential_ref) < abs(x_tangential*tol) and abs(x_tangential_ref) > abs(x_tangential/tol)):
            raise ValueError('Incorrect value for TANGENTIAL_ELASTIC_FORCE_Z: expected value was '+ str(x_tangential_ref) + ' but received ' + str(x_tangential))




    def Finalize(self):
        for node in self.spheres_model_part.Nodes:
            if node.Id == 1:
                x_vel = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_X)
                z_vel = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Z)
                x_force = node.GetSolutionStepValue(KratosMultiphysics.TOTAL_FORCES_X)
                z_force = node.GetSolutionStepValue(KratosMultiphysics.TOTAL_FORCES_Z)

        for node in self.rigid_face_model_part.Nodes:
            if node.Id == 5:
                dem_pressure = node.GetSolutionStepValue(DEM.DEM_PRESSURE)
                z_elastic = node.GetSolutionStepValue(DEM.ELASTIC_FORCES)[2]
                shear= node.GetSolutionStepValue(DEM.SHEAR_STRESS)
                x_tangential = node.GetSolutionStepValue(DEM.TANGENTIAL_ELASTIC_FORCES)[0]

        self.CheckValues(x_vel, z_vel, x_force, z_force, dem_pressure, z_elastic, shear, x_tangential)
        super(DEM3D_ContinuumTestSolution, self).Finalize()


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

        for node in self.spheres_model_part.Nodes:
            if node.Id == 2:
                node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_X, 0.0)
            if node.Id == 1:
                node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_X, 0.1)


        self.rigid_face_model_part.CreateNewNode(3, -5, 5, -1.008)
        self.rigid_face_model_part.CreateNewNode(4, 5, 5, -1.008)

        self.rigid_face_model_part.CreateNewNode(5, -5, -5, -1.008)
        self.rigid_face_model_part.CreateNewNode(6, 5, -5, -1.008)

        condition_name = "RigidFace3D3N"
        self.rigid_face_model_part.CreateNewCondition(condition_name, 7, [5, 6, 3], self.rigid_face_model_part.GetProperties()[0])
        self.rigid_face_model_part.CreateNewCondition(condition_name, 8, [3, 6, 4], self.rigid_face_model_part.GetProperties()[0])


    @classmethod
    def SetHardcodedProperties(self, properties, properties_walls):
        properties[DEM.PARTICLE_DENSITY] = 4000.0
        properties[KratosMultiphysics.YOUNG_MODULUS] = 1.0e9
        properties[KratosMultiphysics.POISSON_RATIO] = 0.20
        properties[DEM.FRICTION] = 0.5
        properties[DEM.PARTICLE_COHESION] = 0.0
        properties[DEM.COEFFICIENT_OF_RESTITUTION] = 0.5
        properties[KratosMultiphysics.PARTICLE_MATERIAL] = 1
        properties[DEM.ROLLING_FRICTION] = 0.0
        properties[DEM.DEM_CONTINUUM_CONSTITUTIVE_LAW_NAME] = "DEM_KDEM"
        properties[DEM.DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME] = "DEM_D_Hertz_viscous_Coulomb"
        properties[DEM.CONTACT_TAU_ZERO] = 0.5e6
        properties[DEM.CONTACT_SIGMA_MIN] = 1e6
        properties[DEM.CONTACT_INTERNAL_FRICC] = 1.0
        properties[DEM.ROTATIONAL_MOMENT_COEFFICIENT] = 0.0

        properties_walls[DEM.FRICTION] = 0.0
        properties_walls[DEM.WALL_COHESION] = 0.0
        properties_walls[DEM.COMPUTE_WEAR] = 0
        properties_walls[DEM.SEVERITY_OF_WEAR] = 0.001
        properties_walls[DEM.IMPACT_WEAR_SEVERITY] = 0.001
        properties_walls[DEM.BRINELL_HARDNESS] = 200.0
        properties_walls[KratosMultiphysics.YOUNG_MODULUS] = 1.0e20
        properties_walls[KratosMultiphysics.POISSON_RATIO] = 0.23


class TestDEM3DContinuum(KratosUnittest.TestCase):

    def setUp(self):
        pass

    @classmethod
    def test_DEM3D_continuum(self):
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "test_DEM_3D_continuum")
        parameters_file_name = os.path.join(path, "ProjectParametersDEM.json")
        model = KratosMultiphysics.Model()
        CreateAndRunStageInOneOpenMPThread(DEM3D_ContinuumTestSolution, model, parameters_file_name)

    def tearDown(self):
        file_to_remove = os.path.join("test_DEM_3D_continuum", "TimesPartialRelease")
        kratos_utils.DeleteFileIfExisting(GetFilePath(file_to_remove))
        os.chdir(this_working_dir_backup)


if __name__ == "__main__":
    Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
    KratosUnittest.main()
