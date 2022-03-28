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

class DEM3D_ContinuumTestSolution(KratosMultiphysics.DEMApplication.DEM_analysis_stage.DEMAnalysisStage, KratosUnittest.TestCase):

    def Initialize(self):
        super().Initialize()
        for node in self.spheres_model_part.Nodes:
            self.initial_normal_vel = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Z)


    def GetMainPath(self):
        return os.path.join(os.path.dirname(os.path.realpath(__file__)), "test_DEM_3D_continuum")

    def GetProblemNameWithPath(self):
        return os.path.join(self.main_path, self.DEM_parameters["problem_name"].GetString())


    def CheckValues(self, x_vel, z_vel, x_force, z_force, dem_pressure, z_elastic, shear, x_tangential):
        tol = 1.0e-8
        # DEM reference values
        x_vel_ref = 0.028907825348927448
        z_vel_ref = -0.8757276957864403
        x_force_ref = -26919.437972831598
        z_force_ref = 1970950.3578554934

        #FEM reference values
        dem_pressure_ref = 21566.85065708402
        z_elastic_ref = -273575.41245014494
        shear_ref = 362.391011482587
        x_tangential_ref = 6039.850191376444

        self.assertAlmostEqual(x_vel, x_vel_ref, delta=tol)
        self.assertAlmostEqual(z_vel, z_vel_ref, delta=tol)
        self.assertAlmostEqual(x_force, x_force_ref, delta=tol)
        self.assertAlmostEqual(z_force, z_force_ref, delta=tol)
        self.assertAlmostEqual(dem_pressure, dem_pressure_ref, delta=tol)
        self.assertAlmostEqual(z_elastic, z_elastic_ref, delta=tol)
        self.assertAlmostEqual(shear, shear_ref, delta=tol)
        self.assertAlmostEqual(x_tangential, x_tangential_ref, delta=tol)

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

        self.rigid_face_model_part.CreateSubModelPart('RigidFacePart')
        self.rigid_face_submpart = self.rigid_face_model_part.GetSubModelPart('RigidFacePart')
        rigid_face_part_nodes_id = [node.Id for node in self.rigid_face_model_part.Nodes]
        self.rigid_face_submpart.AddNodes(rigid_face_part_nodes_id)
        rigid_face_part_elements_id = [elem.Id for elem in self.rigid_face_model_part.Elements]
        self.rigid_face_submpart.AddElements(rigid_face_part_elements_id)
        rigid_face_part_conditions_id = [cond.Id for cond in self.rigid_face_model_part.Conditions]
        self.rigid_face_submpart.AddConditions(rigid_face_part_conditions_id)

        self.rigid_face_submpart.SetValue(DEM.COMPUTE_FORCES_ON_THIS_RIGID_ELEMENT, True)
        self.rigid_face_submpart.SetValue(DEM.RIGID_BODY_OPTION, True)
        self.rigid_face_submpart.SetValue(DEM.RIGID_BODY_MASS, 0.0)
        self.rigid_face_submpart.SetValue(DEM.RIGID_BODY_CENTER_OF_ROTATION_X, 0.0)
        self.rigid_face_submpart.SetValue(DEM.RIGID_BODY_CENTER_OF_ROTATION_Y, 0.0)
        self.rigid_face_submpart.SetValue(DEM.RIGID_BODY_CENTER_OF_ROTATION_Z, 0.0)
        self.rigid_face_submpart.SetValue(DEM.RIGID_BODY_INERTIAS_X, 0.0)
        self.rigid_face_submpart.SetValue(DEM.RIGID_BODY_INERTIAS_Y, 0.0)
        self.rigid_face_submpart.SetValue(DEM.RIGID_BODY_INERTIAS_Z, 0.0)

        orientation = KratosMultiphysics.Quaternion()
        orientation.X= 0.0
        orientation.Y = 0.0
        orientation.Z = 0.0
        orientation.W = 1.0
        self.rigid_face_submpart.SetValue(KratosMultiphysics.ORIENTATION, orientation)


class TestDEM3DContinuum(KratosUnittest.TestCase):

    def setUp(self):
        pass


    def test_DEM3D_continuum(self):
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "test_DEM_3D_continuum")
        parameters_file_name = os.path.join(path, "ProjectParametersDEM.json")
        model = KratosMultiphysics.Model()
        auxiliary_functions_for_tests.CreateAndRunStageInSelectedNumberOfOpenMPThreads(DEM3D_ContinuumTestSolution, model, parameters_file_name, auxiliary_functions_for_tests.GetHardcodedNumberOfThreads())

if __name__ == "__main__":
    Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
    KratosUnittest.main()
