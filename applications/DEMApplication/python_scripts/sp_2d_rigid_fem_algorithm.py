import KratosMultiphysics
from KratosMultiphysics.DEMApplication.DEM_analysis_stage import DEMAnalysisStage
from KratosMultiphysics.DEMApplication import mesh_creator_sphere
import KratosMultiphysics.DemStructuresCouplingApplication as DemFem

class DEMAnalysisStage2DSpRigidFem(DEMAnalysisStage):

    def __init__(self, model, project_parameters):
        super(DEMAnalysisStage2DSpRigidFem, self).__init__(model, project_parameters)

        # TEST NUMBER:
        # 1. CTW16, 2. CTW10, 3. CTW13, 4. CTW12, 5.Blind
        self.test_number = project_parameters["test_number"].GetInt()
        self.gmesh_with_inner_skin = project_parameters["gmesh_with_inner_skin"].GetBool()
        self.compute_skin_factor = project_parameters["compute_skin_factor"].GetDouble()
        self.inner_mesh_diameter = project_parameters["inner_mesh_diameter"].GetDouble() # This depends on the particular GiD mesh (diameter of the finer mesh)
        self.outer_mesh_diameter = project_parameters["outer_mesh_diameter"].GetDouble() # This depends on the particular GiD mesh (diameter of the coarser mesh)
        # The two values that follow may depend on the GiD mesh used. The higher the value, the more skin particles
        self.inner_skin_factor = project_parameters["inner_skin_factor"].GetDouble() # 2.4
        self.outer_skin_factor = project_parameters["outer_skin_factor"].GetDouble() # 0.8

    def Initialize(self):
        super(DEMAnalysisStage2DSpRigidFem, self).Initialize()

        self.SettingGeometricalSPValues()
        self.CreateSPMeasuringRingSubmodelpart(self.spheres_model_part)
        #self.RebuildSkinElements()

    def SettingGeometricalSPValues(self):

        self.center = KratosMultiphysics.Array3()
        self.center[0] = 0; # self.sp_parameters["problem_data"]["center"][0].GetDouble()
        self.center[1] = 0; # self.sp_parameters["problem_data"]["center"][1].GetDouble()
        self.center[2] = 0; # self.sp_parameters["problem_data"]["center"][2].GetDouble()
        self.axis = KratosMultiphysics.Array3()
        self.axis[0] = 0; # self.sp_parameters["problem_data"]["axis"][0].GetDouble()
        self.axis[1] = 0; # self.sp_parameters["problem_data"]["axis"][1].GetDouble()
        self.axis[2] = 1; # self.sp_parameters["problem_data"]["axis"][2].GetDouble()

        self.outer_radius = 0.0508 # For the CTWs
        if self.test_number == 1: # CTW16
            self.inner_radius = 0.00381
            self.zone_radius_to_measure_2d_sp = 0.015
            self.radius_to_delete_sp = 0.0015
        elif self.test_number == 2: # CTW10
            self.inner_radius = 0.0127
            self.zone_radius_to_measure_2d_sp = 0.02
            self.radius_to_delete_sp = 0.005
        elif self.test_number == 3: # CTW13
            self.inner_radius = 0.00635
            self.zone_radius_to_measure_2d_sp = 0.017
            self.radius_to_delete_sp = 0.0025
        elif self.test_number == 4: # CTW12
            self.inner_radius = 0.008
            self.zone_radius_to_measure_2d_sp = 0.018
            self.radius_to_delete_sp = 0.003
        else: # Blind
            self.inner_radius = 0.0381
            self.zone_radius_to_measure_2d_sp = 0.06
            self.radius_to_delete_sp = 0.015
            self.outer_radius = 0.1524

    def RebuildSkinElements(self):

        self.PreUtilities.ResetSkinParticles(self.spheres_model_part)

        if self.gmesh_with_inner_skin:
            self.PreUtilities.ComputeSkinIncludingInnerVoids(self.spheres_model_part, self.compute_skin_factor)
        else:
            self.PreUtilities.SetSkinParticlesInnerBoundary(self.spheres_model_part, self.inner_radius, self.inner_skin_factor * 0.5 * self.inner_mesh_diameter)
            if self.test_number < 5: # CTWs
                self.PreUtilities.SetSkinParticlesOuterBoundary(self.spheres_model_part, self.outer_radius, self.outer_skin_factor * 0.5 * self.outer_mesh_diameter)
            else: # Blind
                self.PreUtilities.SetSkinParticlesOuterBoundaryBlind(self.spheres_model_part, self.outer_radius, self.center, self.outer_skin_factor * 0.5 * self.outer_mesh_diameter)

    def CreateSPMeasuringRingSubmodelpart(self, spheres_model_part):

        if not self.spheres_model_part.HasSubModelPart("RingSubmodelPart"):
            self.spheres_model_part.CreateSubModelPart('RingSubmodelPart')
        self.ring_submodelpart = self.spheres_model_part.GetSubModelPart('RingSubmodelPart')

        nodes_in_zone_radius_list = []
        elements_in_zone_radius_list = []

        for element in self.spheres_model_part.Elements:
            node = element.GetNode(0)
            x = node.X
            y = node.Y

            if (x * x + y * y) < self.zone_radius_to_measure_2d_sp * self.zone_radius_to_measure_2d_sp:
                nodes_in_zone_radius_list.append(node.Id)
                elements_in_zone_radius_list.append(element.Id)

        self.ring_submodelpart.AddNodes(nodes_in_zone_radius_list)
        self.ring_submodelpart.AddElements(elements_in_zone_radius_list)

    def PrintResultsForGid(self, time):
        super(DEMAnalysisStage2DSpRigidFem, self).PrintResultsForGid(time)

        DemFem.DemStructuresCouplingUtilities().MarkBrokenSpheres(self.ring_submodelpart)
        self.creator_destructor.MarkParticlesForErasingGivenCylinder(self.ring_submodelpart, self.center, self.axis, self.radius_to_delete_sp)

        DemFem.DemStructuresCouplingUtilities().ComputeSandProductionWithDepthFirstSearchNonRecursiveImplementation(self.ring_submodelpart, self.rigid_face_model_part, self.time)
        DemFem.DemStructuresCouplingUtilities().ComputeSandProduction(self.ring_submodelpart, self.rigid_face_model_part, self.time)

    def AdditionalFinalizeOperations(self):

        spheres_mp_filename_post = self.problem_name + 'DEM_Post'

        if self.write_mdpa_from_results:
            mesh_creator_sphere.WriteSphereMdpaFromResults(self.problem_name + 'DEM', self.main_path, spheres_mp_filename_post, self.spheres_model_part)

if __name__ == "__main__":
    DEMAnalysisStage2DSpRigidFem(model, project_parameters).Run()
