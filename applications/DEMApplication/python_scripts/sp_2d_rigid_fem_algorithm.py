import KratosMultiphysics
from KratosMultiphysics.DEMApplication.DEM_analysis_stage import DEMAnalysisStage
from KratosMultiphysics.DEMApplication import mesh_creator_sphere
import KratosMultiphysics.DemStructuresCouplingApplication as DemFem

class DEMAnalysisStage2DSpRigidFem(DEMAnalysisStage):

    def __init__(self, model, project_parameters):
        super(DEMAnalysisStage2DSpRigidFem, self).__init__(model, project_parameters)

        # TEST NUMBER:
        # 1. CTW16, 2. CTW10, 3. CTW13, 4. Blind
        self.test_number = project_parameters["test_number"].GetInt()
        self.gmesh_with_inner_skin = project_parameters["gmesh_with_inner_skin"].GetBool()
        self.compute_skin_factor = project_parameters["compute_skin_factor"].GetDouble()

    def Initialize(self):
        super(DEMAnalysisStage2DSpRigidFem, self).Initialize()

        self.CreateSPMeasuringRingSubmodelpart(self.spheres_model_part)
        if self.gmesh_with_inner_skin:
            self.RebuildBlindSkinElements()

    def RebuildBlindSkinElements(self):

        self.PreUtilities.ResetSkinParticles(self.spheres_model_part)
        self.PreUtilities.ComputeSkin(self.spheres_model_part, self.compute_skin_factor)

        # These first two depend on the geometry of the test
        #inner_radius = 0.0381
        #outer_radius = 0.1524
        # These two depend on the particular mesh sizes
        #min_radius = 0.5 * 0.001
        #max_radius = 0.5 * 0.0015
        #center = KratosMultiphysics.Array3()
        #center[0] = center[1] = center[2] = 0.0
        #self.PreUtilities.ResetSkinParticles(self.spheres_model_part)
        #inner_skin_factor = 2.2
        #outer_skin_factor = 1.6
        #self.PreUtilities.SetSkinParticlesInnerBoundary(self.spheres_model_part, inner_radius, inner_skin_factor * min_radius)
        #self.PreUtilities.SetSkinParticlesOuterBoundaryBlind(self.spheres_model_part, outer_radius, center, outer_skin_factor * max_radius)

    def CreateSPMeasuringRingSubmodelpart(self, spheres_model_part):

        if not self.spheres_model_part.HasSubModelPart("RingSubmodelPart"):
            self.spheres_model_part.CreateSubModelPart('RingSubmodelPart')
        self.ring_submodelpart = self.spheres_model_part.GetSubModelPart('RingSubmodelPart')

        if self.test_number == 1: # CTW16
            zone_radius_to_measure_2d_sp = 0.015
        elif self.test_number == 2: # CTW10
            zone_radius_to_measure_2d_sp = 0.02
        elif self.test_number == 3: # CTW13
            zone_radius_to_measure_2d_sp = 0.017
        else: # Blind
            zone_radius_to_measure_2d_sp = 0.06

        nodes_in_zone_radius_list = []
        elements_in_zone_radius_list = []

        for element in self.spheres_model_part.Elements:
            node = element.GetNode(0)
            x = node.X
            y = node.Y

            if (x * x + y * y) < zone_radius_to_measure_2d_sp * zone_radius_to_measure_2d_sp:
                nodes_in_zone_radius_list.append(node.Id)
                elements_in_zone_radius_list.append(element.Id)

        self.ring_submodelpart.AddNodes(nodes_in_zone_radius_list)
        self.ring_submodelpart.AddElements(elements_in_zone_radius_list)

    def PrintResultsForGid(self, time):
        super(DEMAnalysisStage2DSpRigidFem, self).PrintResultsForGid(time)

        DemFem.DemStructuresCouplingUtilities().MarkBrokenSpheres(self.ring_submodelpart)
        self.SettingGeometricalSPValues()
        self.creator_destructor.MarkParticlesForErasingGivenCylinder(self.ring_submodelpart, self.center, self.axis, self.radius)

        DemFem.DemStructuresCouplingUtilities().ComputeSandProductionWithDepthFirstSearchNonRecursiveImplementation(self.ring_submodelpart, self.rigid_face_model_part, self.time)
        DemFem.DemStructuresCouplingUtilities().ComputeSandProduction(self.ring_submodelpart, self.rigid_face_model_part, self.time)

    def SettingGeometricalSPValues(self):

        self.center = KratosMultiphysics.Array3()
        self.center[0] = 0; # self.sp_parameters["problem_data"]["center"][0].GetDouble()
        self.center[1] = 0; # self.sp_parameters["problem_data"]["center"][1].GetDouble()
        self.center[2] = 0; # self.sp_parameters["problem_data"]["center"][2].GetDouble()
        self.axis = KratosMultiphysics.Array3()
        self.axis[0] = 0; # self.sp_parameters["problem_data"]["axis"][0].GetDouble()
        self.axis[1] = 0; # self.sp_parameters["problem_data"]["axis"][1].GetDouble()
        self.axis[2] = 1; # self.sp_parameters["problem_data"]["axis"][2].GetDouble()

        # Zone radius to delete SP
        if self.test_number == 1: # CTW16
            self.radius = 0.0015
        elif self.test_number == 2: # CTW10
            self.radius = 0.005
        elif self.test_number == 3: # CTW13
            self.radius = 0.0025
        else: # Blind
            self.radius = 0.015

    def AdditionalFinalizeOperations(self):

        spheres_mp_filename_post = self.problem_name + 'DEM_Post'

        if self.write_mdpa_from_results:
            mesh_creator_sphere.WriteSphereMdpaFromResults(self.problem_name + 'DEM', self.main_path, spheres_mp_filename_post, self.spheres_model_part)

if __name__ == "__main__":
    DEMAnalysisStage2DSpRigidFem(model, project_parameters).Run()
