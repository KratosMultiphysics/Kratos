import KratosMultiphysics
from KratosMultiphysics.DEMApplication.DEM_analysis_stage import DEMAnalysisStage
from KratosMultiphysics.DEMApplication import mesh_creator_sphere_2D
import KratosMultiphysics.DemStructuresCouplingApplication as DemFem

class DEMAnalysisStage2DSpRigidFem(DEMAnalysisStage):

    def __init__(self, model, project_parameters):
        super(DEMAnalysisStage2DSpRigidFem, self).__init__(model, project_parameters)

    def Initialize(self):
        super(DEMAnalysisStage2DSpRigidFem, self).Initialize()

        self.CreateSPMeasuringRingSubmodelpart(self.spheres_model_part)

    def CreateSPMeasuringRingSubmodelpart(self, spheres_model_part):

        if not self.spheres_model_part.HasSubModelPart("RingSubmodelPart"):
            self.spheres_model_part.CreateSubModelPart('RingSubmodelPart')
        self.ring_submodelpart = self.spheres_model_part.GetSubModelPart('RingSubmodelPart')

        zone_radius_to_measure_2d_sp = 0.02 # In meters. This is a little bit less than half the CTW16 radius
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

        self.radius = 0
        self.test_number = 1
        if self.test_number == 1:
            self.radius = 0.0036195; #0.01; #0.0036195; #95% of the real hole. CTW16 specimen
        elif self.test_number == 2:
            self.radius = 0.012065; #95% of the real hole. CTW10 specimen
        elif self.test_number == 3:
            self.radius = 0.036195; #95% of the real hole. Blind Test

    def AdditionalFinalizeOperations(self):

        spheres_mp_filename_post = self.problem_name + 'DEM_Post'

        if self.write_mdpa_from_results:
            mesh_creator_sphere_2D.WriteSphereMdpaFromResults(self.problem_name + 'DEM', self.main_path, spheres_mp_filename_post, self.file_msh, self.file_res, self.post_path)

if __name__ == "__main__":
    DEMAnalysisStage2DSpRigidFem(model, project_parameters).Run()
