import KratosMultiphysics
from KratosMultiphysics.DEMApplication.DEM_analysis_stage import DEMAnalysisStage
from KratosMultiphysics.DEMApplication import mesh_creator_sphere
import KratosMultiphysics.DemStructuresCouplingApplication as DemFem
import math

class DEMAnalysisStage2DSpRigidFem(DEMAnalysisStage):

    def __init__(self, model, project_parameters):
        super(DEMAnalysisStage2DSpRigidFem, self).__init__(model, project_parameters)

        sp_project_parameters_file_name = "sp_2d_rigid_fem_parameters.json"

        with open(sp_project_parameters_file_name,'r') as parameters_file:
            sp_project_parameters = KratosMultiphysics.Parameters(parameters_file.read())

        # TEST NUMBER:
        # 1. CTW16, 2. CTW10, 3. CTW13, 4. CTW12, 5.Blind
        self.test_number = sp_project_parameters["test_number"].GetInt()
        self.inner_mesh_diameter = sp_project_parameters["inner_mesh_diameter"].GetDouble() # This depends on the particular GiD mesh (diameter of the finer mesh)
        self.outer_mesh_diameter = sp_project_parameters["outer_mesh_diameter"].GetDouble() # This depends on the particular GiD mesh (diameter of the coarser mesh)
        # The two values that follow may depend on the GiD mesh used. The higher the value, the more skin particles
        self.inner_skin_factor = sp_project_parameters["inner_skin_factor"].GetDouble() # 2.4
        self.outer_skin_factor = sp_project_parameters["outer_skin_factor"].GetDouble() # 0.8
        self.respect_preprocessor_marked_skin = sp_project_parameters["respect_preprocessor_marked_skin"].GetBool()

        self.automatic_skin_computation = project_parameters["AutomaticSkinComputation"].GetBool()

    def Initialize(self):
        super(DEMAnalysisStage2DSpRigidFem, self).Initialize()

        # self.SettingGeometricalSPValues()
        # self.CreateSPMeasuringRingSubmodelpart(self.spheres_model_part)
        # self.SetSkinManually()

        # sandstone_target_porosity, actual_porosity, sp_porosity_multiplier = self.ComputePorosityParameters(self.spheres_model_part)
        # porosity_message = "\nPorosity in sandstones should be around %.2f and the obtained value is %.2f" % (sandstone_target_porosity, actual_porosity) + \
        # "\nSand production results should therefore be multiplied by a factor of %.2f\n" % sp_porosity_multiplier
        # print(porosity_message)
        # with open('sp_porosity.txt', 'w') as poro_file:
        #     poro_file.write(porosity_message)

        from KratosMultiphysics.DEMApplication.multiaxial_control_module_generalized_2d_utility import MultiaxialControlModuleGeneralized2DUtility
        self.multiaxial_control_module = MultiaxialControlModuleGeneralized2DUtility(self.spheres_model_part, self.rigid_face_model_part)
        self.multiaxial_control_module.ExecuteInitialize()

        self.times = []
        self.reaction_x = []
        self.reaction_y = []
        self.reaction_z = []
        self.smoothed_reaction_x = []
        self.smoothed_reaction_y = []
        self.smoothed_reaction_z = []
        self.elastic_reaction_x = []
        self.elastic_reaction_y = []
        self.elastic_reaction_z = []
        self.smoothed_elastic_reaction_x = []
        self.smoothed_elastic_reaction_y = []
        self.smoothed_elastic_reaction_z = []
        self.target_x = []
        self.target_y = []
        self.target_z = []
        self.velocity_x = []
        self.velocity_y = []
        self.velocity_z = []

    def InitializeSolutionStep(self):
        super(DEMAnalysisStage2DSpRigidFem, self).InitializeSolutionStep()

        self.multiaxial_control_module.ExecuteInitializeSolutionStep()

    def ComputePorosityParameters(self, spheres_model_part):

        total_cylinders_volume = 0.0 # Remember we are here in 2D with 1m depth
        self.unit_depth = 1.0
        for node in spheres_model_part.Nodes:
            radius = node.GetSolutionStepValue(KratosMultiphysics.RADIUS)
            total_cylinders_volume +=  math.pi * radius * radius * self.unit_depth

        total_spheres_volume = (2.0/3.0) * total_cylinders_volume
        actual_solid_fraction = total_spheres_volume / self.ComputeSpecimenFullVolume()
        actual_porosity = 1.0 - actual_solid_fraction

        sandstone_target_porosity = 0.25 # Porosity values in sandstones are in the range of 25%
        target_solid_fraction = 1.0 - sandstone_target_porosity

        sp_porosity_multiplier = target_solid_fraction / actual_solid_fraction

        return sandstone_target_porosity, actual_porosity, sp_porosity_multiplier

    def ComputeSpecimenFullVolume(self):

        if self.test_number < 5:
            return math.pi * (self.outer_radius * self.outer_radius - self.inner_radius * self.inner_radius) * self.unit_depth
        else:
            return (4.0 * self.outer_radius * self.outer_radius - math.pi * self.inner_radius * self.inner_radius) * self.unit_depth

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

    def SetSkinManually(self):

        if not self.automatic_skin_computation and not self.respect_preprocessor_marked_skin:
            self.PreUtilities.ResetSkinParticles(self.spheres_model_part)
            self.PreUtilities.SetSkinParticlesInnerCircularBoundary(self.spheres_model_part, self.inner_radius, self.inner_skin_factor * 0.5 * self.inner_mesh_diameter)
            if self.test_number < 5: # CTWs
                self.PreUtilities.SetSkinParticlesOuterCircularBoundary(self.spheres_model_part, self.outer_radius, self.outer_skin_factor * 0.5 * self.outer_mesh_diameter)
            else: # Blind
                self.PreUtilities.SetSkinParticlesOuterSquaredBoundary(self.spheres_model_part, self.outer_radius, self.center, self.outer_skin_factor * 0.5 * self.outer_mesh_diameter)

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

    def FinalizeSolutionStep(self):
        super(DEMAnalysisStage2DSpRigidFem, self).FinalizeSolutionStep()
        self.multiaxial_control_module.ExecuteFinalizeSolutionStep()

        self.times.append(self.time)
        # self.reaction_x.append(self.rigid_face_model_part.Nodes[11].GetSolutionStepValue(KratosMultiphysics.DEMApplication.REACTION_STRESS_X))
        # self.reaction_y.append(self.rigid_face_model_part.Nodes[12].GetSolutionStepValue(KratosMultiphysics.DEMApplication.REACTION_STRESS_Y))
        # self.reaction_z.append(self.spheres_model_part.Nodes[9].GetSolutionStepValue(KratosMultiphysics.DEMApplication.REACTION_STRESS_Z))
        # self.target_x.append(self.rigid_face_model_part.Nodes[11].GetSolutionStepValue(KratosMultiphysics.DEMApplication.TARGET_STRESS_X))
        # self.target_y.append(self.rigid_face_model_part.Nodes[12].GetSolutionStepValue(KratosMultiphysics.DEMApplication.TARGET_STRESS_Y))
        # self.target_z.append(self.spheres_model_part.Nodes[9].GetSolutionStepValue(KratosMultiphysics.DEMApplication.TARGET_STRESS_Z))
        # self.velocity_x.append(self.rigid_face_model_part.Nodes[11].GetSolutionStepValue(KratosMultiphysics.DEMApplication.LOADING_VELOCITY_X))
        # self.velocity_y.append(self.rigid_face_model_part.Nodes[12].GetSolutionStepValue(KratosMultiphysics.DEMApplication.LOADING_VELOCITY_Y))
        # self.velocity_z.append(self.spheres_model_part.Nodes[9].GetSolutionStepValue(KratosMultiphysics.DEMApplication.LOADING_VELOCITY_Z))
        self.reaction_x.append(self.rigid_face_model_part.Nodes[3525].GetSolutionStepValue(KratosMultiphysics.DEMApplication.REACTION_STRESS_X))
        self.reaction_y.append(self.rigid_face_model_part.Nodes[3514].GetSolutionStepValue(KratosMultiphysics.DEMApplication.REACTION_STRESS_Y))
        self.reaction_z.append(self.spheres_model_part.Nodes[1].GetSolutionStepValue(KratosMultiphysics.DEMApplication.REACTION_STRESS_Z))
        self.smoothed_reaction_x.append(self.rigid_face_model_part.Nodes[3525].GetSolutionStepValue(KratosMultiphysics.DEMApplication.SMOOTHED_REACTION_STRESS_X))
        self.smoothed_reaction_y.append(self.rigid_face_model_part.Nodes[3514].GetSolutionStepValue(KratosMultiphysics.DEMApplication.SMOOTHED_REACTION_STRESS_Y))
        self.smoothed_reaction_z.append(self.spheres_model_part.Nodes[1].GetSolutionStepValue(KratosMultiphysics.DEMApplication.SMOOTHED_REACTION_STRESS_Z))
        self.elastic_reaction_x.append(self.rigid_face_model_part.Nodes[3525].GetSolutionStepValue(KratosMultiphysics.DEMApplication.ELASTIC_REACTION_STRESS_X))
        self.elastic_reaction_y.append(self.rigid_face_model_part.Nodes[3514].GetSolutionStepValue(KratosMultiphysics.DEMApplication.ELASTIC_REACTION_STRESS_Y))
        self.elastic_reaction_z.append(self.spheres_model_part.Nodes[1].GetSolutionStepValue(KratosMultiphysics.DEMApplication.ELASTIC_REACTION_STRESS_Z))
        self.smoothed_elastic_reaction_x.append(self.rigid_face_model_part.Nodes[3525].GetSolutionStepValue(KratosMultiphysics.DEMApplication.SMOOTHED_ELASTIC_REACTION_STRESS_X))
        self.smoothed_elastic_reaction_y.append(self.rigid_face_model_part.Nodes[3514].GetSolutionStepValue(KratosMultiphysics.DEMApplication.SMOOTHED_ELASTIC_REACTION_STRESS_Y))
        self.smoothed_elastic_reaction_z.append(self.spheres_model_part.Nodes[1].GetSolutionStepValue(KratosMultiphysics.DEMApplication.SMOOTHED_ELASTIC_REACTION_STRESS_Z))
        self.target_x.append(self.rigid_face_model_part.Nodes[3525].GetSolutionStepValue(KratosMultiphysics.DEMApplication.TARGET_STRESS_X))
        self.target_y.append(self.rigid_face_model_part.Nodes[3514].GetSolutionStepValue(KratosMultiphysics.DEMApplication.TARGET_STRESS_Y))
        self.target_z.append(self.spheres_model_part.Nodes[1].GetSolutionStepValue(KratosMultiphysics.DEMApplication.TARGET_STRESS_Z))
        self.velocity_x.append(self.rigid_face_model_part.Nodes[3525].GetSolutionStepValue(KratosMultiphysics.DEMApplication.LOADING_VELOCITY_X))
        self.velocity_y.append(self.rigid_face_model_part.Nodes[3514].GetSolutionStepValue(KratosMultiphysics.DEMApplication.LOADING_VELOCITY_Y))
        self.velocity_z.append(self.spheres_model_part.Nodes[1].GetSolutionStepValue(KratosMultiphysics.DEMApplication.LOADING_VELOCITY_Z))

    def PrintResultsForGid(self, time):
        super(DEMAnalysisStage2DSpRigidFem, self).PrintResultsForGid(time)

        ## TODO
        import matplotlib.pyplot as plt 

        f = plt.figure()

        plt.plot(self.times, self.reaction_x, color='b', label='reaction_x')
        plt.plot(self.times, self.reaction_y, color='r', label='reaction_y')
        plt.plot(self.times, self.reaction_z, color='g', label='reaction_z')
        plt.plot(self.times, self.smoothed_reaction_x, linewidth='2', color='b', label='smoothed_reaction_x')
        plt.plot(self.times, self.smoothed_reaction_y, linewidth='2', color='r', label='smoothed_reaction_y')
        plt.plot(self.times, self.smoothed_reaction_z, linewidth='2', color='g', label='smoothed_reaction_z')
        plt.plot(self.times, self.elastic_reaction_x, linestyle='dotted', color='b', label='elastic_reaction_x')
        plt.plot(self.times, self.elastic_reaction_y, linestyle='dotted', color='r', label='elastic_reaction_y')
        plt.plot(self.times, self.elastic_reaction_z, linestyle='dotted', color='g', label='elastic_reaction_z')
        plt.plot(self.times, self.smoothed_elastic_reaction_x, linestyle='dotted', linewidth='2', color='b', label='smoothed_elastic_reaction_x')
        plt.plot(self.times, self.smoothed_elastic_reaction_y, linestyle='dotted', linewidth='2', color='r', label='smoothed_elastic_reaction_y')
        plt.plot(self.times, self.smoothed_elastic_reaction_z, linestyle='dotted', linewidth='2', color='g', label='smoothed_elastic_reaction_z')
        plt.plot(self.times, self.target_x, '--', linewidth='3', color='b', label='target_x')
        plt.plot(self.times, self.target_y, '--', linewidth='3', color='r', label='target_y')
        plt.plot(self.times, self.target_z, '--', linewidth='3', color='g', label='target_z')

        plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
        # plt.legend()
        max_stress = max(abs(min(self.target_x)),max(self.target_x),
                         abs(min(self.target_y)),max(self.target_y),
                         abs(min(self.target_z)),max(self.target_z))
        plt.ylim(-2*max_stress,2*max_stress)

        # naming the x axis 
        plt.xlabel('Time (s)') 
        # naming the y axis 
        plt.ylabel('Stress (Pa)') 
        # giving a title to my graph 
        plt.title('Reaction vs target stresses') 

        f.savefig("forces.pdf", bbox_inches='tight')
        plt.close()

        f = plt.figure()

        plt.plot(self.times, self.velocity_x, label='velocity_x')
        plt.plot(self.times, self.velocity_y, label='velocity_y')
        plt.plot(self.times, self.velocity_z, label='velocity_z')

        plt.legend()

        # naming the x axis 
        plt.xlabel('Time (s)') 
        # naming the y axis 
        plt.ylabel('Velocity (m/s)') 
        # giving a title to my graph 
        plt.title('Loading velocity') 

        f.savefig("velocities.pdf", bbox_inches='tight')
        plt.close()
        ## TODO

        # DemFem.DemStructuresCouplingUtilities().MarkBrokenSpheres(self.ring_submodelpart)
        # self.creator_destructor.MarkParticlesForErasingGivenCylinder(self.ring_submodelpart, self.center, self.axis, self.radius_to_delete_sp)

        # DemFem.DemStructuresCouplingUtilities().ComputeSandProductionWithDepthFirstSearchNonRecursiveImplementation(self.ring_submodelpart, self.rigid_face_model_part, self.time)
        # DemFem.DemStructuresCouplingUtilities().ComputeSandProduction(self.ring_submodelpart, self.rigid_face_model_part, self.time)

    def AdditionalFinalizeOperations(self):

        spheres_mp_filename_post = self.problem_name + 'DEM_Post'

        if self.write_mdpa_from_results:
            mesh_creator_sphere.WriteSphereMdpaFromResults(self.problem_name + 'DEM', self.main_path, spheres_mp_filename_post, self.spheres_model_part)

if __name__ == "__main__":
    DEMAnalysisStage2DSpRigidFem(model, project_parameters).Run()
