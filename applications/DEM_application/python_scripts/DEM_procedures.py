from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *

import math
import DEM_material_test_script
import os
import shutil
import sys
from glob import glob

def Var_Translator(variable):

    if (variable == "OFF" or variable == "0" or variable == 0):
        variable = 0
    else:
        variable = 1

    return variable


class MdpaCreator(object):

    def __init__(self, path, DEM_parameters):

        self.DEM_parameters = DEM_parameters
        self.current_path = path

        # Creating necessary directories

        self.post_mdpas = os.path.join(str(self.current_path), str(self.DEM_parameters.problem_name) + '_post_mdpas')
        os.chdir(self.current_path)
        if not os.path.isdir(self.post_mdpas):
            os.makedirs(str(self.post_mdpas))

    def WriteMdpa(self, model_part):
        os.chdir(self.post_mdpas)
        time = model_part.ProcessInfo.GetValue(TIME)
        mdpa = open(str(self.DEM_parameters.problem_name) + '_post_' + str(time) + '.mdpa', 'w'+'\n')
        mdpa.write('Begin ModelPartData'+'\n')
        mdpa.write('//  VARIABLE_NAME value')
        mdpa.write('End ModelPartData'+'\n'+'\n'+'\n'+'\n')
        mdpa.write('Begin Nodes'+'\n')

        for node in model_part.Nodes:
            mdpa.write(str(node.Id) + ' ' + str(node.X) + ' ' + str(node.Y) + ' ' + str(node.Z)+'\n')
        mdpa.write('End Nodes'+'\n'+'\n')

        mdpa.write('Begin Elements SphericParticle3D'+'\n')
        for element in model_part.Elements:
            mdpa.write(str(element.Id) + ' ' +'1'+' ' + str(element.GetNode(0).Id )+'\n')
        mdpa.write('End Elements'+'\n'+'\n')

        fixed = 0; #how to read fixed? it can be either a flag or a nodal variable property
        self.WriteVariableData(RADIUS, mdpa, model_part)
        #self.WriteVariableData(VELOCITY_X, mdpa, model_part)
        #self.WriteVariableData(VELOCITY_Y, mdpa, model_part)
        #self.WriteVariableData(VELOCITY_Z, mdpa, model_part)


    def WriteVariableData(self, variable_name, mdpa, model_part):

        mdpa.write('Begin NodalData '+str(variable_name)+'\n')
        for node in model_part.Nodes:
            mdpa.write(str(node.Id) + ' ' + str(0) + ' ' + str(node.GetSolutionStepValue(variable_name))+'\n')
        mdpa.write('End NodalData'+'\n'+'\n')


class SetOfModelParts(object):
    def __init__(self, spheres_model_part, rigid_face_model_part, cluster_model_part, DEM_inlet_model_part, mapping_model_part, contact_model_part):
        self.spheres_model_part    = spheres_model_part
        self.rigid_face_model_part = rigid_face_model_part
        self.cluster_model_part    = cluster_model_part
        self.DEM_inlet_model_part  = DEM_inlet_model_part
        self.mapping_model_part    = mapping_model_part
        self.contact_model_part    = contact_model_part

class GranulometryUtils(object):

    def __init__(self, domain_volume, model_part):

        if (domain_volume <= 0.0):
            raise ValueError("Error: The input domain volume must be strictly positive!")

        self.spheres_model_part = model_part
        self.UpdateData(domain_volume)

    def UpdateData(self, domain_volume):

        self.physics_calculator = SphericElementGlobalPhysicsCalculator(self.spheres_model_part)
        self.number_of_spheres  = self.spheres_model_part.NumberOfElements(0)
        self.solid_volume       = self.physics_calculator.CalculateTotalVolume(self.spheres_model_part)
        self.d_50               = self.physics_calculator.CalculateD50(self.spheres_model_part)

        if (self.number_of_spheres == 0):
            self.spheres_per_area = 0.0
        else:
            self.spheres_per_area = domain_volume / self.number_of_spheres

        self.voids_volume    = domain_volume - self.solid_volume
        self.global_porosity = self.voids_volume / domain_volume

    def PrintCurrentData(self):

        print("number_of_spheres: ", self.number_of_spheres)
        print("solid volume: ", self.solid_volume)
        print("voids volume: ", self.voids_volume)
        print("global porosity: ", self.global_porosity)
        print("D50: ", self.d_50)
        print("spheres per area unit: ", self.spheres_per_area)


class PostUtils(object):

    def __init__(self, DEM_parameters, spheres_model_part):

        self.DEM_parameters = DEM_parameters
        self.spheres_model_part = spheres_model_part
        self.post_utilities = PostUtilities()

        self.vel_trap_graph_counter = 0
        self.vel_trap_graph_frequency = int(DEM_parameters.VelTrapGraphExportFreq/spheres_model_part.ProcessInfo.GetValue(DELTA_TIME))
        if self.vel_trap_graph_frequency < 1:
            self.vel_trap_graph_frequency = 1 #that means it is not possible to print results with a higher frequency than the computations delta time

        self.previous_vector_of_inner_nodes = []
        self.previous_time = 0.0;
        
    def Flush(self,a):
        a.flush()

    def ComputeMeanVelocitiesinTrap(self, file_name, time_dem):

        if (self.DEM_parameters.VelocityTrapOption):
            compute_flow = False

            self.vel_trap_graph_counter += 1

            if (self.vel_trap_graph_counter == self.vel_trap_graph_frequency):
                self.vel_trap_graph_counter = 0
                average_velocity = Array3()
                low_point = Array3()

                low_point[0] = self.DEM_parameters.VelocityTrapMinX
                low_point[1] = self.DEM_parameters.VelocityTrapMinY
                low_point[2] = self.DEM_parameters.VelocityTrapMinZ
                high_point = Array3()
                high_point[0] = self.DEM_parameters.VelocityTrapMaxX
                high_point[1] = self.DEM_parameters.VelocityTrapMaxY
                high_point[2] = self.DEM_parameters.VelocityTrapMaxZ

                average_velocity = self.post_utilities.VelocityTrap(self.spheres_model_part, low_point, high_point)

                if compute_flow == True:
                    vector_of_inner_nodes = []
                    for node in self.spheres_model_part.Nodes:
                        if (node.X > low_point[0]) & (node.Y > low_point[1]) & (node.Z > low_point[2]) & (node.X < high_point[0]) & (node.Y < high_point[1]) & (node.Z < high_point[2]) :
                            vector_of_inner_nodes.append(node)

                    crossing_spheres = 0
                    crossing_volume = 0.0

                    for node in vector_of_inner_nodes:
                        id_found = False
                        for previous_node in self.previous_vector_of_inner_nodes:
                            if node.Id == previous_node.Id:
                                id_found = True
                                break
                        if id_found == False: #This only happens if None of the previous nodes were capable of setting id_found = True.
                            crossing_spheres = crossing_spheres + 1
                            radius = node.GetSolutionStepValue(RADIUS)
                            crossing_volume = crossing_volume + 4.0/3.0 * math.pi * radius*radius*radius

                    time_between_measures = self.spheres_model_part.ProcessInfo.GetValue(TIME) - self.previous_time
                    number_of_spheres_flow = float(crossing_spheres) / time_between_measures
                    net_volume_flow = crossing_volume / time_between_measures

                    self.previous_time = self.spheres_model_part.ProcessInfo.GetValue(TIME)
                    self.previous_vector_of_inner_nodes = vector_of_inner_nodes


                f = open(file_name, 'a')
                tmp = str(time_dem) + "   " + str(average_velocity[0]) + "   " + str(average_velocity[1]) + "   " + str(average_velocity[2])
                if compute_flow == True:
                    tmp = tmp + "   " + str(net_volume_flow)  + "   " + str(number_of_spheres_flow)
                tmp = tmp + "\n"

                f.write(tmp)
                self.Flush(f)

    def PrintEulerAngles(self, spheres_model_part, cluster_model_part):
        PostUtilities().ComputeEulerAngles(spheres_model_part, cluster_model_part)


class DEMEnergyCalculator(object):
    
    def __init__(self, DEM_parameters, spheres_model_part, cluster_model_part, energy_plot):
        
        self.calculate_option = False
        
        if (hasattr(DEM_parameters, "EnergyCalculationOption")):
            if (DEM_parameters.EnergyCalculationOption):
                self.calculate_option = True
                self.DEM_parameters = DEM_parameters
                self.SpheresModelPart = spheres_model_part;
                self.ClusterModelPart = cluster_model_part;                
                self.energy_plot = open(energy_plot, 'w');
                self.SpheresEnergyUtil = SphericElementGlobalPhysicsCalculator(spheres_model_part);
                self.ClusterEnergyUtil = SphericElementGlobalPhysicsCalculator(cluster_model_part);
                self.PotentialEnergyReferencePoint          = Array3();
                self.PotentialEnergyReferencePoint[0]       = self.DEM_parameters.PotentialEnergyReferencePointX;
                self.PotentialEnergyReferencePoint[1]       = self.DEM_parameters.PotentialEnergyReferencePointY;
                self.PotentialEnergyReferencePoint[2]       = self.DEM_parameters.PotentialEnergyReferencePointZ;
                self.translational_kinematic_energy         = 0.0;
                self.rotational_kinematic_energy            = 0.0;
                self.kinematic_energy                       = 0.0;
                self.gravitational_energy                   = 0.0;
                self.elastic_energy                         = 0.0;
                self.inelastic_frictonal_energy             = 0.0;
                self.inelastic_viscodamping_energy          = 0.0;
                self.external_energy                        = 0.0;
                self.total_energy                           = 0.0;
                self.graph_frequency                        = int(self.DEM_parameters.GraphExportFreq/spheres_model_part.ProcessInfo.GetValue(DELTA_TIME));
                self.energy_graph_counter                   = 0;
                self.energy_plot.write(str("Time").rjust(9)+"   "+str("Trans kinematic energy").rjust(22)+"   "+str("Rot kinematic energy").rjust(20)+"   "+str("Kinematic energy").rjust(16)+"   "+str("Gravitational energy").rjust(20)+"   "+str("Elastic energy").rjust(14)+"   "+str("Frictonal energy").rjust(16)+"   "+str("Viscodamping energy").rjust(19)+"   "+str("Total energy").rjust(12)+"\n")

    def CalculateEnergyAndPlot(self, time):
        if self.calculate_option:
            if not hasattr(self.DEM_parameters, "TestType"):
                if (self.energy_graph_counter == self.graph_frequency):
                    self.energy_graph_counter = 0

                    self.CalculateEnergy();
                    self.PlotEnergyGraph(time);

                self.energy_graph_counter += 1

    def CalculateEnergy(self):

        self.translational_kinematic_energy = self.SpheresEnergyUtil.CalculateTranslationalKinematicEnergy(self.SpheresModelPart) + self.ClusterEnergyUtil.CalculateTranslationalKinematicEnergy(self.ClusterModelPart);
        self.rotational_kinematic_energy    = self.SpheresEnergyUtil.CalculateRotationalKinematicEnergy(self.SpheresModelPart) + self.ClusterEnergyUtil.CalculateRotationalKinematicEnergy(self.ClusterModelPart);
        self.kinematic_energy               = self.translational_kinematic_energy + self.rotational_kinematic_energy;
        self.gravitational_energy           = self.SpheresEnergyUtil.CalculateGravitationalPotentialEnergy(self.SpheresModelPart,self.PotentialEnergyReferencePoint) + self.ClusterEnergyUtil.CalculateGravitationalPotentialEnergy(self.ClusterModelPart,self.PotentialEnergyReferencePoint);
        self.elastic_energy                 = self.SpheresEnergyUtil.CalculateElasticEnergy(self.SpheresModelPart) + self.ClusterEnergyUtil.CalculateElasticEnergy(self.ClusterModelPart);
        self.inelastic_frictional_energy    = self.SpheresEnergyUtil.CalculateInelasticFrictionalEnergy(self.SpheresModelPart) + self.ClusterEnergyUtil.CalculateInelasticFrictionalEnergy(self.ClusterModelPart);
        self.inelastic_viscodamping_energy  = self.SpheresEnergyUtil.CalculateInelasticViscodampingEnergy(self.SpheresModelPart) + self.ClusterEnergyUtil.CalculateInelasticViscodampingEnergy(self.ClusterModelPart);
        self.total_energy                   = self.kinematic_energy + self.gravitational_energy + self.elastic_energy + self.inelastic_frictional_energy + self.inelastic_viscodamping_energy;
        
    def PlotEnergyGraph(self,time):

        plot_kinematic               = self.kinematic_energy
        plot_translational_kinematic = self.translational_kinematic_energy
        plot_rotational_kinematic    = self.rotational_kinematic_energy
        plot_gravitational           = self.gravitational_energy
        plot_elastic                 = self.elastic_energy
        plot_inelastic_frictional    = self.inelastic_frictional_energy
        plot_inelastic_viscodamping  = self.inelastic_viscodamping_energy
        plot_total                   = self.total_energy
        self.energy_plot.write( str("%.8g"%time).rjust(9)+"   "+str("%.6g"%plot_translational_kinematic).rjust(22)+"   "+str("%.6g"%plot_rotational_kinematic).rjust(20)+"   "+str("%.6g"%plot_kinematic).rjust(16)+"   "+str("%.6g"%plot_gravitational).rjust(20)+"   "+str("%.6g"%plot_elastic).rjust(14)+"   "+str("%.6g"%plot_inelastic_frictional).rjust(16)+"   "+str("%.6g"%plot_inelastic_viscodamping).rjust(19)+"   "+str("%.6g"%plot_total).rjust(12)+'\n' )
        self.energy_plot.flush()
        
    def FinalizeEnergyPlot(self):
        if self.calculate_option:
            self.energy_plot.close


class Procedures(object):

    def __init__(self, DEM_parameters):

        # GLOBAL VARIABLES OF THE SCRIPT
        # Defining list of skin particles (For a test tube of height 30 cm and diameter 15 cm)

        # Initialization of member variables
        self.DEM_parameters = DEM_parameters

        # SIMULATION FLAGS
        self.rotation_OPTION               = Var_Translator(self.DEM_parameters.RotationOption)
        self.bounding_box_OPTION           = Var_Translator(self.DEM_parameters.BoundingBoxOption)
        self.automatic_bounding_box_OPTION = Var_Translator(self.DEM_parameters.AutomaticBoundingBoxOption)
        self.contact_mesh_OPTION           = Var_Translator(self.DEM_parameters.ContactMeshOption)

        # SIMULATION SETTINGS
        self.b_box_minX = self.DEM_parameters.BoundingBoxMinX
        self.b_box_minY = self.DEM_parameters.BoundingBoxMinY
        self.b_box_minZ = self.DEM_parameters.BoundingBoxMinZ
        self.b_box_maxX = self.DEM_parameters.BoundingBoxMaxX
        self.b_box_maxY = self.DEM_parameters.BoundingBoxMaxY
        self.b_box_maxZ = self.DEM_parameters.BoundingBoxMaxZ
        self.bounding_box_enlargement_factor = self.DEM_parameters.BoundingBoxEnlargementFactor

        # MODEL
        self.domain_size = self.DEM_parameters.Dimension
        
    def SetScheme(self):
        if (self.DEM_parameters.IntegrationScheme == 'Forward_Euler'):
            scheme = ForwardEulerScheme()
        elif (self.DEM_parameters.IntegrationScheme == 'Symplectic_Euler'):
            scheme = SymplecticEulerScheme()
        elif (self.DEM_parameters.IntegrationScheme == 'Taylor_Scheme'):
            scheme = TaylorScheme()
        elif (self.DEM_parameters.IntegrationScheme == 'Newmark_Beta_Method'):
            scheme = NewmarkBetaScheme(0.5, 0.25)
        elif (self.DEM_parameters.IntegrationScheme == 'Verlet_Velocity'):
            scheme = VerletVelocityScheme()
        else:
            self.KRATOSprint('Error: selected scheme not defined. Please select a different scheme')
            sys.exit("\nExecution was aborted.\n")
        return scheme
        
    def AddAllVariablesInAllModelParts(self, solver, scheme, spheres_model_part, cluster_model_part, DEM_inlet_model_part, rigid_face_model_part, DEM_parameters):
        self.solver=solver
        self.scheme=scheme
        self.AddCommonVariables(spheres_model_part, DEM_parameters)
        self.AddSpheresVariables(spheres_model_part, DEM_parameters)
        self.AddMpiVariables(spheres_model_part)
        self.solver.AddAdditionalVariables(spheres_model_part, DEM_parameters)
        #self.scheme.AddSpheresVariables(spheres_model_part)
        self.AddCommonVariables(cluster_model_part, DEM_parameters)
        self.AddClusterVariables(cluster_model_part, DEM_parameters)
        self.AddMpiVariables(cluster_model_part)
        #self.scheme.AddClustersVariables(cluster_model_part)
        self.AddCommonVariables(DEM_inlet_model_part, DEM_parameters)
        self.AddSpheresVariables(DEM_inlet_model_part, DEM_parameters)
        self.solver.AddAdditionalVariables(DEM_inlet_model_part, DEM_parameters)
        #self.scheme.AddSpheresVariables(DEM_inlet_model_part)
        self.AddCommonVariables(rigid_face_model_part, DEM_parameters)
        self.AddRigidFaceVariables(rigid_face_model_part, DEM_parameters)
        self.AddMpiVariables(rigid_face_model_part)

    def AddCommonVariables(self, model_part, DEM_parameters):
        model_part.AddNodalSolutionStepVariable(VELOCITY)
        model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(DELTA_DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(TOTAL_FORCES)

    def AddSpheresVariables(self, model_part, DEM_parameters):

        # KINEMATIC
        model_part.AddNodalSolutionStepVariable(DELTA_ROTATION) #TODO: only if self.DEM_parameters.RotationOption! Check that no one accesses them in c++ without checking the rotation option
        model_part.AddNodalSolutionStepVariable(PARTICLE_ROTATION_ANGLE)  #TODO: only if self.DEM_parameters.RotationOption! Check that no one accesses them in c++ without checking the rotation option
        model_part.AddNodalSolutionStepVariable(ANGULAR_VELOCITY)  #TODO: only if self.DEM_parameters.RotationOption! Check that no one accesses them in c++ without checking the rotation option

        # FORCES
        model_part.AddNodalSolutionStepVariable(ELASTIC_FORCES)
        model_part.AddNodalSolutionStepVariable(LOCAL_CONTACT_FORCE)
        model_part.AddNodalSolutionStepVariable(CONTACT_FORCES)
        model_part.AddNodalSolutionStepVariable(RIGID_ELEMENT_FORCE)
        model_part.AddNodalSolutionStepVariable(DAMP_FORCES)
        model_part.AddNodalSolutionStepVariable(PARTICLE_MOMENT) #TODO: only if self.DEM_parameters.RotationOption! Check that no one accesses them in c++ without checking the rotation option
        model_part.AddNodalSolutionStepVariable(EXTERNAL_APPLIED_FORCE)
        model_part.AddNodalSolutionStepVariable(EXTERNAL_APPLIED_MOMENT)

        # BASIC PARTICLE PROPERTIES
        model_part.AddNodalSolutionStepVariable(RADIUS)
        model_part.AddNodalSolutionStepVariable(NODAL_MASS)
        model_part.AddNodalSolutionStepVariable(REPRESENTATIVE_VOLUME)
        model_part.AddNodalSolutionStepVariable(NEIGHBOUR_SIZE)

        # ROTATION RELATED PROPERTIES
        if (Var_Translator(self.DEM_parameters.RotationOption)):
            model_part.AddNodalSolutionStepVariable(PARTICLE_MOMENT_OF_INERTIA) #TODO: only if self.DEM_parameters.RotationOption! Check that no one accesses them in c++ without checking the rotation option
            model_part.AddNodalSolutionStepVariable(PARTICLE_ROTATION_DAMP_RATIO) #TODO: only if self.DEM_parameters.RotationOption! Check that no one accesses them in c++ without checking the rotation option
            if (Var_Translator(self.DEM_parameters.RollingFrictionOption)):
                model_part.AddNodalSolutionStepVariable(ROLLING_FRICTION)
                model_part.AddNodalSolutionStepVariable(ROLLING_RESISTANCE_MOMENT)

        # OTHER PROPERTIES
        model_part.AddNodalSolutionStepVariable(PARTICLE_MATERIAL)   # Colour defined in GiD
        
        if (hasattr(self.DEM_parameters, "PostSkinSphere")):
            if (Var_Translator(self.DEM_parameters.PostSkinSphere)):
                model_part.AddNodalSolutionStepVariable(SKIN_SPHERE)
        
        # LOCAL AXIS
        if (Var_Translator(self.DEM_parameters.PostEulerAngles)):
            model_part.AddNodalSolutionStepVariable(EULER_ANGLES)


        if ((hasattr(self.DEM_parameters, "PostStressStrainOption")) and self.DEM_parameters.PostStressStrainOption):
            model_part.AddNodalSolutionStepVariable(DEM_STRESS_TENSOR)
        if (self.solver.compute_stress_tensor_option):
            model_part.AddNodalSolutionStepVariable(FORCE_REACTION)
            model_part.AddNodalSolutionStepVariable(MOMENT_REACTION)

        if (self.solver.poisson_ratio_option):
            model_part.AddNodalSolutionStepVariable(POISSON_VALUE)

        # Nano Particle
        if self.DEM_parameters.ElementType == "SwimmingNanoParticle":
            model_part.AddNodalSolutionStepVariable(CATION_CONCENTRATION)
            model_part.AddNodalSolutionStepVariable(DRAG_COEFFICIENT)
        # ONLY VISUALIZATION
        if (Var_Translator(self.DEM_parameters.PostExportId)):
            model_part.AddNodalSolutionStepVariable(EXPORT_ID)

        #model_part.AddNodalSolutionStepVariable(SPRAYED_MATERIAL)

    def AddRigidFaceVariables(self, model_part, DEM_parameters):

        model_part.AddNodalSolutionStepVariable(ELASTIC_FORCES)
        model_part.AddNodalSolutionStepVariable(CONTACT_FORCES)
        model_part.AddNodalSolutionStepVariable(DEM_PRESSURE)
        model_part.AddNodalSolutionStepVariable(TANGENTIAL_ELASTIC_FORCES)
        model_part.AddNodalSolutionStepVariable(SHEAR_STRESS)
        model_part.AddNodalSolutionStepVariable(DEM_NODAL_AREA)
        model_part.AddNodalSolutionStepVariable(NON_DIMENSIONAL_VOLUME_WEAR)
        model_part.AddNodalSolutionStepVariable(IMPACT_WEAR)

    def AddElasticFaceVariables(self, model_part, DEM_parameters): #Only used in CSM coupling
        self.AddRigidFaceVariables(model_part,self.DEM_parameters)
        model_part.AddNodalSolutionStepVariable(TOTAL_FORCES)

    def AddClusterVariables(self, model_part, DEM_parameters):
        # KINEMATIC
        model_part.AddNodalSolutionStepVariable(DELTA_DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(PARTICLE_ROTATION_ANGLE)
        model_part.AddNodalSolutionStepVariable(DELTA_ROTATION)
        model_part.AddNodalSolutionStepVariable(ANGULAR_VELOCITY)
        model_part.AddNodalSolutionStepVariable(LOCAL_ANGULAR_VELOCITY)
        model_part.AddNodalSolutionStepVariable(ORIENTATION_REAL)
        model_part.AddNodalSolutionStepVariable(ORIENTATION_IMAG)
        model_part.AddNodalSolutionStepVariable(ORIENTATION)
        model_part.AddNodalSolutionStepVariable(ANGULAR_MOMENTUM)

        # FORCES
        model_part.AddNodalSolutionStepVariable(TOTAL_FORCES)
        model_part.AddNodalSolutionStepVariable(RIGID_ELEMENT_FORCE)
        model_part.AddNodalSolutionStepVariable(PARTICLE_MOMENT)
        model_part.AddNodalSolutionStepVariable(EXTERNAL_APPLIED_FORCE)
        model_part.AddNodalSolutionStepVariable(EXTERNAL_APPLIED_MOMENT)
        
        # PHYSICAL PROPERTIES
        model_part.AddNodalSolutionStepVariable(PRINCIPAL_MOMENTS_OF_INERTIA)
        model_part.AddNodalSolutionStepVariable(CLUSTER_VOLUME)
        model_part.AddNodalSolutionStepVariable(NODAL_MASS)
        model_part.AddNodalSolutionStepVariable(CHARACTERISTIC_LENGTH)
        model_part.AddNodalSolutionStepVariable(PARTICLE_DENSITY)

        # OTHER PROPERTIES
        model_part.AddNodalSolutionStepVariable(PARTICLE_MATERIAL)   # Colour defined in GiD
        
        # LOCAL AXIS
        if (Var_Translator(self.DEM_parameters.PostEulerAngles)):
            model_part.AddNodalSolutionStepVariable(EULER_ANGLES)

    def AddMpiVariables(self, model_part):
        pass
    
    def SetUpBufferSizeInAllModelParts(self, spheres_model_part, spheres_b_size, cluster_model_part, clusters_b_size, DEM_inlet_model_part, inlet_b_size, rigid_face_model_part, rigid_b_size):
        spheres_model_part.SetBufferSize(spheres_b_size)
        cluster_model_part.SetBufferSize(clusters_b_size)
        DEM_inlet_model_part.SetBufferSize(inlet_b_size)
        rigid_face_model_part.SetBufferSize(rigid_b_size)

    def FindMaxNodeIdAccrossModelParts(self, creator_destructor, all_model_parts):
        
        max_node_Id = creator_destructor.FindMaxNodeIdInModelPart(all_model_parts.spheres_model_part)
        max_elem_Id = creator_destructor.FindMaxElementIdInModelPart(all_model_parts.spheres_model_part)
        max_FEM_node_Id = creator_destructor.FindMaxNodeIdInModelPart(all_model_parts.rigid_face_model_part)
        max_cluster_node_Id = creator_destructor.FindMaxNodeIdInModelPart(all_model_parts.cluster_model_part)
        max_Id = max(max_FEM_node_Id, max_node_Id, max_elem_Id, max_cluster_node_Id)
        
        return max_Id

    def ModelData(self, spheres_model_part, solver):
        
        contact_model_part = solver.contact_model_part
        # Previous Calculations.
        Model_Data = open('Model_Data.txt', 'w')

        # mean radius, and standard deviation:
        i = 0.0
        sum_radi = 0.0
        partial_sum_squared = 0.0
        total_sum_squared = 0.0
        volume = 0.0
        area = 0.0
        mean = 0.0
        var  = 0.0
        rel_std_dev = 0.0

        for node in spheres_model_part.Nodes:

            sum_radi += node.GetSolutionStepValue(RADIUS)
            partial_sum_squared = node.GetSolutionStepValue(RADIUS) ** 2.0
            total_sum_squared += partial_sum_squared
            volume += 4 * 3.141592 / 3 * node.GetSolutionStepValue(RADIUS) ** 3.0
            area += 3.141592 * partial_sum_squared
            i += 1.0

        if (i>0.0):
            mean = sum_radi / i
            var = total_sum_squared / i - mean ** 2.0
        std_dev = 0.0

        if (abs(var) > 1e-9):
            std_dev = var ** 0.5

        if (i>0.0):
            rel_std_dev = std_dev / mean

        Model_Data.write("Radius Mean: " + str(mean) + '\n')
        Model_Data.write("Std Deviation: " + str(std_dev) + '\n')
        Model_Data.write("Relative Std Deviation: " + str(rel_std_dev) + '\n')
        Model_Data.write("Total Particle Volume 3D: " + str(volume) + '\n')
        Model_Data.write("Total Particle Area 2D: " + str(area) + '\n')
        Model_Data.write('\n')

        Total_Particles = len(spheres_model_part.Nodes)

        Total_Contacts = 0

        if solver.continuum_type:
            Coordination_Number = 0.0

            if (self.contact_mesh_OPTION):
                for bar in contact_model_part.Elements:
                    Total_Contacts += 1
                if (Total_Particles):
                    Coordination_Number = 2.0 * Total_Contacts / Total_Particles

            Model_Data.write("Total Number of Particles: " + str(Total_Particles) + '\n')
            Model_Data.write("Total Number of Contacts: " + str(Total_Contacts) + '\n')
            Model_Data.write("Coordination Number NC: " + str(Coordination_Number) + '\n')
            Model_Data.write('\n')
            #Model_Data.write("Volume Elements: " + str(total_volume) + '\n')            
            self.KRATOSprint ("Coordination Number: " + str(Coordination_Number) + "\n")
            
        Model_Data.close()

    def MeasureBOT(self, solver):

        tol = 2.0
        y_mean = 0.0
        counter = 0.0

        for node in self.BOT:
            r = node.GetSolutionStepValue(RADIUS)
            y = node.Y
            y_mean += (y - r) * r
            counter += r

        return (y_mean, counter)

    def MeasureTOP(self, solver):

        tol = 2.0
        y_mean = 0.0
        counter = 0.0

        for node in self.TOP:
            r = node.GetSolutionStepValue(RADIUS)
            y = node.Y

            y_mean += (y + r) * r
            counter += r

        return (y_mean, counter)

    def MonitorPhysicalProperties(self, model_part, physics_calculator, properties_list):

        # This function returns a list of arrays (also lists)
        # Each array contains the values of the physical properties at the current time
        time = model_part.ProcessInfo.GetValue(TIME)
        present_prop = []

        if (len(properties_list) == 0):  # The first array in the list only contains the entries names
            names = []
            names.append("time")
            names.append("mass")
            names.append("gravitational_energy")
            names.append("kinematic_energy")
            #names.append("elastic_energy")
            names.append("momentum")
            names.append("angular_momentum")
            names.append("total_energy")

            properties_list.append(names)

        # Calculating current values
        mass = physics_calculator.CalculateTotalMass(model_part)
        center = physics_calculator.CalculateCenterOfMass(model_part)
        initial_center = physics_calculator.GetInitialCenterOfMass()
        gravity_energy = physics_calculator.CalculateGravitationalPotentialEnergy(model_part, initial_center)
        kinematic_energy = physics_calculator.CalculateKinematicEnergy(model_part)
        #elastic_energy = physics_calculator.CalculateElasticEnergy(model_part)
        momentum = physics_calculator.CalculateTotalMomentum(model_part)
        angular_momentum = physics_calculator.CalulateTotalAngularMomentum(model_part)
        total_energy = gravity_energy + kinematic_energy #+ elastic_energy

        # Filling in the entries values corresponding to the entries names above
        present_prop.append(time)
        present_prop.append(mass)
        present_prop.append(gravity_energy)
        present_prop.append(kinematic_energy)
        #present_prop.append(elastic_energy)
        present_prop.append(momentum)
        present_prop.append(angular_momentum)
        present_prop.append(total_energy)

        properties_list.append(present_prop)

        return properties_list

    # def PlotPhysicalProperties(self, properties_list, path):

    # This function creates one graph for each physical property.
    # properties_list[0][0] = 'time'
    # properties_list[0][j] = 'property_j'
    # properties_list[i][j] = value of property_j at time properties_list[i][0]

        # n_measures     = len(properties_list)
        # entries        = properties_list[0]
        # n_entries      = len(entries)
        # time_vect      = []
        # os.chdir(path)

        # for j in range(1, n_measures):
            # time_vect.append(properties_list[j][0])

        # for i in range(1, n_entries):
            # prop_vect_i = []

            # for j in range(1, n_measures):
                # prop_i_j = properties_list[j][i]

                # if (hasattr(prop_i_j, '__getitem__')): # Checking if it is an iterable object (a vector). If yes, take the modulus
                    # mod_prop_i_j = 0.0

                    # for k in range(len(prop_i_j)):
                        # mod_prop_i_j += prop_i_j[k] * prop_i_j[k]

                    # prop_i_j = sqrt(mod_prop_i_j) # Euclidean norm

                # prop_vect_i.append(prop_i_j)

            # plt.figure(i)
            # plot = plt.plot(time_vect, prop_vect_i)
            # plt.xlabel(entries[0])
            # plt.ylabel(entries[i])
            # plt.title('Evolution of ' + entries[i] + ' in time')
            # plt.savefig(entries[i] + '.pdf')

    def SetCustomSkin(self, spheres_model_part):

        for element in spheres_model_part.Elements:

            x = element.GetNode(0).X
            y = element.GetNode(0).Y
            #z = element.GetNode(0).Z

            if(x>21.1):
                element.GetNode(0).SetSolutionStepValue(SKIN_SPHERE,1)
            if(x<1.25):
                element.GetNode(0).SetSolutionStepValue(SKIN_SPHERE,1)
            if(y>1.9):
                element.GetNode(0).SetSolutionStepValue(SKIN_SPHERE,1)
            if(y<0.1):
                element.GetNode(0).SetSolutionStepValue(SKIN_SPHERE,1)

    def CreateDirectories(self, main_path, problem_name, run_code = ''):

        root             = os.path.join(main_path, problem_name)
        post_path        = root + '_Post_Files' + run_code
        data_and_results = root + '_Results_and_Data'
        graphs_path      = root + '_Graphs'
        MPI_results      = root + '_MPI_results'
        
        '''
        answer = input("\nWarning: If there already exists previous results, they are about to be deleted. Do you want to proceed (y/n)? ")
        if answer=='y':
            shutil.rmtree(os.path.join(main_path, problem_name + '_Post_Files'), ignore_errors = True)
            shutil.rmtree(os.path.join(main_path, problem_name + '_Graphs'    ), ignore_errors = True)
        else:
            sys.exit("\nExecution was aborted.\n")
        '''

        shutil.rmtree(os.path.join(main_path, problem_name + '_Post_Files'), ignore_errors = True)
        shutil.rmtree(os.path.join(main_path, problem_name + '_Graphs'    ), ignore_errors = True)

        for directory in [post_path, data_and_results, graphs_path, MPI_results]:
            if not os.path.isdir(directory):
                os.makedirs(str(directory))

        return [post_path, data_and_results, graphs_path, MPI_results]

    def FindMaxNodeIdInModelPart(self, model_part):

        maxid = 0

        for node in model_part.Nodes:
            if (node.Id > maxid):
                maxid = node.Id

        return maxid
    
    def SetBoundingBoxLimits(self, all_model_parts, creator_destructor):
        
        bounding_box_time_limits = []
        if (self.DEM_parameters.BoundingBoxOption == "ON"):
            self.SetBoundingBox(all_model_parts.spheres_model_part, all_model_parts.cluster_model_part, all_model_parts.rigid_face_model_part, creator_destructor)
            bounding_box_time_limits = [self.solver.bounding_box_start_time, self.solver.bounding_box_stop_time]
            return bounding_box_time_limits

    def SetBoundingBox(self, spheres_model_part, clusters_model_part, rigid_faces_model_part, creator_destructor):

        b_box_low = Array3()
        b_box_high = Array3()
        b_box_low[0] = self.b_box_minX
        b_box_low[1] = self.b_box_minY
        b_box_low[2] = self.b_box_minZ
        b_box_high[0] = self.b_box_maxX
        b_box_high[1] = self.b_box_maxY
        b_box_high[2] = self.b_box_maxZ

        creator_destructor.SetLowNode(b_box_low)
        creator_destructor.SetHighNode(b_box_high)
        creator_destructor.CalculateSurroundingBoundingBox(spheres_model_part, clusters_model_part, rigid_faces_model_part, self.bounding_box_enlargement_factor, self.automatic_bounding_box_OPTION)

    def DeleteFiles(self):
        files_to_delete_list = glob('*.time')
        for to_erase_file in files_to_delete_list:
            os.remove(to_erase_file)

    def PreProcessModel(self, DEM_parameters):
        pass
    
    def CheckVariableType(self, var, expected_type, msg):                  
        actual_type = type(var)
        if actual_type is int and expected_type is float:
            return
        if actual_type is not expected_type:
            self.KRATOSprint("**************************************************************************")
            self.KRATOSprint("ERROR: Input parameter of wrong type in file 'DEM_explicit_solver_var.py'." )
            a = str(expected_type)
            b = str(var)
            self.KRATOSprint("The type expected was "+ a + " but " + b +" was read.")                        
            self.KRATOSprint("**************************************************************************")                        
            sys.exit()
    
    def CheckInputParameters(self, DEM_parameters):        
        self.CheckVariableType(DEM_parameters.Dimension, int, "")
        self.CheckVariableType(DEM_parameters.BoundingBoxEnlargementFactor, float, "")
        self.CheckVariableType(DEM_parameters.BoundingBoxMaxX, float, "")
        self.CheckVariableType(DEM_parameters.BoundingBoxMaxY, float, "")
        self.CheckVariableType(DEM_parameters.BoundingBoxMaxZ, float, "")
        self.CheckVariableType(DEM_parameters.BoundingBoxMinX, float, "")
        self.CheckVariableType(DEM_parameters.BoundingBoxMinY, float, "")
        self.CheckVariableType(DEM_parameters.BoundingBoxMinZ, float, "")        
        self.CheckVariableType(DEM_parameters.GravityX, float, "")
        self.CheckVariableType(DEM_parameters.GravityY, float, "")
        self.CheckVariableType(DEM_parameters.GravityZ, float, "")    
        if hasattr(self.DEM_parameters, "AmplifiedSearchRadiusExtension"):
            self.CheckVariableType(DEM_parameters.AmplifiedSearchRadiusExtension, float, "")
        if hasattr(self.DEM_parameters, "MaxAmplificationRatioOfSearchRadius"):
            self.CheckVariableType(DEM_parameters.MaxAmplificationRatioOfSearchRadius, int, "")
        self.CheckVariableType(DEM_parameters.VirtualMassCoefficient, float, "")
        self.CheckVariableType(DEM_parameters.DeltaTimeSafetyFactor, float, "")
        self.CheckVariableType(DEM_parameters.MaxTimeStep, float, "")
        self.CheckVariableType(DEM_parameters.FinalTime, float, "")
        self.CheckVariableType(DEM_parameters.ControlTime, float, "")
        self.CheckVariableType(DEM_parameters.NeighbourSearchFrequency, int, "")
        self.CheckVariableType(DEM_parameters.GraphExportFreq, float, "")
        self.CheckVariableType(DEM_parameters.VelTrapGraphExportFreq, float, "")
        self.CheckVariableType(DEM_parameters.OutputTimeStep, float, "")

    def Flush(self,a):
        a.flush()
    
    def KRATOSprint(self,message):
        print(message)
        self.Flush(sys.stdout)


# #~CHARLIE~# Aixo no ho entenc
class DEMFEMProcedures(object):

    def __init__(self, DEM_parameters, graphs_path, spheres_model_part, RigidFace_model_part):

        # GLOBAL VARIABLES OF THE SCRIPT
        self.DEM_parameters = DEM_parameters

        if not hasattr(DEM_parameters, "TestType"):
            self.TestType = "None"
        # self.TestType = self.DEM_parameters.TestType

        # Initialization of member variables
        # SIMULATION FLAGS
        self.rotation_OPTION     = Var_Translator(self.DEM_parameters.RotationOption)
        self.bounding_box_OPTION = Var_Translator(self.DEM_parameters.BoundingBoxOption)
        self.contact_mesh_OPTION           = 0
        if (hasattr(DEM_parameters, "contact_mesh_OPTION")):
            self.contact_mesh_OPTION           = Var_Translator(self.DEM_parameters.ContactMeshOption)

        self.graphs_path = graphs_path
        self.spheres_model_part = spheres_model_part
        self.RigidFace_model_part = RigidFace_model_part
        #self.solver = solver
        self.aux = AuxiliaryUtilities()

        self.fem_mesh_nodes = []

        self.graph_counter = 0
        self.balls_graph_counter = 0

        self.graph_frequency        = int(self.DEM_parameters.GraphExportFreq/spheres_model_part.ProcessInfo.GetValue(DELTA_TIME))
        if self.graph_frequency < 1:
            self.graph_frequency = 1 #that means it is not possible to print results with a higher frequency than the computations delta time
        os.chdir(self.graphs_path)
        #self.graph_forces = open(self.DEM_parameters.problem_name +"_force_graph.grf", 'w')
        self.mesh_motion = DEMFEMUtilities()
        
        def Flush(self,a):
            a.flush()
                    
        def open_graph_files(self, RigidFace_model_part):
            #os.chdir(self.graphs_path)
            for mesh_number in range(0, self.RigidFace_model_part.NumberOfSubModelParts()):
                if (self.aux.GetIthSubModelPartData(self.RigidFace_model_part, mesh_number, FORCE_INTEGRATION_GROUP)):
                    identifier = self.aux.GetIthSubModelPartData(self.RigidFace_model_part, mesh_number, IDENTIFIER)
                    self.graph_forces[identifier] = open(str(self.DEM_parameters.problem_name) + "_" + str(identifier) + "_force_graph.grf", 'w')
                    self.graph_forces[identifier].write(str("#time").rjust(12)+" "+str("total_force[0]").rjust(13)+" "+str("total_force[1]").rjust(13)+" "+str("total_force[2]").rjust(13)+" "+str("total_moment[0]").rjust(13)+" "+str("total_moment[1]").rjust(13)+" "+str("total_moment[2]").rjust(13)+"\n")

        self.graph_forces = {}

        def open_balls_graph_files(self, spheres_model_part):
            #os.chdir(self.graphs_path)
            for mesh_number in range(0, self.spheres_model_part.NumberOfSubModelParts()):
                if (self.aux.GetIthSubModelPartData(self.spheres_model_part, mesh_number, FORCE_INTEGRATION_GROUP)):
                    identifier = self.aux.GetIthSubModelPartData(self.spheres_model_part, mesh_number, IDENTIFIER)
                    self.particle_graph_forces[identifier] = open(str(self.DEM_parameters.problem_name) + "_" + str(identifier) + "_particle_force_graph.grf", 'w');
                    self.particle_graph_forces[identifier].write(str("#time").rjust(12)+" "+str("total_force_x").rjust(13)+" "+str("total_force_y").rjust(13)+" "+str("total_force_z").rjust(13)+" "+str("total_disp_y").rjust(13)+"\n")

        def evaluate_computation_of_fem_results():

            self.spheres_model_part.ProcessInfo.SetValue(COMPUTE_FEM_RESULTS_OPTION, 0)
            elastic_forces            = Var_Translator(self.DEM_parameters.PostElasticForces)
            tangential_elastic_forces = Var_Translator(self.DEM_parameters.PostTangentialElasticForces)
            dem_pressure              = Var_Translator(self.DEM_parameters.PostPressure)
            if not (hasattr(self.DEM_parameters, "PostContactForces")):
                contact_forces = 0
            else:
                contact_forces = Var_Translator(self.DEM_parameters.PostShearStress)
            if not (hasattr(self.DEM_parameters, "PostShearStress")):
                shear_stress = 0
            else:
                shear_stress = Var_Translator(self.DEM_parameters.PostNodalArea)
            if not (hasattr(self.DEM_parameters, "PostNodalArea")):
                dem_nodal_area = 0
            else:
                dem_nodal_area = Var_Translator(self.DEM_parameters.PostNodalArea)
            integration_groups        = False
            if self.RigidFace_model_part.NumberOfSubModelParts() > 1:
                for mesh_number in range(0, self.RigidFace_model_part.NumberOfSubModelParts()):
                    if (self.aux.GetIthSubModelPartData(self.RigidFace_model_part, mesh_number, FORCE_INTEGRATION_GROUP)):
                        integration_groups = True
                        break
            if (elastic_forces or contact_forces or dem_pressure or tangential_elastic_forces or shear_stress or dem_nodal_area or integration_groups):
                self.spheres_model_part.ProcessInfo.SetValue(COMPUTE_FEM_RESULTS_OPTION, 1)

        self.particle_graph_forces = {}                    

        if not hasattr(self.DEM_parameters, "TestType"):
            open_graph_files(self, RigidFace_model_part)
            open_balls_graph_files(self,spheres_model_part)

        # SIMULATION SETTINGS
        self.bounding_box_enlargement_factor = self.DEM_parameters.BoundingBoxEnlargementFactor

        # MODEL
        self.domain_size = self.DEM_parameters.Dimension
        evaluate_computation_of_fem_results()
        
    def MoveAllMeshes(self, rigid_face_model_part, spheres_model_part, DEM_inlet_model_part, time, dt):
        
        self.mesh_motion.MoveAllMeshes(rigid_face_model_part, time, dt)
        self.mesh_motion.MoveAllMeshes(spheres_model_part, time, dt)
        self.mesh_motion.MoveAllMeshes(DEM_inlet_model_part, time, dt)
    
    def UpdateTimeInModelParts(self, spheres_model_part, rigid_face_model_part, cluster_model_part, time,dt,step):       
        spheres_model_part.ProcessInfo[TIME]          = time
        spheres_model_part.ProcessInfo[DELTA_TIME]    = dt
        spheres_model_part.ProcessInfo[TIME_STEPS]    = step

        rigid_face_model_part.ProcessInfo[TIME]       = time
        rigid_face_model_part.ProcessInfo[DELTA_TIME] = dt
        rigid_face_model_part.ProcessInfo[TIME_STEPS] = step

        cluster_model_part.ProcessInfo[TIME]          = time
        cluster_model_part.ProcessInfo[DELTA_TIME]    = dt
        cluster_model_part.ProcessInfo[TIME_STEPS]    = step  

    def close_graph_files(self, RigidFace_model_part):
        
        for mesh_number in range(0, self.RigidFace_model_part.NumberOfSubModelParts()):
            if (self.aux.GetIthSubModelPartData(self.RigidFace_model_part, mesh_number, FORCE_INTEGRATION_GROUP)):
                identifier = self.aux.GetIthSubModelPartData(self.RigidFace_model_part, mesh_number, IDENTIFIER)
                self.graph_forces[identifier].close()

    def close_balls_graph_files(self, spheres_model_part):
        
        for mesh_number in range(0, self.spheres_model_part.NumberOfSubModelParts()):
            if (self.aux.GetIthSubModelPartData(self.spheres_model_part, mesh_number, FORCE_INTEGRATION_GROUP)):
                identifier = self.aux.GetIthSubModelPartData(self.spheres_model_part, mesh_number, IDENTIFIER)
                self.particle_graph_forces[identifier].close()

    def PrintPoisson(self, model_part, DEM_parameters, filename, time):
        
        if (DEM_parameters.Dimension == 3):
            poisson, dummy, _ = PostUtilities().ComputePoisson(model_part)
        else:
            poisson, dummy, _ = PostUtilities().ComputePoisson2D(model_part)
            
        file_open = open(filename, 'a')
        data = str(time) + "  " + str(poisson) + "\n"
        file_open.write(data)
    
    
    def PrintGraph(self, time):

        if not hasattr(self.DEM_parameters, "TestType"):
            if (self.graph_counter == self.graph_frequency):
                self.graph_counter = 0
                if self.RigidFace_model_part.NumberOfSubModelParts() > 1:

                    for mesh_number in range(0, self.RigidFace_model_part.NumberOfSubModelParts()):
                        if (self.aux.GetIthSubModelPartData(self.RigidFace_model_part, mesh_number, FORCE_INTEGRATION_GROUP)):
                            mesh_nodes = self.aux.GetIthSubModelPartNodes(self.RigidFace_model_part,mesh_number)

                            total_force = Array3()

                            total_force[0] = 0.0
                            total_force[1] = 0.0
                            total_force[2] = 0.0

                            total_moment = Array3()

                            total_moment[0] = 0.0
                            total_moment[1] = 0.0
                            total_moment[2] = 0.0

                            rotation_center = self.aux.GetIthSubModelPartData(self.RigidFace_model_part, mesh_number, ROTATION_CENTER)

                            PostUtilities().IntegrationOfForces(mesh_nodes, total_force, rotation_center, total_moment)

                            identifier = self.aux.GetIthSubModelPartData(self.RigidFace_model_part, mesh_number, IDENTIFIER)
                            self.graph_forces[identifier].write(str("%.8g"%time).rjust(12)+" "+str("%.6g"%total_force[0]).rjust(13)+" "+str("%.6g"%total_force[1]).rjust(13)+" "+str("%.6g"%total_force[2]).rjust(13)+" "+str("%.6g"%total_moment[0]).rjust(13)+" "+str("%.6g"%total_moment[1]).rjust(13)+" "+str("%.6g"%total_moment[2]).rjust(13)+"\n")
                            self.graph_forces[identifier].flush()

            self.graph_counter += 1

    def FinalizeGraphs(self,RigidFace_model_part):

        if not hasattr(self.DEM_parameters, "TestType"):
            self.close_graph_files(RigidFace_model_part)

    def PrintBallsGraph(self, time):

        if not hasattr(self.DEM_parameters, "TestType"):
            if (self.balls_graph_counter == self.graph_frequency):
                self.balls_graph_counter = 0

                if self.spheres_model_part.NumberOfSubModelParts() > 1:
                    for mesh_number in range(0, self.spheres_model_part.NumberOfSubModelParts()):

                        if(self.aux.GetIthSubModelPartData(self.spheres_model_part, mesh_number, FORCE_INTEGRATION_GROUP)):
                            self.mesh_nodes = self.aux.GetIthSubModelPartNodes(self.spheres_model_part,mesh_number)
                            self.total_force_x = 0.0
                            self.total_force_y = 0.0
                            self.total_force_z = 0.0
                            self.total_disp_x = 0.0
                            self.total_disp_y = 0.0
                            self.total_disp_z = 0.0

                            for node in self.mesh_nodes:
                                force_node_x = node.GetSolutionStepValue(ELASTIC_FORCES)[0]
                                force_node_y = node.GetSolutionStepValue(ELASTIC_FORCES)[1]
                                force_node_z = node.GetSolutionStepValue(ELASTIC_FORCES)[2]
                                disp_node_x = node.GetSolutionStepValue(DISPLACEMENT)[0]
                                disp_node_y = node.GetSolutionStepValue(DISPLACEMENT)[1]
                                disp_node_z = node.GetSolutionStepValue(DISPLACEMENT)[2]
                                self.total_force_x += force_node_x
                                self.total_force_y += force_node_y
                                self.total_force_z += force_node_z
                                self.total_disp_x += disp_node_x
                                self.total_disp_y += disp_node_y
                                self.total_disp_z += disp_node_z

                            self.particle_graph_forces[self.spheres_model_part.GetMesh((mesh_number))[TOP]].write(str("%.8g"%time).rjust(12)+
                            " " + str("%.6g"%self.total_force_x).rjust(13) + " " + str("%.6g"%self.total_force_y).rjust(13)+
                            " " + str("%.6g"%self.total_force_z).rjust(13) + " " + str("%.6g"%self.total_disp_y).rjust(13)+
                            " " + str("%.6g"%self.total_disp_y).rjust(13)  + " " + str("%.6g"%self.total_disp_z).rjust(13) + "\n")
                            self.Flush(self.particle_graph_forces[self.spheres_model_part.GetMesh((mesh_number))[IDENTIFIER]])

            self.balls_graph_counter += 1

    def FinalizeBallsGraphs(self,spheres_model_part):

        if not hasattr(self.DEM_parameters, "TestType"):
            self.close_balls_graph_files(spheres_model_part)


    def ApplyNodalRotation(self,time):

            #if (time < 0.5e-2 ) :
            if (time < 3.8e-5):

                #while ( time < self.DEM_parameters.FinalTime):
                #print("TIME STEP BEGINS.  STEP:"+str(time)+"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")

                d0 = 1.694
                avance = 0.50100
                distance = (d0 - avance)
                w = 62.8
                distance = 2

                vx = -distance * w * math.sin(w * time)
                #vx = distance * rpm * math.cos(1.0)
                vy = distance * w * math.cos(w * time)
                #vz = - distance * rpm * math.sin(1.0)

                for mesh_number in range(0, self.spheres_model_part.NumberOfSubModelParts()):
                    if(self.aux.GetIthSubModelPartData(self.spheres_model_part, mesh_number, FORCE_INTEGRATION_GROUP)):
                        self.mesh_nodes = self.aux.GetIthSubModelPartNodes(self.spheres_model_part,mesh_number)

                        for node in self.mesh_nodes:

                            node.SetSolutionStepValue(VELOCITY_X, vx)
                            node.SetSolutionStepValue(VELOCITY_Y, vy)
                            node.Fix(VELOCITY_X)
                            node.Fix(VELOCITY_Y)
            else:

                 d0 = 1.694
                 avance = 0.50100
                 distance = (d0 - avance)
                 w = 62.8
                 distance = 2

                 vx = -distance * w * math.sin(w * time)
                 #vx = distance * rpm * math.cos(1.0)
                 vy = distance * w * math.cos(w * time)
                 #vz = - distance * rpm * math.sin(1.0)
                 radius = 1.0001

                 for mesh_number in range(0, self.spheres_model_part.NumberOfSubModelParts()):
                     if(self.aux.GetIthSubModelPartData(self.spheres_model_part, mesh_number, FORCE_INTEGRATION_GROUP)):
                        self.mesh_nodes = self.aux.GetIthSubModelPartNodes(self.spheres_model_part,mesh_number)

                        for node in self.mesh_nodes:

                            node.SetSolutionStepValue(RADIUS, radius)

                 for mesh_number in range(0, self.spheres_model_part.NumberOfSubModelParts()):
                     if(self.aux.GetIthSubModelPartData(self.spheres_model_part, mesh_number, FORCE_INTEGRATION_GROUP)):
                         self.mesh_nodes = self.aux.GetIthSubModelPartNodes(self.spheres_model_part,mesh_number)

                         for node in self.mesh_nodes:

                             node.SetSolutionStepValue(VELOCITY_X, vx)
                             node.SetSolutionStepValue(VELOCITY_Y, vy)
                             node.Fix(VELOCITY_X)
                             node.Fix(VELOCITY_Y)        

class Report(object):

    def __init__(self):
        pass

    def Prepare(self,timer,control_time):
        self.initial_pr_time      = timer.clock()
        self.initial_re_time      = timer.time()
        self.prev_time            = 0.0
        self.total_steps_expected = 0
        self.control_time         = control_time
        self.first_print          = True

    def BeginReport(self, timer):

        report = "Main loop starting..." + "\n" + "Total number of TIME STEPs expected in the calculation: " + str(self.total_steps_expected) + "\n"

        return report

    def StepiReport(self, timer, time, step):

        incremental_time = (timer.time() - self.initial_re_time) - self.prev_time

        report = ""

        if (incremental_time > self.control_time):

            percentage = 100 * (float(step) / self.total_steps_expected)
            elapsed_time = timer.time() - self.initial_re_time

            report = report + "Real time calculation: " + str(elapsed_time) + " seconds" + "\n"\
                            + "In minutes: " + str(elapsed_time / 60.0) + " minutes" + "\n"\
                            + "In hours: " + str(elapsed_time / 3600.0) + " hours" + "\n"\
                            + "Simulation time: " + str(time) + " seconds" + "\n"\
                            + "%s %.5f %s"%("Percentage Completed: ", percentage,"%") + "\n"\
                            + "Computed time steps: " + str(step) + " out of " + str(self.total_steps_expected) + "\n"
            
            self.prev_time  = (timer.time() - self.initial_re_time)

        if ((timer.time() - self.initial_re_time > 60) and self.first_print == True and step != 0):

            self.first_print = False
            estimated_sim_duration = 60.0 * (self.total_steps_expected / step) # seconds

            report = report + "The total estimated computation time is " + str(estimated_sim_duration) + " seconds" + "\n"\
                + "In minutes: " + str(estimated_sim_duration / 60.0)    + " minutes" + "\n"\
                + "In hours:   " + str(estimated_sim_duration / 3600.0)  + " hours" + "\n"\
                + "In days:    " + str(estimated_sim_duration / 86400.0) + " days" + "\n"

            if ((estimated_sim_duration / 86400.0) > 2.0):
                report = report +"WARNING: VERY LONG CALCULATION......!!!!!!" + "\n"

        return report

    def FinalReport(self, timer):
        elapsed_pr_time = timer.clock() - self.initial_pr_time
        elapsed_re_time = timer.time()  - self.initial_re_time

        report = "Calculation ends at instant: "              + str(timer.time()) + "\n"\
            + "Calculation ends at processing time instant: " + str(timer.clock()) + "\n"\
            + "Elapsed processing time: "                     + str(elapsed_pr_time) + "\n"\
            + "Elapsed real time: "                           + str(elapsed_re_time) +"\n"

        report = report + "ANALYSIS COMPLETED" + "\n"

        return report


class MaterialTest(object):

    def __init__(self):
        pass

    def Initialize(self, DEM_parameters, procedures, solver, graphs_path, post_path, spheres_model_part, rigid_face_model_part):

        if not hasattr(DEM_parameters, "TestType"):
            self.TestType = "None"
        else:
            self.TestType = DEM_parameters.TestType

        if (self.TestType != "None"):
            self.script = DEM_material_test_script.MaterialTest(DEM_parameters, procedures, solver, graphs_path, post_path, spheres_model_part, rigid_face_model_part)
            self.script.Initialize()

    def PrepareDataForGraph(self):
        if (self.TestType != "None"):
            self.script.PrepareDataForGraph()

    def MeasureForcesAndPressure(self):
        if (self.TestType != "None"):
            self.script.MeasureForcesAndPressure()

    def PrintGraph(self, time):
        if (self.TestType != "None"):
            self.script.PrintGraph(time)

    def FinalizeGraphs(self):
        if (self.TestType != "None"):
            self.script.FinalizeGraphs()


    def PrintChart(self):
        if (self.TestType != "None"):
            self.script.PrintChart()


class MultifileList(object):

    def __init__(self, post_path, name, step):
        os.chdir(post_path)
        self.index = 0
        self.step = step
        self.name = name
        self.file = open("_list_" + self.name + "_" + str(step) + ".post.lst","w")


class DEMIo(object):

    def __init__(self, DEM_parameters, post_path):

        self.post_path = post_path
        self.mixed_model_part                                     = ModelPart("Mixed_Part")
        self.mixed_spheres_and_clusters_model_part                = ModelPart("MixedSpheresAndClustersPart")
        self.mixed_spheres_not_in_cluster_and_clusters_model_part = ModelPart("MixedSpheresNotInClusterAndClustersPart")
        
        # Printing variables
        self.DEM_parameters = DEM_parameters
        self.global_variables                          = []
        self.spheres_and_clusters_variables            = []
        self.spheres_and_clusters_local_axis_variables = []
        self.spheres_not_in_cluster_and_clusters_variables = []
        self.spheres_not_in_cluster_and_clusters_local_axis_variables = []
        self.spheres_variables                         = []
        self.spheres_local_axis_variables              = []
        self.fem_boundary_variables                    = []
        self.clusters_variables                        = []
        self.contact_variables                         = []
        self.multifilelists                            = []

        # Reading Post options from DEM_parameters
        self.PostDisplacement             = getattr(self.DEM_parameters, "PostDisplacement", 0)
        self.PostVelocity                 = getattr(self.DEM_parameters, "PostVelocity", 0)
        self.PostTotalForces              = getattr(self.DEM_parameters, "PostTotalForces", 0)
        self.PostNonDimensionalVolumeWear = getattr(self.DEM_parameters, "PostNonDimensionalVolumeWear", 0)
        self.PostAppliedForces            = getattr(self.DEM_parameters, "PostAppliedForces", 0)
        self.PostDampForces               = getattr(self.DEM_parameters, "PostDampForces", 0)
        self.PostRadius                   = getattr(self.DEM_parameters, "PostRadius", 0)
        self.PostExportId                 = getattr(self.DEM_parameters, "PostExportId", 0)
        self.PostSkinSphere               = getattr(self.DEM_parameters, "PostSkinSphere", 0)
        self.PostAngularVelocity          = getattr(self.DEM_parameters, "PostAngularVelocity", 0)
        self.PostParticleMoment           = getattr(self.DEM_parameters, "PostParticleMoment", 0)
        self.PostEulerAngles              = getattr(self.DEM_parameters, "PostEulerAngles", 0)
        self.PostRollingResistanceMoment  = getattr(self.DEM_parameters, "PostRollingResistanceMoment", 0)
        self.PostLocalContactForce        = getattr(self.DEM_parameters, "PostLocalContactForce", 0)
        self.PostFailureCriterionState    = getattr(self.DEM_parameters, "PostFailureCriterionState", 0)
        self.PostContactFailureId         = getattr(self.DEM_parameters, "PostContactFailureId", 0)
        self.PostContactTau               = getattr(self.DEM_parameters, "PostContactTau", 0)
        self.PostContactSigma             = getattr(self.DEM_parameters, "PostContactSigma", 0)
        self.PostMeanContactArea          = getattr(self.DEM_parameters, "PostMeanContactArea", 0)
        self.PostElasticForces            = getattr(self.DEM_parameters, "PostElasticForces", 0)
        self.PostContactForces            = getattr(self.DEM_parameters, "PostContactForces", 0)
        self.PostRigidElementForces       = getattr(self.DEM_parameters, "PostRigidElementForces", 0)
        self.PostPressure                 = getattr(self.DEM_parameters, "PostPressure", 0)
        self.PostTangentialElasticForces  = getattr(self.DEM_parameters, "PostTangentialElasticForces", 0)
        self.PostShearStress              = getattr(self.DEM_parameters, "PostShearStress", 0)
        self.PostNodalArea                = getattr(self.DEM_parameters, "PostNodalArea", 0)
        self.VelTrapGraphExportFreq       = getattr(self.DEM_parameters, "VelTrapGraphExportFreq", 0)
        self.PostTemperature              = getattr(self.DEM_parameters, "PostTemperature", 0)
        self.PostHeatFlux                 = getattr(self.DEM_parameters, "PostHeatFlux", 0)
        self.PostNeighbourSize            = getattr(self.DEM_parameters, "PostNeighbourSize", 0)

        if not (hasattr(self.DEM_parameters, "PostBoundingBox")):
            self.PostBoundingBox = 0
        else:
            self.PostBoundingBox = getattr(self.DEM_parameters, "PostBoundingBox", 0)
            
        self.automatic_bounding_box_option = Var_Translator(self.DEM_parameters.AutomaticBoundingBoxOption)
        self.b_box_minX = self.DEM_parameters.BoundingBoxMinX
        self.b_box_minY = self.DEM_parameters.BoundingBoxMinY
        self.b_box_minZ = self.DEM_parameters.BoundingBoxMinZ
        self.b_box_maxX = self.DEM_parameters.BoundingBoxMaxX
        self.b_box_maxY = self.DEM_parameters.BoundingBoxMaxY
        self.b_box_maxZ = self.DEM_parameters.BoundingBoxMaxZ

        self.continuum_element_types = ["SphericContPartDEMElement3D", "CylinderContPartDEMElement2D", "IceContPartDEMElement3D"]
        
        self.multifiles = (
            MultifileList(self.post_path, DEM_parameters.problem_name, 1),
            MultifileList(self.post_path, DEM_parameters.problem_name, 2),
            MultifileList(self.post_path, DEM_parameters.problem_name, 5),
            MultifileList(self.post_path, DEM_parameters.problem_name,10),
            MultifileList(self.post_path, DEM_parameters.problem_name,20),
            MultifileList(self.post_path, DEM_parameters.problem_name,50),
            )
            
        self.SetMultifileLists(self.multifiles)
        
        # Ice
        if (hasattr(self.DEM_parameters, "PostVirtualSeaSurfaceX1")):
            self.SeaSurfaceX1 = self.DEM_parameters.PostVirtualSeaSurfaceX1
            self.SeaSurfaceY1 = self.DEM_parameters.PostVirtualSeaSurfaceY1
            self.SeaSurfaceX2 = self.DEM_parameters.PostVirtualSeaSurfaceX2
            self.SeaSurfaceY2 = self.DEM_parameters.PostVirtualSeaSurfaceY2
            self.SeaSurfaceX3 = self.DEM_parameters.PostVirtualSeaSurfaceX3
            self.SeaSurfaceY3 = self.DEM_parameters.PostVirtualSeaSurfaceY3
            self.SeaSurfaceX4 = self.DEM_parameters.PostVirtualSeaSurfaceX4
            self.SeaSurfaceY4 = self.DEM_parameters.PostVirtualSeaSurfaceY4

    def KRATOSprint(self,message):
        print(message)
        self.Flush(sys.stdout)
        
    def Flush(self,a):
        a.flush()
        
    def ShowPrintingResultsOnScreen(self, all_model_parts):
        self.KRATOSprint("*******************  PRINTING RESULTS FOR GID  ***************************")
        self.KRATOSprint("                        ("+ str(all_model_parts.spheres_model_part.NumberOfElements(0)) + " elements)")
        self.KRATOSprint("                        ("+ str(all_model_parts.spheres_model_part.NumberOfNodes(0)) + " nodes)")
        self.KRATOSprint("")
        
    def Initialize(self, DEM_parameters):
        self.AddGlobalVariables()
        self.AddSpheresVariables()
        self.AddSpheresAndClustersVariables()
        self.AddSpheresNotInClusterAndClustersVariables()
        self.AddFEMBoundaryVariables()
        self.AddClusterVariables()
        self.AddContactVariables()
        self.AddMpiVariables()
        self.Configure(DEM_parameters.problem_name, DEM_parameters.OutputFileType, DEM_parameters.Multifile, DEM_parameters.ContactMeshOption)
        self.SetOutputName(DEM_parameters.problem_name)

    def PushPrintVar(self, variable, name, print_list):
        if (Var_Translator(variable)):
            print_list.append(name)

    def AddGlobalVariables(self):
        self.PushPrintVar(self.PostDisplacement,     DISPLACEMENT,                 self.global_variables)
        self.PushPrintVar(self.PostVelocity,         VELOCITY,                     self.global_variables)
        self.PushPrintVar(self.PostTotalForces,      TOTAL_FORCES,                 self.global_variables)

    def AddSpheresAndClustersVariables(self):  # variables common to spheres and clusters
        self.PushPrintVar(self.PostAppliedForces,       EXTERNAL_APPLIED_FORCE,  self.spheres_and_clusters_variables)
        self.PushPrintVar(self.PostAppliedForces,       EXTERNAL_APPLIED_MOMENT, self.spheres_and_clusters_variables)
        self.PushPrintVar(self.PostRigidElementForces,  RIGID_ELEMENT_FORCE,     self.spheres_and_clusters_variables)
        #if (Var_Translator(self.DEM_parameters.RotationOption)):  # xapuza
        if (Var_Translator(self.DEM_parameters.PostAngularVelocity)):
            self.PushPrintVar(self.PostAngularVelocity, ANGULAR_VELOCITY,        self.spheres_and_clusters_variables)
        if (Var_Translator(self.DEM_parameters.PostParticleMoment)):       
            self.PushPrintVar(self.PostParticleMoment,  PARTICLE_MOMENT,         self.spheres_and_clusters_variables)
            
    def AddSpheresNotInClusterAndClustersVariables(self):  # variables common to spheres and clusters
        #if (Var_Translator(self.DEM_parameters.RotationOption)):  # xapuza
        if (Var_Translator(self.DEM_parameters.PostEulerAngles)):
            self.PushPrintVar(self.PostEulerAngles,     EULER_ANGLES,            self.spheres_not_in_cluster_and_clusters_local_axis_variables)

    def AddSpheresVariables(self):
        self.PushPrintVar(self.PostDampForces,       DAMP_FORCES,                  self.spheres_variables)
        self.PushPrintVar(self.PostRadius,           RADIUS,                       self.spheres_variables)
        self.PushPrintVar(self.PostExportId,         EXPORT_ID,                    self.spheres_variables)
        self.PushPrintVar(self.PostTemperature,      TEMPERATURE,                  self.spheres_variables)
        self.PushPrintVar(self.PostHeatFlux,         HEATFLUX,                     self.spheres_variables)

        #self.PushPrintVar(                        1, DELTA_DISPLACEMENT,           self.spheres_variables)  # Debugging
        #self.PushPrintVar(                        1, PARTICLE_ROTATION_ANGLE,      self.spheres_variables)  # Debugging
        if (hasattr(self.DEM_parameters, "PostRollingResistanceMoment")):
            if (Var_Translator(self.DEM_parameters.RotationOption)):
                if (Var_Translator(self.DEM_parameters.RollingFrictionOption)):
                    self.PushPrintVar( self.PostRollingResistanceMoment, ROLLING_RESISTANCE_MOMENT, self.spheres_variables)

        if (hasattr(self.DEM_parameters, "PostSkinSphere")):
            if (Var_Translator(self.DEM_parameters.PostSkinSphere)):
                self.PushPrintVar(self.PostSkinSphere,       SKIN_SPHERE,              self.spheres_variables)

        if (hasattr(self.DEM_parameters, "PostNeighbourSize")):
            if (Var_Translator(self.DEM_parameters.PostNeighbourSize)):
                self.PushPrintVar(self.PostNeighbourSize,       NEIGHBOUR_SIZE,              self.spheres_variables)

        # NANO
        if self.DEM_parameters.ElementType == "SwimmingNanoParticle":
            self.PushPrintVar(self.PostHeatFlux, CATION_CONCENTRATION, self.spheres_variables)

        if (hasattr(self.DEM_parameters, "PostStressStrainOption")):
            if (Var_Translator(self.DEM_parameters.PostStressStrainOption)):
                self.PushPrintVar(1, REPRESENTATIVE_VOLUME, self.spheres_variables)
                self.PushPrintVar(1, DEM_STRESS_TENSOR,     self.spheres_variables)
                self.PushPrintVar(1, FORCE_REACTION,        self.spheres_variables)
                self.PushPrintVar(1, MOMENT_REACTION,       self.spheres_variables)

        if (hasattr(self.DEM_parameters, "PostPoissonRatio")):
            if Var_Translator(self.DEM_parameters.PostPoissonRatio):
                self.PushPrintVar(self.DEM_parameters.PostPoissonRatio, POISSON_VALUE, self.spheres_variables)

    def AddFEMBoundaryVariables(self):
        self.PushPrintVar(self.PostElasticForces,            ELASTIC_FORCES, self.fem_boundary_variables)
        self.PushPrintVar(self.PostContactForces,            CONTACT_FORCES, self.fem_boundary_variables)
        self.PushPrintVar(self.PostPressure,                 DEM_PRESSURE, self.fem_boundary_variables)
        self.PushPrintVar(self.PostTangentialElasticForces,  TANGENTIAL_ELASTIC_FORCES, self.fem_boundary_variables)
        self.PushPrintVar(self.PostShearStress,              SHEAR_STRESS, self.fem_boundary_variables)
        self.PushPrintVar(self.PostNodalArea,                DEM_NODAL_AREA, self.fem_boundary_variables)
        if (Var_Translator(self.PostNonDimensionalVolumeWear)):
            self.PushPrintVar(1,                             NON_DIMENSIONAL_VOLUME_WEAR, self.fem_boundary_variables)
            self.PushPrintVar(1,                             IMPACT_WEAR,                 self.fem_boundary_variables)

    def AddClusterVariables(self):
        self.PushPrintVar(self.PostRadius, CHARACTERISTIC_LENGTH, self.clusters_variables)

    def AddContactVariables(self):
        # Contact Elements Variables
        if (self.DEM_parameters.ElementType in self.continuum_element_types):
            if (Var_Translator(self.DEM_parameters.ContactMeshOption)):
                self.PushPrintVar(self.PostLocalContactForce,     LOCAL_CONTACT_FORCE,     self.contact_variables)
                self.PushPrintVar(self.PostFailureCriterionState, FAILURE_CRITERION_STATE, self.contact_variables)
                self.PushPrintVar(self.PostContactFailureId,      CONTACT_FAILURE,         self.contact_variables)
                self.PushPrintVar(self.PostContactTau,            CONTACT_TAU,             self.contact_variables)
                self.PushPrintVar(self.PostContactSigma,          CONTACT_SIGMA,           self.contact_variables)
                self.PushPrintVar(self.PostMeanContactArea,       MEAN_CONTACT_AREA,       self.contact_variables)

    def AddMpiVariables(self):
        pass

    def Configure(self, problem_name, encoding, file_system, contact_mesh_option):
        self.problem_name = problem_name

        if (encoding == "Binary"):
            self.encoding = GiDPostMode.GiD_PostBinary
        else:
            self.encoding = GiDPostMode.GiD_PostAscii

        if (self.DEM_parameters.Multifile == "multiple_files"):
            self.filesystem = MultiFileFlag.MultipleFiles
        else:
            self.filesystem = MultiFileFlag.SingleFile

        self.deformed_mesh_flag  = WriteDeformedMeshFlag.WriteDeformed
        self.write_conditions    = WriteConditionsFlag.WriteConditions
        self.contact_mesh_option = contact_mesh_option

        self.gid_io = GidIO(self.problem_name,
                            self.encoding,
                            self.filesystem,
                            self.deformed_mesh_flag,
                            self.write_conditions)

        self.post_utility = PostUtilities()

    def SetOutputName(self,name):
        self.gid_io.ChangeOutputName(name)

    def SetMultifileLists(self,multifile_list):
        for mfilelist in multifile_list:
            self.multifilelists.append(mfilelist)

        for mfilelist in self.multifilelists:
            mfilelist.file.write("Multiple\n")
            mfilelist.index = 1

    def PrintMultifileLists(self, time, post_path):
        for mfilelist in self.multifilelists:

            if mfilelist.index == mfilelist.step:
                
                if (self.encoding == GiDPostMode.GiD_PostBinary):
                    mfilelist.file.write(self.GetMultiFileListName(mfilelist.name)+"_"+"%.12g"%time+".post.bin\n")
                else:
                    mfilelist.file.write(self.GetMultiFileListName(mfilelist.name)+"_"+"%.12g"%time+".post.msh\n")
                    mfilelist.file.write(self.GetMultiFileListName(mfilelist.name)+"_"+"%.12g"%time+".post.res\n")
                self.Flush(mfilelist.file)
                mfilelist.index = 0
            mfilelist.index += 1
            
    def GetMultiFileListName(self, name):
        return name

    def CloseMultifiles(self):
        for mfilelist in self.multifilelists:
            mfilelist.file.close()

    def InitializeMesh(self, all_model_parts): #MIQUEL MAPPING
        if (self.filesystem == MultiFileFlag.SingleFile):

            self.post_utility.AddModelPartToModelPart(self.mixed_model_part, spheres_model_part)

            if (self.contact_mesh_option == "ON"):
                self.post_utility.AddModelPartToModelPart(self.mixed_model_part, contact_model_part)

            self.post_utility.AddModelPartToModelPart(self.mixed_model_part, rigid_face_model_part)
            self.post_utility.AddModelPartToModelPart(self.mixed_model_part, cluster_model_part)
            self.post_utility.AddModelPartToModelPart(self.mixed_spheres_and_clusters_model_part, spheres_model_part)
            self.post_utility.AddModelPartToModelPart(self.mixed_spheres_and_clusters_model_part, cluster_model_part)
            
            self.post_utility.AddSpheresNotBelongingToClustersToMixModelPart(self.mixed_spheres_not_in_cluster_and_clusters_model_part, spheres_model_part)
            self.post_utility.AddModelPartToModelPart(self.mixed_spheres_not_in_cluster_and_clusters_model_part, cluster_model_part)

            self.gid_io.InitializeMesh(0.0)
            self.gid_io.WriteMesh(rigid_face_model_part.GetCommunicator().LocalMesh())
            self.gid_io.WriteClusterMesh(cluster_model_part.GetCommunicator().LocalMesh())
            if (self.DEM_parameters.ElementType == "CylinderContPartDEMElement2D"):
                self.gid_io.WriteCircleMesh(spheres_model_part.GetCommunicator().LocalMesh())
            else:
                self.gid_io.WriteSphereMesh(spheres_model_part.GetCommunicator().LocalMesh())

            if (self.contact_mesh_option == "ON"):
                self.gid_io.WriteMesh(contact_model_part.GetCommunicator().LocalMesh())

            self.gid_io.FinalizeMesh()
            self.gid_io.InitializeResults(0.0, self.mixed_model_part.GetCommunicator().LocalMesh())
            #self.gid_io.InitializeResults(0.0, mixed_spheres_and_clusters_model_part.GetCommunicator().LocalMesh())

    def InitializeResults(self, spheres_model_part, rigid_face_model_part, cluster_model_part, contact_model_part, mapping_model_part, creator_destructor, dem_fem_search, time, bounding_box_time_limits): #MIQUEL MAPPING

        if (self.filesystem == MultiFileFlag.MultipleFiles):
            self.mixed_model_part.Elements.clear()
            self.mixed_model_part.Nodes.clear()
            self.mixed_spheres_and_clusters_model_part.Elements.clear()
            self.mixed_spheres_and_clusters_model_part.Nodes.clear()
            self.mixed_spheres_not_in_cluster_and_clusters_model_part.Elements.clear()
            self.mixed_spheres_not_in_cluster_and_clusters_model_part.Nodes.clear()

            self.post_utility.AddModelPartToModelPart(self.mixed_model_part, spheres_model_part)
            if (self.contact_mesh_option == "ON"):
                self.post_utility.AddModelPartToModelPart(self.mixed_model_part, contact_model_part)
            self.post_utility.AddModelPartToModelPart(self.mixed_model_part, rigid_face_model_part)
            self.post_utility.AddModelPartToModelPart(self.mixed_model_part, cluster_model_part)
            self.post_utility.AddModelPartToModelPart(self.mixed_spheres_and_clusters_model_part, spheres_model_part)
            self.post_utility.AddModelPartToModelPart(self.mixed_spheres_and_clusters_model_part, cluster_model_part)
            
            self.post_utility.AddSpheresNotBelongingToClustersToMixModelPart(self.mixed_spheres_not_in_cluster_and_clusters_model_part, spheres_model_part)
            self.post_utility.AddModelPartToModelPart(self.mixed_spheres_not_in_cluster_and_clusters_model_part, cluster_model_part)

            self.gid_io.InitializeMesh(time)
            if (self.DEM_parameters.ElementType == "CylinderContPartDEMElement2D"):
                self.gid_io.WriteCircleMesh(spheres_model_part.GetCommunicator().LocalMesh())
            else:
                self.gid_io.WriteSphereMesh(spheres_model_part.GetCommunicator().LocalMesh())
            if (self.contact_mesh_option == "ON"):
                #We overwrite the Id of the properties 0 not to overlap with other entities that use layer 0 for PRINTING
                contact_model_part.GetProperties(0)[0].Id = 9184
                self.gid_io.WriteMesh(contact_model_part.GetCommunicator().LocalMesh())

            self.gid_io.WriteMesh(rigid_face_model_part.GetCommunicator().LocalMesh())
            self.gid_io.WriteClusterMesh(cluster_model_part.GetCommunicator().LocalMesh())
            
            #Compute and print the graphical bounding box if active in time
            if ((getattr(self.DEM_parameters, "BoundingBoxOption", 0) == "ON") and (time >= bounding_box_time_limits[0] and time <= bounding_box_time_limits[1])):
                self.ComputeAndPrintBoundingBox(spheres_model_part, rigid_face_model_part, contact_model_part, creator_destructor)
                        
            # Ice. Printing a virtual sea surface
            if (hasattr(self.DEM_parameters, "PostVirtualSeaSurfaceX1")):
                self.ComputeAndPrintSeaSurface(spheres_model_part, rigid_face_model_part)
                
            #self.ComputeAndPrintDEMFEMSearchBinBoundingBox(spheres_model_part, rigid_face_model_part, dem_fem_search)#MSIMSI

            self.gid_io.FinalizeMesh()                
            self.gid_io.InitializeResults(time, self.mixed_model_part.GetCommunicator().LocalMesh())                        
            #self.gid_io.InitializeResults(time, mixed_spheres_and_clusters_model_part.GetCommunicator().LocalMesh())

    def FinalizeMesh(self):
        if (self.filesystem == MultiFileFlag.SingleFile):
            self.gid_io.FinalizeResults()

    def FinalizeResults(self):
        if (self.filesystem == MultiFileFlag.MultipleFiles):
            self.gid_io.FinalizeResults()

    def PrintingGlobalVariables(self, export_model_part, time):
        for variable in self.global_variables:
            self.gid_io.WriteNodalResults(variable, export_model_part.Nodes, time, 0)

    def PrintingSpheresAndClustersVariables(self, export_model_part, time):
        for variable in self.spheres_and_clusters_variables:
            self.gid_io.WriteNodalResults(variable, export_model_part.Nodes, time, 0)
        for variable in self.spheres_and_clusters_local_axis_variables:
            self.gid_io.WriteLocalAxesOnNodes(variable, export_model_part.Nodes, time, 0)
            
    def PrintingSpheresNotInClusterAndClustersVariables(self, export_model_part, time):
        for variable in self.spheres_not_in_cluster_and_clusters_variables:
            self.gid_io.WriteNodalResults(variable, export_model_part.Nodes, time, 0)
        for variable in self.spheres_not_in_cluster_and_clusters_local_axis_variables:
            self.gid_io.WriteLocalAxesOnNodes(variable, export_model_part.Nodes, time, 0)

    def PrintingSpheresVariables(self, export_model_part, time):
        for variable in self.spheres_variables:
            self.gid_io.WriteNodalResults(variable, export_model_part.Nodes, time, 0)
        for variable in self.spheres_local_axis_variables:
            self.gid_io.WriteLocalAxesOnNodes(variable, export_model_part.Nodes, time, 0)

    def PrintingFEMBoundaryVariables(self, export_model_part, time):
        for variable in self.fem_boundary_variables:
            self.gid_io.WriteNodalResults(variable, export_model_part.Nodes, time, 0)

    def PrintingClusterVariables(self, export_model_part, time):
        for variable in self.clusters_variables:
            self.gid_io.WriteNodalResults(variable, export_model_part.Nodes, time, 0)

    def PrintingContactElementsVariables(self, export_model_part, time):
        if (self.contact_mesh_option == "ON"):
            for variable in self.contact_variables:
                self.gid_io.PrintOnGaussPoints(variable, export_model_part, time)

    def PrintResults(self, spheres_model_part, rigid_face_model_part, cluster_model_part, contact_model_part, mapping_model_part, creator_destructor, dem_fem_search, time, bounding_box_time_limits):
        if (self.filesystem == MultiFileFlag.MultipleFiles):
            self.InitializeResults(spheres_model_part,
                                   rigid_face_model_part,
                                   cluster_model_part,
                                   contact_model_part,
                                   mapping_model_part,
                                   creator_destructor,
                                   dem_fem_search,
                                   time,
                                   bounding_box_time_limits)

        self.PrintingGlobalVariables(self.mixed_model_part, time)
        self.PrintingSpheresAndClustersVariables(self.mixed_spheres_and_clusters_model_part, time)
        self.PrintingSpheresNotInClusterAndClustersVariables(self.mixed_spheres_not_in_cluster_and_clusters_model_part, time)
        self.PrintingSpheresVariables(spheres_model_part, time)
        self.PrintingFEMBoundaryVariables(rigid_face_model_part, time)
        self.PrintingClusterVariables(cluster_model_part, time)
        self.PrintingContactElementsVariables(contact_model_part, time)
        
        self.mixed_model_part.Elements.clear() #to remove the shared pointers that remain and prevent objects from being removed
        self.mixed_model_part.Nodes.clear() #to remove the shared pointers that remain and prevent objects from being removed
        self.mixed_spheres_and_clusters_model_part.Elements.clear() #to remove the shared pointers that remain and prevent objects from being removed
        self.mixed_spheres_and_clusters_model_part.Nodes.clear() #to remove the shared pointers that remain and prevent objects from being removed
        self.mixed_spheres_not_in_cluster_and_clusters_model_part.Elements.clear() #to remove the shared pointers that remain and prevent objects from being removed
        self.mixed_spheres_not_in_cluster_and_clusters_model_part.Nodes.clear() #to remove the shared pointers that remain and prevent objects from being removed

        if (self.filesystem == MultiFileFlag.MultipleFiles):
            self.FinalizeResults()

    def ComputeAndPrintBoundingBox(self, spheres_model_part, rigid_face_model_part, contact_model_part, creator_destructor):

        bounding_box_model_part = ModelPart("BoundingBoxPart") # Creation of bounding box's model part

        max_node_Id        = ParticleCreatorDestructor().FindMaxNodeIdInModelPart(spheres_model_part)
        max_FEM_node_Id    = ParticleCreatorDestructor().FindMaxNodeIdInModelPart(rigid_face_model_part)
        max_element_Id     = ParticleCreatorDestructor().FindMaxElementIdInModelPart(spheres_model_part)
        max_FEM_element_Id = ParticleCreatorDestructor().FindMaxElementIdInModelPart(rigid_face_model_part)
        max_contact_element_Id = ParticleCreatorDestructor().FindMaxElementIdInModelPart(contact_model_part)

        if (max_FEM_node_Id > max_node_Id):
            max_node_Id = max_FEM_node_Id

        if (max_FEM_element_Id > max_element_Id):
            max_element_Id = max_FEM_element_Id
        
        if (max_contact_element_Id > max_element_Id):
            max_element_Id = max_contact_element_Id

        BBMaxX = creator_destructor.GetHighNode()[0]
        BBMaxY = creator_destructor.GetHighNode()[1]
        BBMaxZ = creator_destructor.GetHighNode()[2]
        BBMinX = creator_destructor.GetLowNode()[0]
        BBMinY = creator_destructor.GetLowNode()[1]
        BBMinZ = creator_destructor.GetLowNode()[2]

        # BB Nodes:
        node1 = bounding_box_model_part.CreateNewNode(max_node_Id + 1, BBMinX, BBMinY, BBMinZ)
        node2 = bounding_box_model_part.CreateNewNode(max_node_Id + 2, BBMaxX, BBMinY, BBMinZ)
        node3 = bounding_box_model_part.CreateNewNode(max_node_Id + 3, BBMaxX, BBMaxY, BBMinZ)
        node4 = bounding_box_model_part.CreateNewNode(max_node_Id + 4, BBMinX, BBMaxY, BBMinZ)
        node5 = bounding_box_model_part.CreateNewNode(max_node_Id + 5, BBMinX, BBMinY, BBMaxZ)
        node6 = bounding_box_model_part.CreateNewNode(max_node_Id + 6, BBMaxX, BBMinY, BBMaxZ)
        node7 = bounding_box_model_part.CreateNewNode(max_node_Id + 7, BBMaxX, BBMaxY, BBMaxZ)
        node8 = bounding_box_model_part.CreateNewNode(max_node_Id + 8, BBMinX, BBMaxY, BBMaxZ)

        props = Properties(10000) # Property 10000 corresponds to black colour
        
        # BB Elements:
        bounding_box_model_part.CreateNewCondition("RigidEdge3D", max_element_Id +  1, [node1.Id, node4.Id], props)
        bounding_box_model_part.CreateNewCondition("RigidEdge3D", max_element_Id +  2, [node4.Id, node8.Id], props)
        bounding_box_model_part.CreateNewCondition("RigidEdge3D", max_element_Id +  3, [node8.Id, node5.Id], props)
        bounding_box_model_part.CreateNewCondition("RigidEdge3D", max_element_Id +  4, [node5.Id, node1.Id], props)
        bounding_box_model_part.CreateNewCondition("RigidEdge3D", max_element_Id +  5, [node1.Id, node2.Id], props)
        bounding_box_model_part.CreateNewCondition("RigidEdge3D", max_element_Id +  6, [node3.Id, node4.Id], props)
        bounding_box_model_part.CreateNewCondition("RigidEdge3D", max_element_Id +  7, [node7.Id, node8.Id], props)
        bounding_box_model_part.CreateNewCondition("RigidEdge3D", max_element_Id +  8, [node5.Id, node6.Id], props)
        bounding_box_model_part.CreateNewCondition("RigidEdge3D", max_element_Id +  9, [node6.Id, node2.Id], props)
        bounding_box_model_part.CreateNewCondition("RigidEdge3D", max_element_Id + 10, [node2.Id, node3.Id], props)
        bounding_box_model_part.CreateNewCondition("RigidEdge3D", max_element_Id + 11, [node3.Id, node7.Id], props)
        bounding_box_model_part.CreateNewCondition("RigidEdge3D", max_element_Id + 12, [node7.Id, node6.Id], props)

        if (self.PostBoundingBox):
            self.gid_io.WriteMesh(bounding_box_model_part.GetCommunicator().LocalMesh())

    def ComputeAndPrintSeaSurface(self, spheres_model_part, rigid_face_model_part):

        sea_surface_model_part = ModelPart("SeaSurfacePart") # Creation of sea surface model part

        max_node_Id        = ParticleCreatorDestructor().FindMaxNodeIdInModelPart(spheres_model_part)
        max_FEM_node_Id    = ParticleCreatorDestructor().FindMaxNodeIdInModelPart(rigid_face_model_part)
        max_element_Id     = ParticleCreatorDestructor().FindMaxElementIdInModelPart(spheres_model_part)
        max_FEM_element_Id = ParticleCreatorDestructor().FindMaxElementIdInModelPart(rigid_face_model_part)

        if (max_FEM_node_Id > max_node_Id):
            max_node_Id = max_FEM_node_Id

        if (max_FEM_element_Id > max_element_Id):
            max_element_Id = max_FEM_element_Id

        node1 = sea_surface_model_part.CreateNewNode(max_node_Id +  9, self.SeaSurfaceX1, self.SeaSurfaceY1, 0.0) # Z = 0.0 as sea level. We will always assume this value
        node2 = sea_surface_model_part.CreateNewNode(max_node_Id + 10, self.SeaSurfaceX2, self.SeaSurfaceY2, 0.0)
        node3 = sea_surface_model_part.CreateNewNode(max_node_Id + 11, self.SeaSurfaceX3, self.SeaSurfaceY3, 0.0)
        node4 = sea_surface_model_part.CreateNewNode(max_node_Id + 12, self.SeaSurfaceX4, self.SeaSurfaceY4, 0.0)

        ''' Properties colours: 0 -> grey,        1 -> dark blue, 2 -> pink,       3 -> light blue,       4 -> dark red,    5 -> light green
                                6 -> light brown, 7 -> red-brown, 8 -> dark brown, 9 -> dark green/blue, 10 -> dark purple'''

        # Sea Surface Element, consisting in a quadrilateral. Property 3 corresponds to a light blue for water
        sea_surface_model_part.CreateNewCondition("RigidFace3D4N", max_element_Id + 1, [node1.Id, node2.Id, node3.Id, node4.Id], Properties(3))

        self.gid_io.WriteMesh(sea_surface_model_part.GetCommunicator().LocalMesh())

    def ComputeAndPrintDEMFEMSearchBinBoundingBox(self, spheres_model_part, rigid_face_model_part, dem_fem_search):

        bounding_box_model_part = ModelPart("BoundingBoxPart")

        max_node_Id        = ParticleCreatorDestructor().FindMaxNodeIdInModelPart(spheres_model_part)
        max_FEM_node_Id    = ParticleCreatorDestructor().FindMaxNodeIdInModelPart(rigid_face_model_part)
        max_element_Id     = ParticleCreatorDestructor().FindMaxElementIdInModelPart(spheres_model_part)
        max_FEM_element_Id = ParticleCreatorDestructor().FindMaxElementIdInModelPart(rigid_face_model_part)

        if (max_FEM_node_Id > max_node_Id):
            max_node_Id = max_FEM_node_Id

        if (max_FEM_element_Id > max_element_Id):
            max_element_Id = max_FEM_element_Id

        BBMaxX = dem_fem_search.GetBBHighPoint()[0]
        BBMaxY = dem_fem_search.GetBBHighPoint()[1]
        BBMaxZ = dem_fem_search.GetBBHighPoint()[2]
        BBMinX = dem_fem_search.GetBBLowPoint()[0]
        BBMinY = dem_fem_search.GetBBLowPoint()[1]
        BBMinZ = dem_fem_search.GetBBLowPoint()[2]

        DX = (BBMaxX-BBMinX)
        DY = (BBMaxY-BBMinY)
        DZ = (BBMaxZ-BBMinZ)

        #The cases with 0 thickness in one direction, a 10% of the shortest other two is given to the 0-thickness direction.
        if (DX == 0):
          height = min(DY,DZ)
          BBMinX = BBMinX - 0.05*height
          BBMaxX = BBMaxX + 0.05*height
        if (DY == 0):
          height = min(DX,DZ)
          BBMinY = BBMinY - 0.05*height
          BBMaxY = BBMaxY + 0.05*height
        if (DZ == 0):
          height = min(DX,DY)
          BBMinZ = BBMinZ - 0.05*height
          BBMaxZ = BBMaxZ + 0.05*height

        volume = DX*DY*DZ

        if (abs(volume) > 1e21) :
          BBMaxX = 0.0
          BBMaxY = 0.0
          BBMaxZ = 0.0
          BBMinX = 0.0
          BBMinY = 0.0
          BBMinZ = 0.0

        # BB Nodes:
        node1 = bounding_box_model_part.CreateNewNode(max_node_Id + 1, BBMinX, BBMinY, BBMinZ)
        node2 = bounding_box_model_part.CreateNewNode(max_node_Id + 2, BBMaxX, BBMinY, BBMinZ)
        node3 = bounding_box_model_part.CreateNewNode(max_node_Id + 3, BBMaxX, BBMaxY, BBMinZ)
        node4 = bounding_box_model_part.CreateNewNode(max_node_Id + 4, BBMinX, BBMaxY, BBMinZ)
        node5 = bounding_box_model_part.CreateNewNode(max_node_Id + 5, BBMinX, BBMinY, BBMaxZ)
        node6 = bounding_box_model_part.CreateNewNode(max_node_Id + 6, BBMaxX, BBMinY, BBMaxZ)
        node7 = bounding_box_model_part.CreateNewNode(max_node_Id + 7, BBMaxX, BBMaxY, BBMaxZ)
        node8 = bounding_box_model_part.CreateNewNode(max_node_Id + 8, BBMinX, BBMaxY, BBMaxZ)

        # BB Elements:
        props = Properties(10000)
        bounding_box_model_part.CreateNewCondition("RigidEdge3D", max_element_Id +  1, [node1.Id, node4.Id], props)
        bounding_box_model_part.CreateNewCondition("RigidEdge3D", max_element_Id +  2, [node4.Id, node8.Id], props)
        bounding_box_model_part.CreateNewCondition("RigidEdge3D", max_element_Id +  3, [node8.Id, node5.Id], props)
        bounding_box_model_part.CreateNewCondition("RigidEdge3D", max_element_Id +  4, [node5.Id, node1.Id], props)
        bounding_box_model_part.CreateNewCondition("RigidEdge3D", max_element_Id +  5, [node1.Id, node2.Id], props)
        bounding_box_model_part.CreateNewCondition("RigidEdge3D", max_element_Id +  6, [node3.Id, node4.Id], props)
        bounding_box_model_part.CreateNewCondition("RigidEdge3D", max_element_Id +  7, [node7.Id, node8.Id], props)
        bounding_box_model_part.CreateNewCondition("RigidEdge3D", max_element_Id +  8, [node5.Id, node6.Id], props)
        bounding_box_model_part.CreateNewCondition("RigidEdge3D", max_element_Id +  9, [node6.Id, node2.Id], props)
        bounding_box_model_part.CreateNewCondition("RigidEdge3D", max_element_Id + 10, [node2.Id, node3.Id], props)
        bounding_box_model_part.CreateNewCondition("RigidEdge3D", max_element_Id + 11, [node3.Id, node7.Id], props)
        bounding_box_model_part.CreateNewCondition("RigidEdge3D", max_element_Id + 12, [node7.Id, node6.Id], props)

        #self.gid_io.WriteMesh(bounding_box_model_part.GetCommunicator().LocalMesh()) #BOUNDING BOX IMPLEMENTATION


class ParallelUtils(object):

    def __init__(self):
        pass

    def Repart(self, spheres_model_part):
        pass

    def CalculateModelNewIds(self, spheres_model_part):
        pass

    def PerformInitialPartition(self, model_part):
        pass
    
    def SetCommunicator(self, spheres_model_part, model_part_io_spheres, spheres_mp_filename):
        MPICommSetup = 0
        return [model_part_io_spheres, spheres_model_part, MPICommSetup]

    def GetSearchStrategy(self, solver, model_part):
        return solver.search_strategy
