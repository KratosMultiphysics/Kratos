#  Kratos Multi-Physics - DEM Application
#
#  License:       BSD License
#                 Kratos default license: kratos/license.txt
#
#  Main authors:  Miguel Angel Celigueta (maceli@cimne.upc.edu)
#                 Chengshun Shang (cshang@cimne.upc.edu)
#

import os
from pyevtk import hl
from pyevtk import vtk
import numpy as np
import weakref
import shutil
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *

def GetBoolParameterIfItExists(set_of_parameters, parameter_key):
    if parameter_key in set_of_parameters.keys():
        return set_of_parameters[parameter_key].GetBool()
    return False

class VtkOutput():
    def __init__(self, main_path, problem_name, spheres_model_part, contact_model_part, rigid_face_model_part, DEM_parameters):
        self.problem_name = problem_name
        self.DEM_parameters = DEM_parameters

        # Reading Post options from DEM_parameters
        self.PostDisplacement = self.DEM_parameters["PostDisplacement"].GetBool()
        self.PostVelocity = self.DEM_parameters["PostVelocity"].GetBool()
        self.PostTotalForces = self.DEM_parameters["PostTotalForces"].GetBool()
        self.PostNonDimensionalVolumeWear = self.DEM_parameters["PostNonDimensionalVolumeWear"].GetBool()
        self.PostAppliedForces = self.DEM_parameters["PostAppliedForces"].GetBool()
        self.PostDampForces = self.DEM_parameters["PostDampForces"].GetBool()
        self.PostRadius = self.DEM_parameters["PostRadius"].GetBool()
        self.PostGroupId = GetBoolParameterIfItExists(self.DEM_parameters, "PostGroupId")
        self.PostExportId = self.DEM_parameters["PostExportId"].GetBool()
        self.PostSkinSphere = GetBoolParameterIfItExists(self.DEM_parameters, "PostSkinSphere")
        self.PostGluedSphere = GetBoolParameterIfItExists(self.DEM_parameters, "PostGluedSphere")
        self.PostAngularVelocity = self.DEM_parameters["PostAngularVelocity"].GetBool()
        self.PostParticleMoment = self.DEM_parameters["PostParticleMoment"].GetBool()
        self.PostEulerAngles = self.DEM_parameters["PostEulerAngles"].GetBool()
        self.PostRollingResistanceMoment = self.DEM_parameters["PostRollingResistanceMoment"].GetBool()
        self.PostNeighbourSize = GetBoolParameterIfItExists(self.DEM_parameters, "PostNeighbourSize")
        self.PostDamageRatio = GetBoolParameterIfItExists(self.DEM_parameters, "PostDamageRatio")
        
        # for contacts
        self.PostLocalContactForce = GetBoolParameterIfItExists(self.DEM_parameters, "PostLocalContactForce")
        self.PostFailureCriterionState = GetBoolParameterIfItExists(self.DEM_parameters, "PostFailureCriterionState")
        self.PostContactFailureId = GetBoolParameterIfItExists(self.DEM_parameters, "PostContactFailureId")
        self.PostContactTau = GetBoolParameterIfItExists(self.DEM_parameters, "PostContactTau")
        self.PostContactSigma = GetBoolParameterIfItExists(self.DEM_parameters, "PostContactSigma")
        self.PostContactRadius = GetBoolParameterIfItExists(self.DEM_parameters, "PostContactRadius")
        
        #TODO: those variables need to be added
        ''' 
        self.PostElasticForces = self.DEM_parameters["PostElasticForces"].GetBool()
        self.PostContactForces = self.DEM_parameters["PostContactForces"].GetBool()
        self.PostRigidElementForces = self.DEM_parameters["PostRigidElementForces"].GetBool()
        self.PostPressure = self.DEM_parameters["PostPressure"].GetBool()
        self.PostTangentialElasticForces = self.DEM_parameters["PostTangentialElasticForces"].GetBool()
        self.PostShearStress = self.DEM_parameters["PostShearStress"].GetBool()
        self.PostNodalArea = self.DEM_parameters["PostNodalArea"].GetBool()
        self.PostNormalImpactVelocity = GetBoolParameterIfItExists(self.DEM_parameters, "PostNormalImpactVelocity")
        self.PostTangentialImpactVelocity = GetBoolParameterIfItExists(self.DEM_parameters, "PostTangentialImpactVelocity")
        self.PostControlModule = GetBoolParameterIfItExists(self.DEM_parameters, "PostControlModule")
        self.VelTrapGraphExportFreq = self.DEM_parameters["VelTrapGraphExportFreq"].GetDouble()
        self.PostDeltaDisplacement = GetBoolParameterIfItExists(self.DEM_parameters, "PostDeltaDisplacement")
        self.PostCharacteristicLength = GetBoolParameterIfItExists(self.DEM_parameters, "PostCharacteristicLength")
        self.PostBoundingBox = GetBoolParameterIfItExists(self.DEM_parameters, "PostBoundingBox")
        '''
        self.spheres_model_part = weakref.proxy(spheres_model_part)
        self.contact_model_part = weakref.proxy(contact_model_part)
        self.rigid_face_model_part = weakref.proxy(rigid_face_model_part)

        self.vtk_post_path_directory = os.path.join(main_path, problem_name + "_Post_VTK_Files")
        shutil.rmtree(self.vtk_post_path_directory, ignore_errors=True)
        os.makedirs(str(self.vtk_post_path_directory))

        self.counter = 0

    def ConvertParticlesToNumpyArrays(self):
        number_of_nodes = self.spheres_model_part.NumberOfNodes(0)

        if number_of_nodes == 0:
            return False

        self.particles_X = np.empty(number_of_nodes)
        self.particles_Y = np.empty(number_of_nodes)
        self.particles_Z = np.empty(number_of_nodes)
        self.particles_material = np.empty(number_of_nodes)

        if self.PostRadius:
            self.particles_R = np.empty(number_of_nodes)

        if self.PostVelocity:
            self.velocities_X = np.empty(number_of_nodes)
            self.velocities_Y = np.empty(number_of_nodes)
            self.velocities_Z = np.empty(number_of_nodes)
        
        if self.PostDisplacement:
            self.displacement_X = np.empty(number_of_nodes)
            self.displacement_Y = np.empty(number_of_nodes)
            self.displacement_Z = np.empty(number_of_nodes)

        if self.PostTotalForces:
            self.total_force_X = np.empty(number_of_nodes)
            self.total_force_Y = np.empty(number_of_nodes)
            self.total_force_Z = np.empty(number_of_nodes)

        if self.PostNonDimensionalVolumeWear:
            self.non_dimensional_volume_wear = np.empty(number_of_nodes)
        
        if self.PostAppliedForces:
            self.applied_force_X = np.empty(number_of_nodes)
            self.applied_force_Y = np.empty(number_of_nodes)
            self.applied_force_Z = np.empty(number_of_nodes)

        if self.PostDampForces:
            self.damp_force_X = np.empty(number_of_nodes)
            self.damp_force_Y = np.empty(number_of_nodes)
            self.damp_force_Z = np.empty(number_of_nodes)

        if self.PostGroupId:
            self.group_id = np.empty(number_of_nodes)

        if self.PostExportId:
            self.export_id = np.empty(number_of_nodes)

        if self.PostSkinSphere:
            self.skin_sphere = np.empty(number_of_nodes)

        if self.PostGluedSphere:
            self.glued_sphere = np.empty(number_of_nodes)
        
        if self.PostAngularVelocity:
            self.angular_velocity_X = np.empty(number_of_nodes)
            self.angular_velocity_Y = np.empty(number_of_nodes)
            self.angular_velocity_Z = np.empty(number_of_nodes)
        
        if self.PostParticleMoment:
            self.particle_moment_X = np.empty(number_of_nodes)
            self.particle_moment_Y = np.empty(number_of_nodes)
            self.particle_moment_Z = np.empty(number_of_nodes)
        
        if self.PostEulerAngles:
            self.euler_angles_X = np.empty(number_of_nodes)
            self.euler_angles_Y = np.empty(number_of_nodes)
            self.euler_angles_Z = np.empty(number_of_nodes)
        
        if self.PostRollingResistanceMoment:
            self.rolling_resistance_moment_X = np.empty(number_of_nodes)
            self.rolling_resistance_moment_Y = np.empty(number_of_nodes)
            self.rolling_resistance_moment_Z = np.empty(number_of_nodes)
        
        if self.PostNeighbourSize:
            self.neighbour_size = np.empty(number_of_nodes)
        
        if self.PostDamageRatio:
            self.damage_ratio = np.empty(number_of_nodes)

        i = 0
        for node in self.spheres_model_part.Nodes:
            self.particles_X[i] = node.X
            self.particles_Y[i] = node.Y
            self.particles_Z[i] = node.Z
            self.particles_material[i] = node.GetSolutionStepValue(PARTICLE_MATERIAL)

            if self.PostRadius:
                self.particles_R[i] = node.GetSolutionStepValue(RADIUS)

            if self.PostVelocity:
                self.velocities_X[i] = node.GetSolutionStepValue(VELOCITY_X)
                self.velocities_Y[i] = node.GetSolutionStepValue(VELOCITY_Y)
                self.velocities_Z[i] = node.GetSolutionStepValue(VELOCITY_Z)
            
            if self.PostDisplacement:
                self.displacement_X[i] = node.GetSolutionStepValue(DISPLACEMENT_X)
                self.displacement_Y[i] = node.GetSolutionStepValue(DISPLACEMENT_Y)
                self.displacement_Z[i] = node.GetSolutionStepValue(DISPLACEMENT_Z)

            if self.PostTotalForces:
                self.total_force_X[i] = node.GetSolutionStepValue(TOTAL_FORCES)[0]
                self.total_force_Y[i] = node.GetSolutionStepValue(TOTAL_FORCES)[1]
                self.total_force_Z[i] = node.GetSolutionStepValue(TOTAL_FORCES)[2]

            if self.PostNonDimensionalVolumeWear:
                self.non_dimensional_volume_wear[i] = node.GetSolutionStepValue(NON_DIMENSIONAL_VOLUME_WEAR)
            
            if self.PostAppliedForces:
                self.applied_force_X[i] = node.GetSolutionStepValue(EXTERNAL_APPLIED_FORCE)[0]
                self.applied_force_Y[i] = node.GetSolutionStepValue(EXTERNAL_APPLIED_FORCE)[1]
                self.applied_force_Z[i] = node.GetSolutionStepValue(EXTERNAL_APPLIED_FORCE)[2]

            if self.PostDampForces:
                self.damp_force_X[i] = node.GetSolutionStepValue(DAMP_FORCES)[0]
                self.damp_force_Y[i] = node.GetSolutionStepValue(DAMP_FORCES)[1]
                self.damp_force_Z[i] = node.GetSolutionStepValue(DAMP_FORCES)[2]

            if self.PostGroupId:
                self.group_id[i] = node.GetSolutionStepValue(GROUP_ID)

            if self.PostExportId:
                self.export_id[i] = node.GetSolutionStepValue(EXPORT_ID)

            if self.PostSkinSphere:
                self.skin_sphere[i] = node.GetSolutionStepValue(SKIN_SPHERE)

            if self.PostGluedSphere:
                self.glued_sphere[i] = node.GetSolutionStepValue(IS_STICKY)
            
            if self.PostAngularVelocity:
                self.angular_velocity_X[i] = node.GetSolutionStepValue(ANGULAR_VELOCITY)[0]
                self.angular_velocity_Y[i] = node.GetSolutionStepValue(ANGULAR_VELOCITY)[1]
                self.angular_velocity_Z[i] = node.GetSolutionStepValue(ANGULAR_VELOCITY)[2]
            
            if self.PostParticleMoment:
                self.particle_moment_X[i] = node.GetSolutionStepValue(PARTICLE_MOMENT)[0]
                self.particle_moment_Y[i] = node.GetSolutionStepValue(PARTICLE_MOMENT)[1]
                self.particle_moment_Z[i] = node.GetSolutionStepValue(PARTICLE_MOMENT)[2]
            
            if self.PostEulerAngles:
                self.euler_angles_X[i] = node.GetSolutionStepValue(EULER_ANGLES)[0]
                self.euler_angles_Y[i] = node.GetSolutionStepValue(EULER_ANGLES)[1]
                self.euler_angles_Z[i] = node.GetSolutionStepValue(EULER_ANGLES)[2]
            
            if self.PostRollingResistanceMoment:
                self.rolling_resistance_moment_X[i] = node.GetSolutionStepValue(ROLLING_RESISTANCE_MOMENT)[0]
                self.rolling_resistance_moment_Y[i] = node.GetSolutionStepValue(ROLLING_RESISTANCE_MOMENT)[1]
                self.rolling_resistance_moment_Z[i] = node.GetSolutionStepValue(ROLLING_RESISTANCE_MOMENT)[2]
            
            if self.PostNeighbourSize:
                self.neighbour_size[i] = node.GetSolutionStepValue(EXTERNAL_APPLIED_FORCE)
            
            if self.PostDamageRatio:
                self.damage_ratio[i] = node.GetSolutionStepValue(DAMAGE_RATIO)

            i += 1

        return True

    def ConvertContactsToNumpyArrays(self):
        
        max_radius = 0.0
        for node in self.spheres_model_part.Nodes:
            r = node.GetSolutionStepValue(RADIUS)
            if r > max_radius:
                max_radius = r
        
        number_of_elements = 0

        for element in self.contact_model_part.Elements:
            
            X1 = element.GetNode(0).X
            X2 = element.GetNode(1).X
            Y1 = element.GetNode(0).Y
            Y2 = element.GetNode(1).Y
            Z1 = element.GetNode(0).Z
            Z2 = element.GetNode(1).Z

            element_length = ((X2-X1)**2 + (Y2-Y1)**2 + (Z2-Z1)**2)**0.5

            if element_length < 3 * max_radius:
                number_of_elements += 1
        
        if number_of_elements == 0:
            return False
        
        #number_of_elements = self.contact_model_part.NumberOfElements(0)
        number_of_nodes = number_of_elements * 2

        self.start_and_end_points_X = np.empty(number_of_nodes)
        self.start_and_end_points_Y = np.empty(number_of_nodes)
        self.start_and_end_points_Z = np.empty(number_of_nodes)

        if self.PostLocalContactForce:
            self.local_contact_force_X = np.empty(number_of_elements)
            self.local_contact_force_Y = np.empty(number_of_elements)
            self.local_contact_force_Z = np.empty(number_of_elements)
            self.local_contact_force_X_point = np.empty(number_of_nodes)
            self.local_contact_force_Y_point = np.empty(number_of_nodes)
            self.local_contact_force_Z_point = np.empty(number_of_nodes)

            self.local_contact_force_mag = np.empty(number_of_elements)
            self.local_contact_force_mag_point = np.empty(number_of_nodes)

        if self.PostFailureCriterionState:
            self.failure_criterion_state = np.empty(number_of_elements)

        if self.PostContactFailureId:
            self.contact_failure_id = np.empty(number_of_elements)

        if self.PostContactTau:
            self.contact_tau = np.empty(number_of_elements)
            self.contact_tau_point = np.empty(number_of_nodes)

        if self.PostContactSigma:
            self.contact_sigma = np.empty(number_of_elements)
            self.contact_sigma_point = np.empty(number_of_nodes)

        if self.PostContactRadius:
            self.contact_radius = np.empty(number_of_elements)
            self.contact_radius_point = np.empty(number_of_nodes)

        i_point= 0
        i_contact = 0
        for element in self.contact_model_part.Elements:
            
            X1 = element.GetNode(0).X
            X2 = element.GetNode(1).X
            Y1 = element.GetNode(0).Y
            Y2 = element.GetNode(1).Y
            Z1 = element.GetNode(0).Z
            Z2 = element.GetNode(1).Z

            element_length = ((X2-X1)**2 + (Y2-Y1)**2 + (Z2-Z1)**2)**0.5

            if element_length < 3 * max_radius:

                self.start_and_end_points_X[i_point] = X1
                self.start_and_end_points_X[i_point+1] = X2
                self.start_and_end_points_Y[i_point] = Y1
                self.start_and_end_points_Y[i_point+1] = Y2
                self.start_and_end_points_Z[i_point] = Z1
                self.start_and_end_points_Z[i_point+1] = Z2

                if self.PostLocalContactForce:
                    self.local_contact_force_X[i_contact] = element.GetValue(LOCAL_CONTACT_FORCE)[0]
                    self.local_contact_force_Y[i_contact] = element.GetValue(LOCAL_CONTACT_FORCE)[1]
                    self.local_contact_force_Z[i_contact] = element.GetValue(LOCAL_CONTACT_FORCE)[2]
                    self.local_contact_force_X_point[i_point] = element.GetValue(LOCAL_CONTACT_FORCE)[0]
                    self.local_contact_force_X_point[i_point+1] = element.GetValue(LOCAL_CONTACT_FORCE)[0]
                    self.local_contact_force_Y_point[i_point] = element.GetValue(LOCAL_CONTACT_FORCE)[1]
                    self.local_contact_force_Y_point[i_point+1] = element.GetValue(LOCAL_CONTACT_FORCE)[1]
                    self.local_contact_force_Z_point[i_point] = element.GetValue(LOCAL_CONTACT_FORCE)[2]
                    self.local_contact_force_Z_point[i_point+1] = element.GetValue(LOCAL_CONTACT_FORCE)[2]

                    local_contact_force_mag_value = (self.local_contact_force_X[i_contact]**2 + self.local_contact_force_Y[i_contact]**2 + self.local_contact_force_Z[i_contact]**2)**0.5
                    self.local_contact_force_mag[i_contact] = local_contact_force_mag_value
                    self.local_contact_force_mag_point[i_point] = local_contact_force_mag_value
                    self.local_contact_force_mag_point[i_point + 1] = local_contact_force_mag_value

                if self.PostFailureCriterionState:
                    self.failure_criterion_state[i_contact] = element.GetValue(FAILURE_CRITERION_STATE)

                if self.PostContactFailureId:
                    self.contact_failure_id[i_contact] = element.GetValue(CONTACT_FAILURE)

                if self.PostContactTau:
                    self.contact_tau[i_contact] = element.GetValue(CONTACT_TAU)
                    self.contact_tau_point[i_point] = element.GetValue(CONTACT_TAU)
                    self.contact_tau_point[i_point+1] = element.GetValue(CONTACT_TAU)

                if self.PostContactSigma:
                    self.contact_sigma[i_contact] = element.GetValue(CONTACT_SIGMA)
                    self.contact_sigma_point[i_point] = element.GetValue(CONTACT_SIGMA)
                    self.contact_sigma_point[i_point+1] = element.GetValue(CONTACT_SIGMA)

                if self.PostContactRadius:
                    self.contact_radius[i_contact] = element.GetValue(CONTACT_RADIUS)
                    self.contact_radius_point[i_point] = element.GetValue(CONTACT_RADIUS)
                    self.contact_radius_point[i_point+1] = element.GetValue(CONTACT_RADIUS)

                i_point += 2
                i_contact += 1

        return True

    def ConvertWallsToNumpyArrays(self):
        number_of_nodes = self.rigid_face_model_part.NumberOfNodes(0)

        if number_of_nodes == 0:
            return False
        
        self.walls_X = np.empty(number_of_nodes)
        self.walls_Y = np.empty(number_of_nodes)
        self.walls_Z = np.empty(number_of_nodes)

        position_of_each_id = {}
        i = 0
        for node in self.rigid_face_model_part.Nodes:
            self.walls_X[i] = node.X
            self.walls_Y[i] = node.Y
            self.walls_Z[i] = node.Z
            position_of_each_id[node.Id] = i
            i += 1

        number_of_conditions = self.rigid_face_model_part.NumberOfConditions(0)
        self.walls_connectivity = np.empty(number_of_conditions * 4) #TODO: here could be problematic if 'condition' is a triangle.
        self.walls_offsets = np.empty(number_of_conditions)
        self.walls_cell_types = np.empty(number_of_conditions)
        i = 0
        j = 0
        for cond in self.rigid_face_model_part.Conditions:
            list_of_vertices = cond.GetNodes()
            number_of_vertices = len(list_of_vertices)
            for k in range(number_of_vertices):
                index_in_array = position_of_each_id[list_of_vertices[k].Id]
                self.walls_connectivity[i] = index_in_array
                self.walls_offsets[j] = i+1
                i += 1

            self.walls_cell_types[j] = vtk.VtkTriangle.tid
            j += 1

        return True

    def WriteResults(self, time):

        #---------------Partciles part---------------------
        if self.ConvertParticlesToNumpyArrays():
            particles_filename = self.problem_name + "_Particles_" + str(self.counter)
            path = os.path.join(self.vtk_post_path_directory, particles_filename)
            
            particles_output_dict = {'material':self.particles_material}
            
            if self.PostRadius:
                particles_output_dict['radius'] = self.particles_R

            if self.PostVelocity:
                particles_output_dict['velocity'] = (self.velocities_X, self.velocities_Y, self.velocities_Z)
            
            if self.PostDisplacement:
                particles_output_dict['displacement'] = (self.displacement_X, self.displacement_Y, self.displacement_Z)

            if self.PostTotalForces:
                particles_output_dict['total_forces'] = (self.total_force_X, self.total_force_Y, self.total_force_Z)

            if self.PostNonDimensionalVolumeWear:
                particles_output_dict['non_dimensional_volume_wear'] = self.non_dimensional_volume_wear
            
            if self.PostAppliedForces:
                particles_output_dict['applied_forces'] = (self.applied_force_X, self.applied_force_Y, self.applied_force_Z)

            if self.PostDampForces:
                particles_output_dict['damp_forces'] = (self.damp_force_X, self.damp_force_Y, self.damp_force_Z)

            if self.PostGroupId:
                particles_output_dict['group_id'] = self.group_id

            if self.PostExportId:
                particles_output_dict['group_id'] = self.export_id

            if self.PostSkinSphere:
                particles_output_dict['skin_sphere'] = self.skin_sphere
            
            if self.PostAngularVelocity:
                particles_output_dict['angular_velocity'] = (self.angular_velocity_X, self.angular_velocity_Y, self.angular_velocity_Z)
            
            if self.PostParticleMoment:
                particles_output_dict['particle_moment'] = (self.particle_moment_X, self.particle_moment_Y, self.particle_moment_Z)
            
            if self.PostEulerAngles:
                particles_output_dict['euler_angles'] = (self.euler_angles_X, self.euler_angles_Y, self.euler_angles_Z)
            
            if self.PostRollingResistanceMoment:
                particles_output_dict['rolling_resistance_moment'] = (self.rolling_resistance_moment_X, self.rolling_resistance_moment_Y, self.rolling_resistance_moment_Z)
            
            if self.PostNeighbourSize:
                particles_output_dict['neighbour_size'] = self.neighbour_size
            
            if self.PostDamageRatio:
                particles_output_dict['damage_ratio'] = self.damage_ratio

            hl.pointsToVTK(path, self.particles_X, self.particles_Y, self.particles_Z, particles_output_dict)

        #------------------------------Contacts part---------------------------
        # Tips:
        #   To display force chain with tube of different radius,
        #       1. Import the 'Contacts'.vtu file to Paraview;
        #       2. Apply 'Extract Surface' filter;
        #       3. Apply 'Tube' filter, and select 'Vary radius'->'By Ansolute Scalar', and set the 'Radius factor' as 1.
        #----------------------------------------------------------------------
        if self.ConvertContactsToNumpyArrays():
            contacts_filename = self.problem_name + "_Contacts_" + str(self.counter)
            path = os.path.join(self.vtk_post_path_directory, contacts_filename)

            contacts_output_dict_cell = {}
            contacts_output_dict_point = {}

            if self.PostLocalContactForce:
                contacts_output_dict_cell['local_contact_force'] = (self.local_contact_force_X, self.local_contact_force_Y, self.local_contact_force_Z) 
                contacts_output_dict_point['local_contact_force'] = (self.local_contact_force_X_point, self.local_contact_force_Y_point, self.local_contact_force_Z_point)

                contacts_output_dict_cell['local_contact_force_mag'] = self.local_contact_force_mag
                contacts_output_dict_point['local_contact_force_mag'] = self.local_contact_force_mag_point

            if self.PostFailureCriterionState:
                contacts_output_dict_cell['failure_criterion_state'] = self.failure_criterion_state

            if self.PostContactFailureId:
                contacts_output_dict_cell['contact_failure_id'] = self.contact_failure_id

            if self.PostContactTau:
                contacts_output_dict_cell['contact_tau'] = self.contact_tau
                contacts_output_dict_point['contact_tau'] = self.contact_tau_point

            if self.PostContactSigma:
                contacts_output_dict_cell['contact_sigma'] = self.contact_sigma
                contacts_output_dict_point['contact_sigma'] = self.contact_sigma_point

            if self.PostContactRadius:
                contacts_output_dict_cell['contact_radius'] = self.contact_radius
                contacts_output_dict_point['contact_radius'] = self.contact_radius_point

            hl.linesToVTK(path, self.start_and_end_points_X, self.start_and_end_points_Y, self.start_and_end_points_Z, contacts_output_dict_cell, contacts_output_dict_point)
        
        #-------------------------------Walls part------------------------------
        if self.ConvertWallsToNumpyArrays():
            walls_filename = self.problem_name + "_Walls_" + str(self.counter)
            path = os.path.join(self.vtk_post_path_directory, walls_filename)
            hl.unstructuredGridToVTK(path, self.walls_X, self.walls_Y, self.walls_Z, self.walls_connectivity, self.walls_offsets, self.walls_cell_types, cellData=None, pointData=None)

        self.counter += 1

