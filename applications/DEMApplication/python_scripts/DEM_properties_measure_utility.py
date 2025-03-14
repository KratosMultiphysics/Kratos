import math
import numpy as np

from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *

class DEMPropertiesMeasureUtility:

    def __init__(self, DEM_parameters, spheres_model_part, contact_model_part, graphs_path):
        self.spheres_model_part = spheres_model_part
        self.contact_model_part = contact_model_part
        self.DEM_parameters = DEM_parameters
        self.graphs_path = graphs_path
        self.ContactElementGlobalPhysicsCalculator = ContactElementGlobalPhysicsCalculator()

    def MeasureSphereForGettingPackingProperties(self, radius, center_x, center_y, center_z, type, domain_size=[1,1,1]):        
        '''
        This is a function to establish a sphere to measure local packing properties
        The type could be "porosity", "averaged_coordination_number", "fabric_tensor", "stress_tensor" or "strain" 
        This funtion is only valid for 3D model now
        '''
        if type == "porosity":

            measure_sphere_volume = 4.0 / 3.0 * math.pi * radius * radius * radius
            sphere_volume_inside_range = 0.0
            measured_porosity = 0.0

            for node in self.spheres_model_part.Nodes:
                
                r = node.GetSolutionStepValue(RADIUS)
                x = node.X
                y = node.Y
                z = node.Z

                center_to_sphere_distance = ((x - center_x)**2 + (y - center_y)**2 + (z - center_z)**2)**0.5

                if center_to_sphere_distance < (radius - r):

                    sphere_volume_inside_range += 4/3 * math.pi * r * r * r

                elif center_to_sphere_distance <= (radius + r):

                    other_part_d = radius - (radius * radius + center_to_sphere_distance * center_to_sphere_distance - r * r) / (center_to_sphere_distance * 2)

                    my_part_d = r - (r * r + center_to_sphere_distance * center_to_sphere_distance - radius * radius) / (center_to_sphere_distance * 2)
                    
                    cross_volume = math.pi * other_part_d * other_part_d * (radius - 1/3 * other_part_d) + math.pi * my_part_d * my_part_d * (r - 1/3 * my_part_d)
                    
                    sphere_volume_inside_range += cross_volume
            
            measured_porosity = 1.0 - (sphere_volume_inside_range / measure_sphere_volume)

            return measured_porosity
        
        if type == "averaged_coordination_number":
            
            measured_coordination_number = 0
            if self.DEM_parameters["ContactMeshOption"].GetBool():
                
                total_particle_number = 0
                total_contact_number  = 0
                for element in self.contact_model_part.Elements:
            
                    x_0 = element.GetNode(0).X
                    x_1 = element.GetNode(1).X
                    y_0 = element.GetNode(0).Y
                    y_1 = element.GetNode(1).Y
                    z_0 = element.GetNode(0).Z
                    z_1 = element.GetNode(1).Z
                    r_0 = element.GetNode(0).GetSolutionStepValue(RADIUS)
                    r_1 = element.GetNode(1).GetSolutionStepValue(RADIUS)

                    center_to_sphere_distance_0 = ((x_0 - center_x)**2 + (y_0 - center_y)**2 + (z_0 - center_z)**2)**0.5
                    center_to_sphere_distance_1 = ((x_1 - center_x)**2 + (y_1 - center_y)**2 + (z_1 - center_z)**2)**0.5

                    if center_to_sphere_distance_0 < (radius - r_0):
                        total_contact_number += 1

                    if center_to_sphere_distance_1 < (radius - r_1):
                        total_contact_number += 1

                for node in self.spheres_model_part.Nodes:

                    r = node.GetSolutionStepValue(RADIUS)
                    x = node.X
                    y = node.Y
                    z = node.Z

                    center_to_sphere_distance = ((x - center_x)**2 + (y - center_y)**2 + (z - center_z)**2)**0.5

                    if center_to_sphere_distance < (radius - r):
                        total_particle_number += 1    
                
                if total_particle_number:
                    measured_coordination_number = total_contact_number / total_particle_number
                
                return measured_coordination_number

            else:
                raise Exception('The \"ContactMeshOption\" in the [ProjectParametersDEM.json] should be [True].')
        
        if type == "fabric_tensor":

            if self.DEM_parameters["ContactMeshOption"].GetBool():
                
                total_tensor = np.empty((3, 3))
                total_tensor[:] = 0.0
                total_contact_number  = 0

                for element in self.contact_model_part.Elements:
            
                    x_0 = element.GetNode(0).X
                    x_1 = element.GetNode(1).X
                    y_0 = element.GetNode(0).Y
                    y_1 = element.GetNode(1).Y
                    z_0 = element.GetNode(0).Z
                    z_1 = element.GetNode(1).Z
                    r_0 = element.GetNode(0).GetSolutionStepValue(RADIUS)
                    r_1 = element.GetNode(1).GetSolutionStepValue(RADIUS)

                    center_to_sphere_distance_0 = ((x_0 - center_x)**2 + (y_0 - center_y)**2 + (z_0 - center_z)**2)**0.5
                    center_to_sphere_distance_1 = ((x_1 - center_x)**2 + (y_1 - center_y)**2 + (z_1 - center_z)**2)**0.5

                    if center_to_sphere_distance_0 < (radius - r_0) or center_to_sphere_distance_1 < (radius - r_1):

                        vector1 = np.array([x_1 - x_0 , y_1 - y_0, z_1 - z_0])
                        v1_norm = np.linalg.norm(vector1)
                        if v1_norm:
                            vector1_unit = vector1 / v1_norm
                        tensor = np.outer(vector1_unit, vector1_unit)
                        total_tensor += tensor
                        total_contact_number += 1
                
                if total_contact_number:
                    measured_fabric_tensor = total_tensor / total_contact_number
                else:
                    measured_fabric_tensor = np.empty((3, 3))

                deviatoric_tensor = 4 * (measured_fabric_tensor - 1/3 * np.eye(3)) 

                second_invariant_of_deviatoric_tensor = (0.5 * np.sum(deviatoric_tensor * deviatoric_tensor))**0.5

                eigenvalues, eigenvectors = np.linalg.eig(measured_fabric_tensor)
                
                return eigenvalues, second_invariant_of_deviatoric_tensor, measured_fabric_tensor

            else:
                raise Exception('The \"ContactMeshOption\" in the [ProjectParametersDEM.json] should be [True].')
            
        if type == "conductivity_tensor":

            if self.DEM_parameters["ContactMeshOption"].GetBool():
                measure_sphere_volume = 4.0 / 3.0 * math.pi * radius * radius * radius
                total_tensor = np.empty((3, 3))
                total_tensor[:] = 0.0
                total_contact_number = 0
                angles_xy = []
                angles_xz = []
                angles_yz = []

                for element in self.contact_model_part.Elements:

                    x_0 = element.GetNode(0).X
                    x_1 = element.GetNode(1).X
                    y_0 = element.GetNode(0).Y
                    y_1 = element.GetNode(1).Y
                    z_0 = element.GetNode(0).Z
                    z_1 = element.GetNode(1).Z
                    r_0 = element.GetNode(0).GetSolutionStepValue(RADIUS)
                    r_1 = element.GetNode(1).GetSolutionStepValue(RADIUS)

                    center_to_sphere_distance_0 = ((x_0 - center_x)**2 + (y_0 - center_y)**2 + (z_0 - center_z)**2)**0.5
                    center_to_sphere_distance_1 = ((x_1 - center_x)**2 + (y_1 - center_y)**2 + (z_1 - center_z)**2)**0.5

                    if center_to_sphere_distance_0 < (radius - r_0) or center_to_sphere_distance_1 < (radius - r_1):
                        d = ((x_0 - x_1)**2 + (y_0 - y_1)**2 + (z_0 - z_1)**2)**0.5
                        if d <= (r_0 + r_1):
                            a = np.sqrt(4*r_0**2*d**2-(r_0**2+d**2-r_1**2)**2)/(2*d)
                            contact_area = np.pi*a**2
                            vector1 = np.array([x_1 - x_0 , y_1 - y_0, z_1 - z_0])
                            v1_norm = np.linalg.norm(vector1)
                            if v1_norm:
                                vector1_unit = vector1 / v1_norm
                            theta_xz = np.arctan2((x_0 - x_1),(z_0 - z_1))*180/np.pi
                            theta_xy = np.arctan2((x_0 - x_1),(y_0 - y_1))*180/np.pi
                            theta_yz = np.arctan2((y_0 - y_1),(z_0 - z_1))*180/np.pi
                            if theta_yz < 0: theta_yz += 360
                            if theta_xz < 0: theta_xz += 360
                            if theta_xy < 0: theta_xy += 360
                            tensor = contact_area*d*np.outer(vector1_unit, vector1_unit)
                            angles_xy.append(theta_xy)
                            angles_xz.append(theta_xz)
                            angles_yz.append(theta_yz)
                            total_tensor += tensor
                            total_contact_number += 1
                        else:
                            continue

                particle_is_inside = 0
                for node in self.spheres_model_part.Nodes:
                    x = node.X
                    y = node.Y
                    z = node.Z
                    r = node.GetSolutionStepValue(RADIUS)
                    center_to_sphere_distance = ((x - center_x)**2 + (y - center_y)**2 + (z - center_z)**2)**0.5
                    if center_to_sphere_distance < (radius - r):
                        particle_is_inside += 1

                if total_contact_number:
                    measured_non_homogenized_conductivity_tensor = total_tensor/measure_sphere_volume
                else:
                    measured_non_homogenized_conductivity_tensor = np.empty((3, 3))

                conductivity_tensor_trace = (measured_non_homogenized_conductivity_tensor[0][0] + measured_non_homogenized_conductivity_tensor[1][1] + measured_non_homogenized_conductivity_tensor[2][2])/3

                deviatoric_tensor = 0.0

                second_invariant_of_deviatoric_tensor = 0.0

                eigenvalues = 0.0

                return particle_is_inside, [measured_non_homogenized_conductivity_tensor[0][0],measured_non_homogenized_conductivity_tensor[1][1],measured_non_homogenized_conductivity_tensor[2][2]], \
                        conductivity_tensor_trace, angles_xy, angles_xz, angles_yz
                
        if type == "voronoi_input_data":

            particle_id_positions_and_radius = np.empty((0, 5))
            particle_id_positions_and_radius[:] = 0.0
            particle_number_count = 0
            for node in self.spheres_model_part.Nodes:

                r = node.GetSolutionStepValue(RADIUS)
                x = node.X
                y = node.Y
                z = node.Z

                center_to_sphere_distance = ((x - center_x)**2 + (y - center_y)**2 + (z - center_z)**2)**0.5

                if center_to_sphere_distance < (radius - r):

                    particle_number_count += 1
                    this_particle_info = np.array([particle_number_count,x, y, z, r])
                    particle_id_positions_and_radius = np.vstack((particle_id_positions_and_radius, this_particle_info))

            output_file_name = "voronoi_input_data_of_size_" + str(radius) +".txt"
            fmt_list = ['%d', '%.6f', '%.6f', '%.6f', '%.6f']
            np.savetxt(os.path.join(self.graphs_path, output_file_name), particle_id_positions_and_radius, fmt=fmt_list, delimiter='\t', comments='')
        
        if type == "stress_tensor":
            
            if self.DEM_parameters["PostStressStrainOption"].GetBool() and self.DEM_parameters["ContactMeshOption"].GetBool():
                
                measure_sphere_volume = 4.0 / 3.0 * math.pi * radius * radius * radius
                total_tensor        = np.empty((3, 3))
                total_tensor[:]     = 0.0
                stress_tensor_modulus = 0.0

                for element in self.contact_model_part.Elements:
            
                    x_0 = element.GetNode(0).X
                    x_1 = element.GetNode(1).X
                    y_0 = element.GetNode(0).Y
                    y_1 = element.GetNode(1).Y
                    z_0 = element.GetNode(0).Z
                    z_1 = element.GetNode(1).Z
                    r_0 = element.GetNode(0).GetSolutionStepValue(RADIUS)
                    r_1 = element.GetNode(1).GetSolutionStepValue(RADIUS)
                    r   = 0.5 * (r_0 + r_1)

                    center_to_sphere_distance_0 = ((x_0 - center_x)**2 + (y_0 - center_y)**2 + (z_0 - center_z)**2)**0.5
                    center_to_sphere_distance_1 = ((x_1 - center_x)**2 + (y_1 - center_y)**2 + (z_1 - center_z)**2)**0.5

                    if (center_to_sphere_distance_0 < (radius - r)) or (center_to_sphere_distance_1 < (radius - r)):
                        
                        local_contact_force_X = element.GetValue(GLOBAL_CONTACT_FORCE)[0]
                        local_contact_force_Y = element.GetValue(GLOBAL_CONTACT_FORCE)[1]
                        local_contact_force_Z = element.GetValue(GLOBAL_CONTACT_FORCE)[2]
                        contact_force_vector = np.array([local_contact_force_X , local_contact_force_Y, local_contact_force_Z])
                        vector_l = np.array([x_0 - x_1 , y_0 - y_1, z_0 - z_1])
                        tensor = np.outer(contact_force_vector, vector_l)
                        total_tensor += tensor
                
                total_tensor = total_tensor / measure_sphere_volume

                return total_tensor
            
            else:
                
                raise Exception('The \"PostStressStrainOption\" and \"ContactMeshOption\" in the [ProjectParametersDEM.json] should be [True].')
            
        if type == "unbalanced_force":

            total_particle_force_tensor_modulus_square = 0.0
            averaged_total_particle_force_tensor_modulus_square = 0.0
            particle_number_count = 0

            for node in self.spheres_model_part.Nodes:

                r = node.GetSolutionStepValue(RADIUS)
                x = node.X
                y = node.Y
                z = node.Z

                center_to_sphere_distance = ((x - center_x)**2 + (y - center_y)**2 + (z - center_z)**2)**0.5

                if center_to_sphere_distance < (radius - r):
                    total_force_x = node.GetSolutionStepValue(TOTAL_FORCES)[0]
                    total_force_y = node.GetSolutionStepValue(TOTAL_FORCES)[1]
                    total_force_z = node.GetSolutionStepValue(TOTAL_FORCES)[2]
                    total_force_vector = np.array([total_force_x, total_force_y, total_force_z])
                    total_force_vector_modulus = np.linalg.norm(total_force_vector)
                    total_particle_force_tensor_modulus_square += total_force_vector_modulus**2
                    particle_number_count += 1
            
            if particle_number_count:
                averaged_total_particle_force_tensor_modulus_square = total_particle_force_tensor_modulus_square / particle_number_count
            
            if self.DEM_parameters["ContactMeshOption"].GetBool():
                
                total_contact_force_tensor_modulus_square = 0.0
                averaged_contact_force_modulus_square = 0.0
                total_contact_number  = 0

                for element in self.contact_model_part.Elements:
            
                    x_0 = element.GetNode(0).X
                    x_1 = element.GetNode(1).X
                    y_0 = element.GetNode(0).Y
                    y_1 = element.GetNode(1).Y
                    z_0 = element.GetNode(0).Z
                    z_1 = element.GetNode(1).Z
                    r_0 = element.GetNode(0).GetSolutionStepValue(RADIUS)
                    r_1 = element.GetNode(1).GetSolutionStepValue(RADIUS)
                    r   = 0.5 * (r_0 + r_1)

                    center_to_sphere_distance_0 = ((x_0 - center_x)**2 + (y_0 - center_y)**2 + (z_0 - center_z)**2)**0.5
                    center_to_sphere_distance_1 = ((x_1 - center_x)**2 + (y_1 - center_y)**2 + (z_1 - center_z)**2)**0.5

                    if (center_to_sphere_distance_0 < (radius - r)) or (center_to_sphere_distance_1 < (radius - r)):
                        
                        local_contact_force_X = element.GetValue(LOCAL_CONTACT_FORCE)[0]
                        local_contact_force_Y = element.GetValue(LOCAL_CONTACT_FORCE)[1]
                        local_contact_force_Z = element.GetValue(LOCAL_CONTACT_FORCE)[2]
                        contact_force_vector = np.array([local_contact_force_X , local_contact_force_Y, local_contact_force_Z])
                        contact_force_vector_modulus = np.linalg.norm(contact_force_vector)
                        total_contact_force_tensor_modulus_square += contact_force_vector_modulus**2
                        total_contact_number += 1
                
                if total_contact_number:
                    averaged_contact_force_modulus_square = total_contact_force_tensor_modulus_square / total_contact_number
                                             
                if averaged_contact_force_modulus_square:
                    return (averaged_total_particle_force_tensor_modulus_square / averaged_contact_force_modulus_square)**0.5
                else:
                    return 0.0
            else:

                raise Exception('The \"ContactMeshOption\" in the [ProjectParametersDEM.json] should be [True].')

        if type == "strain":
            pass
    
    def MeasureSphereForGettingGlobalStressTensor(self, Lx, Ly, Lz):

        if self.DEM_parameters["PostStressStrainOption"].GetBool() and self.DEM_parameters["ContactMeshOption"].GetBool():
            

            bounding_box_volume = Lx * Ly * Lz
            total_tensor = np.zeros((3, 3))
            temp_total_tensor = self.ContactElementGlobalPhysicsCalculator.CalculateTotalStressTensor(self.contact_model_part, Lx, Ly, Lz)
            
            for i in range(3):
                for j in range(3):
                    total_tensor[i][j] = temp_total_tensor[i][j]

            averaged_total_tensor = total_tensor / bounding_box_volume

            return averaged_total_tensor
        
        else:
            
            raise Exception('The \"PostStressStrainOption\" and \"ContactMeshOption\" in the [ProjectParametersDEM.json] should be [True].')
    
    def MeasureSphereForGettingRadialDistributionFunction(self, radius, center_x, center_y, center_z, delta_r, d_mean):
        
        min_reference_particle_to_center_distance = 1e10
        particle_positions = np.empty((0, 3))
        particle_positions[:] = 0.0
        IsTheFirstParticle = True
        TotalParticleNumber = 0
        reference_particle = np.empty((0, 3))
        reference_particle[:] = 0.0
        
        for node in self.spheres_model_part.Nodes:

            r = node.GetSolutionStepValue(RADIUS)
            x = node.X
            y = node.Y
            z = node.Z

            center_to_sphere_distance = ((x - center_x)**2 + (y - center_y)**2 + (z - center_z)**2)**0.5

            if center_to_sphere_distance < radius:

                this_particle_position = np.array([x, y, z])

                if center_to_sphere_distance < min_reference_particle_to_center_distance:
                    min_reference_particle_to_center_distance = center_to_sphere_distance
                    if not IsTheFirstParticle:
                        particle_positions = np.vstack((particle_positions, reference_particle))
                        TotalParticleNumber += 1
                    reference_particle = this_particle_position
                    if IsTheFirstParticle:
                        IsTheFirstParticle = False
                else:
                    particle_positions = np.vstack((particle_positions, this_particle_position))
                    TotalParticleNumber += 1
        
        distances = np.linalg.norm(particle_positions - reference_particle, axis=1)

        max_distance = radius
        num_bins = int(max_distance // delta_r)
        bin_edges = np.linspace(0, max_distance, num_bins + 1)

        hist, _ = np.histogram(distances, bins=bin_edges)

        bin_width = bin_edges[1] - bin_edges[0]
        measure_sphere_volume = 4/3 * math.pi * radius * radius * radius
        rdf = hist / (4 * np.pi * bin_edges[1:]**2 * bin_width * TotalParticleNumber / measure_sphere_volume)

        data_to_save = np.column_stack((bin_edges[1:] / d_mean, rdf))
        output_file_name = "rdf_data_of_size_" + str(radius * 2) +".txt"
        np.savetxt(os.path.join(self.graphs_path, output_file_name), data_to_save, fmt='%.6f', delimiter='\t', comments='')