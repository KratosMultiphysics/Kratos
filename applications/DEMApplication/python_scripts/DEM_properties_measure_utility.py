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
        self.SphericElementGlobalPhysicsCalculator = SphericElementGlobalPhysicsCalculator(self.spheres_model_part)
        self.ContactElementGlobalPhysicsCalculator = ContactElementGlobalPhysicsCalculator()

    def MeasureSphereForGettingPackingProperties(self, radius, center_x, center_y, center_z, type, domain_size=[1,1,1]):        
        '''
        This is a function to establish a sphere to measure local packing properties
        The type could be "porosity", "averaged_coordination_number", "fabric_tensor", "stress_tensor" or "strain" 
        This funtion is only valid for 3D model now
        '''
        if type == "porosity":

            measured_porosity = self.SphericElementGlobalPhysicsCalculator.CalculatePorosityWithinSphere(self.spheres_model_part, radius, [center_x, center_y, center_z])
            return measured_porosity
        
        if type == "averaged_coordination_number":
            
            if self.DEM_parameters["ContactMeshOption"].GetBool():
                averaged_coordination_number = self.ContactElementGlobalPhysicsCalculator.CalculateAveragedCoordinationNumberWithinSphere(self.spheres_model_part, self.contact_model_part, radius, [center_x, center_y, center_z])
                return averaged_coordination_number
            else:
                raise Exception('The \"ContactMeshOption\" in the [ProjectParametersDEM.json] should be [True].')
        
        if type == "fabric_tensor":

            if self.DEM_parameters["ContactMeshOption"].GetBool():
                measured_fabric_tensor = self.ContactElementGlobalPhysicsCalculator.CalculateFabricTensorWithinSphere(self.contact_model_part, radius, [center_x, center_y, center_z])
                measured_fabric_tensor = np.array(measured_fabric_tensor)
                deviatoric_tensor = 4 * (measured_fabric_tensor - 1/3 * np.eye(3)) 
                second_invariant_of_deviatoric_tensor = (0.5 * np.sum(deviatoric_tensor * deviatoric_tensor))**0.5
                eigenvalues, eigenvectors = np.linalg.eig(measured_fabric_tensor)
                return eigenvalues, second_invariant_of_deviatoric_tensor, measured_fabric_tensor
            else:
                raise Exception('The \"ContactMeshOption\" in the [ProjectParametersDEM.json] should be [True].')
            
        if type == "stress_tensor":
            
            if self.DEM_parameters["PostStressStrainOption"].GetBool() and self.DEM_parameters["ContactMeshOption"].GetBool():
                measure_sphere_volume = 4.0 / 3.0 * math.pi * radius * radius * radius
                total_stress_tensor = self.ContactElementGlobalPhysicsCalculator.CalculateTotalStressTensorWithinSphere(self.contact_model_part, radius, [center_x, center_y, center_z])
                averaged_stress_tensor = np.array(total_stress_tensor) / measure_sphere_volume
                return averaged_stress_tensor
            else:
                raise Exception('The \"PostStressStrainOption\" and \"ContactMeshOption\" in the [ProjectParametersDEM.json] should be [True].')
            
        if type == "unbalanced_force":

            if self.DEM_parameters["ContactMeshOption"].GetBool():
                measured_unbalanced_force = self.ContactElementGlobalPhysicsCalculator.CalculateUnbalancedForceWithinSphere(self.spheres_model_part, self.contact_model_part, radius, [center_x, center_y, center_z])
                return measured_unbalanced_force
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

                particle_number_inside = self.SphericElementGlobalPhysicsCalculator.CalculateSumOfParticlesWithinSphere(self.spheres_model_part, radius, [center_x, center_y, center_z])

                if total_contact_number:
                    measured_non_homogenized_conductivity_tensor = total_tensor/measure_sphere_volume
                else:
                    measured_non_homogenized_conductivity_tensor = np.empty((3, 3))

                conductivity_tensor_trace = (measured_non_homogenized_conductivity_tensor[0][0] + measured_non_homogenized_conductivity_tensor[1][1] + measured_non_homogenized_conductivity_tensor[2][2])/3

                return particle_number_inside, [measured_non_homogenized_conductivity_tensor[0][0],measured_non_homogenized_conductivity_tensor[1][1],measured_non_homogenized_conductivity_tensor[2][2]], conductivity_tensor_trace, angles_xy, angles_xz, angles_yz
                
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

        if type == "strain":
            pass
    
    def MeasureSphereForGettingGlobalStressTensor(self, Lx, Ly, Lz):

        if self.DEM_parameters["PostStressStrainOption"].GetBool() and self.DEM_parameters["ContactMeshOption"].GetBool():
            bounding_box_volume = Lx * Ly * Lz
            total_tensor = self.ContactElementGlobalPhysicsCalculator.CalculateTotalStressTensor(self.contact_model_part, Lx, Ly, Lz)
            averaged_total_tensor = np.array(total_tensor) / bounding_box_volume
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