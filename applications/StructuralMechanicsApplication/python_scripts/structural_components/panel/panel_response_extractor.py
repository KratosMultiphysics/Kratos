from KratosMultiphysics.StructuralMechanicsApplication.structural_components.panel.panel_data import PanelCompositeResponse
from KratosMultiphysics.StructuralMechanicsApplication.structural_components.panel.panel_helper import NumberOfPlies 
from KratosMultiphysics.StructuralMechanicsApplication.structural_components.panel.panel_data import PanelCompositeResponse, PanelLoads
import KratosMultiphysics.StructuralMechanicsApplication as SMA
import numpy as np

class PanelResponseExtractor:

    def ExtractCompositeResponse(self, sub_model_part, material):
               
        stresses = SMA.SHELL_ORTHOTROPIC_STRESS_THROUGH_THICKNESS
        force = SMA.SHELL_FORCE_GLOBAL
        moment = SMA.SHELL_MOMENT_GLOBAL
        process_info = sub_model_part.ProcessInfo      

# extracting the total number of layers          
        total_plies = 0
        for element in sub_model_part.Elements:
            total_plies += NumberOfPlies(element.Properties)     
        
        if total_plies == 0:
            raise RuntimeError("No plies defined in any element of the submodelpart")
        
# initializing empty arrays, double length according to top and bottom values
        sigma_1 = np.zeros(2*total_plies)
        sigma_2 = np.zeros(2*total_plies)
        tau_21  = np.zeros(2*total_plies)

        num_elements = sub_model_part.NumberOfElements()

        force_x = np.zeros(num_elements)  
        force_y = np.zeros(num_elements)
        force_xy  = np.zeros(num_elements)

        moment_x = np.zeros(num_elements)
        moment_y = np.zeros(num_elements)
        moment_xy  = np.zeros(num_elements)


        plies_element = total_plies/num_elements # calculating the plies for each element, assuming each element has the same ply size 
        num_elements = num_elements  
        
# start extracting esponse data        
        idx_element = 0 # index for counting the elements
        idx = 0 # global index for mapping values to their respective rows across all elements
# iterate through elements using a loop         
        for element in sub_model_part.Elements:

    # determine the number of Gaussian points per element
            integration_points = element.GetIntegrationPoints()
            number_integrationpoints = len(integration_points)

            stresses_sum = None
            force_sum = None
            moment_sum = None         

# calculate the average: force
            for G_punkt in range(number_integrationpoints):
                force_gp = element.CalculateOnIntegrationPoints(force, process_info)[G_punkt]  # value at current Gaussian point
                
                if force_sum is None:
                    force_sum = force_gp   # gaussian point 0 as starting value
                else:
                    force_sum = force_sum + force_gp   
            # average
            force_gp_elementaxis = force_sum / number_integrationpoints           

            # extracting forces from tensor (3x3 matrix) 
            force_x[idx_element] = force_gp_elementaxis[0, 0]
            force_y[idx_element] = force_gp_elementaxis[1, 1]
            force_xy[idx_element]  = force_gp_elementaxis[0, 1] 
            
# calculate the average: moment
            for G_punkt in range(number_integrationpoints):
                moment_gp = element.CalculateOnIntegrationPoints(moment, process_info)[G_punkt]  # value at current Gaussian point

                if moment_sum is None:
                    moment_sum = moment_gp   # gaussian point 0 as starting value
                else:
                    moment_sum = moment_sum + moment_gp  
            # average
            moment_gp_elementaxis = moment_sum / number_integrationpoints           

            # extracting forces from tensor (3x3 matrix) 
            moment_x[idx_element] = moment_gp_elementaxis[0, 0]
            moment_y[idx_element] = moment_gp_elementaxis[1, 1]
            moment_xy[idx_element]  = moment_gp_elementaxis[0, 1] 

            idx_element += 1
 
# calculate the average: stresses
            for G_punkt in range(number_integrationpoints):
                stresses_gp = element.CalculateOnIntegrationPoints(stresses, process_info)[G_punkt]  # value at current Gaussian point

                if stresses_sum is None:
                    stresses_sum = stresses_gp   # gaussian point 0 as starting value
                else:
                    stresses_sum = stresses_sum + stresses_gp   
            # average
            stresses_gp_elementaxis = stresses_sum / number_integrationpoints

           # stresses_gp_elementaxis = Kratos.Matrix, transferring the entries to NumPy-Vector
            rows = stresses_gp_elementaxis.Size1()            
            cols = stresses_gp_elementaxis.Size2() 

# iterate through plies using a loop 
            for PlyId in range(rows): 

                # row -> NumPy-Vector
                current_stress_data = np.zeros(cols)
               
                for col in range(cols):
                    current_stress_data[col] = stresses_gp_elementaxis[PlyId, col]

                    current_stress_data_2D = current_stress_data[:3] # shortening to 3 entries because ShellThick returns, for example, 5 entries where the last 2 are irrelevant or unknown and essentially 0

# transformation of the stress into a ply coordinate system
                ply_angle_rad = np.radians(material.ply_angle[idx])

                cos = np.cos(ply_angle_rad)
                sin = np.sin(ply_angle_rad)

                T = np.array([
                    [cos**2, sin**2,  2*sin*cos],
                    [sin**2, cos**2, -2*sin*cos],
                    [-sin*cos, sin*cos, cos**2 - sin**2]
                ])
               
                stress_state_plyaxis = T @ current_stress_data_2D 

# filling the vectors
                sigma_1[idx] = stress_state_plyaxis[0] 
                sigma_2[idx] = stress_state_plyaxis[1]
                tau_21[idx] = stress_state_plyaxis[2]

                idx += 1 

#filling the container
                composite_response = PanelCompositeResponse(sigma_1,
                                                            sigma_2,
                                                            tau_21,

                                                            force_x, 
                                                            force_y,
                                                            force_xy, 

                                                            moment_x,
                                                            moment_y, 
                                                            moment_xy, 

                                                            plies_element,  
                                                            num_elements)


        return composite_response

    def ExtractRunningLoads(self, sub_model_part, geometry) -> PanelLoads:
        """Running loads N, M per element in the component frame (tension positive).

        Args:
            sub_model_part: Shell elements of the component.
            geometry: PanelGeometry with the component axes and dimensions a, b.
        """
        R = geometry.rotation_matrix
        process_info = sub_model_part.ProcessInfo
        elements = list(sub_model_part.Elements)

        N = np.array([self._InPlaneComponents(self._AverageOverIntegrationPoints(element, SMA.SHELL_FORCE_GLOBAL, process_info), R)
                      for element in elements])
        M = np.array([self._InPlaneComponents(self._AverageOverIntegrationPoints(element, SMA.SHELL_MOMENT_GLOBAL, process_info), R)
                      for element in elements])

        return PanelLoads(element_ids=np.array([element.Id for element in elements]),
                          element_areas=np.array([element.GetGeometry().Area() for element in elements]),
                          centroids=self._NormalizedCentroids(sub_model_part, elements, geometry),
                          N=N, M=M)

    @staticmethod
    def _AverageOverIntegrationPoints(element, variable, process_info) -> np.ndarray:
        values = element.CalculateOnIntegrationPoints(variable, process_info)
        return np.mean([[[value[i, j] for j in range(3)] for i in range(3)] for value in values], axis=0)

    @staticmethod
    def _InPlaneComponents(tensor_global: np.ndarray, R: np.ndarray) -> np.ndarray:
        """Rotate a 3x3 tensor from global into the component frame and return (xx, yy, xy)."""
        tensor_local = R @ tensor_global @ R.T
        return np.array([tensor_local[0, 0], tensor_local[1, 1], tensor_local[0, 1]])

    @staticmethod
    def _NormalizedCentroids(sub_model_part, elements, geometry) -> np.ndarray:
        """Element centroids as (xi, eta) in [0, 1], measured from the panel corner in the panel frame."""
        R = geometry.rotation_matrix
        nodes = np.array([[node.X0, node.Y0, node.Z0] for node in sub_model_part.Nodes]) @ R.T
        corner = nodes.min(axis=0)[:2]
        centroids = np.array([np.mean([[node.X0, node.Y0, node.Z0] for node in element.GetGeometry()], axis=0)
                              for element in elements]) @ R.T
        return (centroids[:, :2] - corner) / np.array([geometry.a, geometry.b])
