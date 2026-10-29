import KratosMultiphysics
from KratosMultiphysics.StructuralMechanicsApplication.structural_components.panel.panel_data import PanelCompositeMaterial
from KratosMultiphysics.StructuralMechanicsApplication.structural_components.panel.panel_helper import NumberOfPlies 
from KratosMultiphysics.StructuralMechanicsApplication.structural_components.panel.laminate import Laminate

import KratosMultiphysics.StructuralMechanicsApplication as SMA
import numpy as np

class PanelMaterialExtractor:

    def ExtractCompositeMaterial(self, sub_model_part, metadata):

        #TODO: Implement this function
# check if SubModelPart include elements
        if sub_model_part.NumberOfElements() == 0:
            raise RuntimeError("Panel submodelpart contains no elements")
        
        layers_variable = SMA.SHELL_ORTHOTROPIC_LAYERS

# extracting the total number of layers   
        total_plies = 0 
        for element in sub_model_part.Elements:
            total_plies += NumberOfPlies(element.Properties) 
        if total_plies == 0: 
            raise RuntimeError("No plies defined in any element of the submodelpart")

            

# initializing empty arrays, double length according to top and bottom values
        element_id       = np.zeros(2*total_plies, dtype=int)
        ply_id           = np.zeros(2*total_plies, dtype=int)
        ply_side         = np.zeros(2*total_plies, dtype=int)
        ply_thickness    = np.zeros(total_plies) # No duplicate values are required for “ply_thickness”
        ply_angle        = np.zeros(2*total_plies)

        strength_R_pa_t = np.zeros(2*total_plies)
        strength_R_pa_c = np.zeros(2*total_plies)
        strength_R_tr_t = np.zeros(2*total_plies)
        strength_R_tr_c = np.zeros(2*total_plies)
        strength_R_trpa = np.zeros(2*total_plies)

        inclination_p_trtr_c = np.zeros(2*total_plies)
        inclination_p_trpa_t = np.zeros(2*total_plies)
        inclination_p_trpa_c = np.zeros(2*total_plies)

        youngs_modul_E_pa    = np.zeros(2*total_plies)
        youngs_modul_E_tr    = np.zeros(2*total_plies)
        poissons_ratio_nu_12 = np.zeros(2*total_plies)
        shear_modulus_G_trpa = np.zeros(2*total_plies)

# extracting user inputs
        # inclination factor
        if metadata is not None and metadata.Has("puck_parameters"):
            if metadata["puck_parameters"].Has("p_trtr_c"):
                p_trtr_c = metadata["puck_parameters"]["p_trtr_c"].GetDouble()
            else: 
                p_trtr_c = 0.25

            if metadata["puck_parameters"].Has("p_trpa_t"):
                p_trpa_t = metadata["puck_parameters"]["p_trpa_t"].GetDouble()
            else:
                p_trpa_t = 0.25

            if metadata["puck_parameters"].Has("p_trpa_c"):
                p_trpa_c = metadata["puck_parameters"]["p_trpa_c"].GetDouble()
            else:
                p_trpa_c = 0.25

        else:
            KratosMultiphysics.Logger.PrintInfo("::[Puck Analysis]::", "No puck parameters found in metadata. Using default values.")
            p_trtr_c = 0.25
            p_trpa_t = 0.25
            p_trpa_c = 0.25

        # degradation factor
        if metadata is not None and metadata.Has("degradation_factors"):
            if metadata["degradation_factors"].Has("E_tr_A"):
                E_tr_A = metadata["degradation_factors"]["E_tr_A"].GetDouble()
            else: 
                E_tr_A = 0

            if metadata["degradation_factors"].Has("G_patr_A"):
                G_patr_A = metadata["degradation_factors"]["G_patr_A"].GetDouble()
            else:
                G_patr_A = 0

            if metadata["degradation_factors"].Has("E_tr_B"):
                E_tr_B = metadata["degradation_factors"]["E_tr_B"].GetDouble()
            else: 
                E_tr_B = 0.5

            if metadata["degradation_factors"].Has("G_patr_B"):
                G_patr_B = metadata["degradation_factors"]["G_patr_B"].GetDouble()
            else:
                G_patr_B = 0.5
        else:
            KratosMultiphysics.Logger.PrintInfo("::[Puck Analysis]::", "No Degradation Factors found in metadata. Using default values.")
            E_tr_A   = 0
            G_patr_A = 0
            E_tr_B   = 0.5
            G_patr_B = 0.5

        degradationfactor_E_tr_A   = E_tr_A
        degradationfactor_G_patr_A = G_patr_A
        degradationfactor_E_tr_B   = E_tr_B
        degradationfactor_G_patr_B = G_patr_B
  

# start extracting material data                
        idx = 0 # global index for mapping values to their respective rows across all elements
# iterate through elements using a loop  
        for element in sub_model_part.Elements: 

            properties = element.Properties  

            if not properties.Has(layers_variable): # check whether SHELL_ORTHOTROPIC_LAYERS exists as a property in the element
                print("Warning: Element ", element.Id, " has no SHELL_ORTHOTROPIC_LAYERS property defined. Skipping material extraction for this element.")
                continue
            
            layers_matrix = properties.GetValue(layers_variable) # a general variable that now contains the material data for the entire element
            
            rows = layers_matrix.Size1() # test-model: 16 Zeilen 
            cols = layers_matrix.Size2() # test-model: 8 Spalten
                        
# layers_matrix =  Kratos.Matrix, transferring the entries to NumPy-Vector
# iterate through plies using a loop 
            for PlyId in range(rows):

                # row -> NumPy-Vector
                current_layer_data = np.zeros(cols)
                
                for col in range(cols):
                    current_layer_data[col] = layers_matrix[PlyId, col]


                if current_layer_data.size < 16:
                    print(f"Fehler: zu wenig Daten für Ply {PlyId+1} in Element {element.Id}")
                    continue
# filling the vectors                
                thickness  = current_layer_data[0]
                alpha      = current_layer_data[1]

                R_pa_t    = current_layer_data[9]
                R_pa_c    = current_layer_data[10]
                R_tr_t    = current_layer_data[11]
                R_tr_c    = current_layer_data[12]
                R_trpa    = current_layer_data[13]


                E_pa    = current_layer_data[3]
                E_tr    = current_layer_data[4]
                nu_12  = current_layer_data[5]
                G_trpa  = current_layer_data[6]


                element_id[idx:idx+2] = element.Id
                ply_id[idx:idx+2]     = PlyId + 1
                ply_side[idx:idx+2]   = [0, 1]

                ply_angle[idx:idx+2] = alpha
                ply_thickness[idx //2] = thickness

                strength_R_pa_t[idx:idx+2] = R_pa_t
                strength_R_pa_c[idx:idx+2] = R_pa_c
                strength_R_tr_t[idx:idx+2] = R_tr_t
                strength_R_tr_c[idx:idx+2] = R_tr_c
                strength_R_trpa[idx:idx+2] = R_trpa

                inclination_p_trtr_c[idx:idx+2] = p_trtr_c
                inclination_p_trpa_t[idx:idx+2] = p_trpa_t
                inclination_p_trpa_c[idx:idx+2] = p_trpa_c

                youngs_modul_E_pa[idx:idx+2]    = E_pa
                youngs_modul_E_tr[idx:idx+2]    = E_tr
                poissons_ratio_nu_12[idx:idx+2] = nu_12
                shear_modulus_G_trpa[idx:idx+2] = G_trpa

                idx += 2

# filling the container
        composite_material = PanelCompositeMaterial(element_id, 
                                                                    ply_id, 
                                                                    ply_side,
                                                                    ply_thickness, 
                                                                    ply_angle, 
                                                                    strength_R_pa_t, 
                                                                    strength_R_pa_c,
                                                                    strength_R_tr_t,
                                                                    strength_R_tr_c, 
                                                                    strength_R_trpa, 
                                                                    inclination_p_trtr_c, 
                                                                    inclination_p_trpa_t, 
                                                                    inclination_p_trpa_c, 
                                                                    youngs_modul_E_pa, 
                                                                    youngs_modul_E_tr, 
                                                                    poissons_ratio_nu_12, 
                                                                    shear_modulus_G_trpa,
                                                                    degradationfactor_E_tr_A, 
                                                                    degradationfactor_G_patr_A, 
                                                                    degradationfactor_E_tr_B, 
                                                                    degradationfactor_G_patr_B, 
                                                                    )


        return composite_material

    def ExtractLaminate(self, sub_model_part, material_type: str, material_axis_offset: float) -> Laminate:
        """Laminate (CLT) of a shell component, ply angles in the component frame.

        Args:
            sub_model_part: Shell elements of the component; all must share one Properties.
            material_type: "isotropic" or "composite" (user input, checked against the Properties).
            material_axis_offset: Angle [rad] from the component x_axis to LOCAL_MATERIAL_AXIS_1.
        """
        properties = self._GetUniqueProperties(sub_model_part)
        has_layers = properties.Has(SMA.SHELL_ORTHOTROPIC_LAYERS)

        if properties.Has(SMA.SHELL_OFFSET) and properties[SMA.SHELL_OFFSET] != 0.0:
            raise RuntimeError(f"'{sub_model_part.Name}': SHELL_OFFSET is not supported, CLT is evaluated about the mid-plane.")

        if material_type == "composite":
            if not has_layers:
                raise RuntimeError(f"'{sub_model_part.Name}': material_type 'composite' requires SHELL_ORTHOTROPIC_LAYERS.")
            matrix = properties[SMA.SHELL_ORTHOTROPIC_LAYERS]
            layers = np.array([[matrix[i, j] for j in range(matrix.Size2())] for i in range(matrix.Size1())])
            laminate = Laminate.FromOrthotropicLayers(layers, np.degrees(material_axis_offset))
        else:
            if has_layers:
                raise RuntimeError(f"'{sub_model_part.Name}': material_type 'isotropic' but SHELL_ORTHOTROPIC_LAYERS is defined.")
            laminate = Laminate.FromIsotropic(properties[KratosMultiphysics.YOUNG_MODULUS],
                                              properties[KratosMultiphysics.POISSON_RATIO],
                                              properties[KratosMultiphysics.THICKNESS])

        KratosMultiphysics.Logger.PrintInfo("Panel",
            f"'{sub_model_part.Name}': t={laminate.thickness:.4g}, symmetric={laminate.IsSymmetric()}")
        return laminate

    @staticmethod
    def _GetUniqueProperties(sub_model_part):
        """A component is one laminate: all elements must share the same Properties."""
        properties = {element.Properties.Id: element.Properties for element in sub_model_part.Elements}
        if len(properties) != 1:
            raise RuntimeError(f"'{sub_model_part.Name}' has {len(properties)} different Properties, expected one laminate.")
        return next(iter(properties.values()))
