import KratosMultiphysics
from KratosMultiphysics.StructuralMechanicsApplication.structural_components.panel.panel_data import PanelCompositeMaterial
from KratosMultiphysics.StructuralMechanicsApplication.structural_components.panel.panel_helper import NumberOfPlies 

import KratosMultiphysics.StructuralMechanicsApplication as SMA
import numpy as np

class PanelMaterialExtractor:

    def ExtractCompositeMaterial(self, sub_model_part, metadata):
        #TODO: Implement this function
        
        # Check ob SubModelPart Elemente enthält
        if sub_model_part.NumberOfElements() == 0:
            raise RuntimeError("Panel submodelpart contains no elements")
        
        layers_variable = SMA.SHELL_ORTHOTROPIC_LAYERS

            
        total_plies = 0 #Gesamtzahl der Plies auf 0 setzen 
        for element in sub_model_part.Elements:
            total_plies += NumberOfPlies(element.Properties) # RUFT METHODE auf ! - Gesamtzahl der Plies über alle Elemente hinweg, um die Größe der Arrays zu bestimmen Hier 128     

        if total_plies == 0: #Error fals keine Plies definiert sind, da sonst die Arrays nicht initialisiert werden können  
            raise RuntimeError("No plies defined in any element of the submodelpart")

            

        # Leere Arrays anlegen in der zuvor ermittelten Länge. Doppelte Länge das sie später mit Spannungen je Ply (Top und Bottom) zusammennpassen
        element_id       = np.zeros(2*total_plies, dtype=int)
        ply_id           = np.zeros(2*total_plies, dtype=int)
        ply_side         = np.zeros(2*total_plies, dtype=int)
        ply_thickness    = np.zeros(total_plies) #PlyThickness brauche ich das doppelte Format nicht daher nur ein Ply ein eintrag kein top und Bottom
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

#einmalig werden die Neigungsparameter eingelesen da für gesamtes LAminat gleich, erst später dann in Vektor in richtiger diomension geschrieben.
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

#Einmalig werden die Degradationsfaktoren eingelesen
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
  

#Start Matrialdaten auslesen                
        idx = 0 #globaler Index für die Zuordnung der Werte zu den jeweiligen Plies über alle Elemente hinweg
        for element in sub_model_part.Elements: #durch die Elemente iterieren, um die Materialdaten der Plies zu extrahieren

            properties = element.Properties  

            if not properties.Has(layers_variable): # Check ob SHELL_ORTHOTROPIC_LAYERS as Property im Element existiert
                print("Warning: Element ", element.Id, " has no SHELL_ORTHOTROPIC_LAYERS property defined. Skipping material extraction for this element.")
                continue

            layers_matrix = properties.GetValue(layers_variable) # eine allgemeine Variable wo jetzt die Materialdaten des gesamten Elements beinhaltet
            
            rows = layers_matrix.Size1() #16 Zeilen            print("layers_matrix for Element",layers_matrix )
            cols = layers_matrix.Size2() #8 Spalten
                        
#Layers_matrix ist eine Kratos.Matrix und daher kann ich nicht einfach auf zeile 1 zugreifen. Gelöst durch Zeile nehmen und durch Spalten cols itterieren und so meine erwartete Matrix füllen die ich dann auslesen kann. 
            for PlyId in range(rows):

                # Zeile -> NumPy-Vektor
                current_layer_data = np.zeros(cols)
                
                for col in range(cols):
                    current_layer_data[col] = layers_matrix[PlyId, col]


                if current_layer_data.size < 16:
                    print(f"Fehler: zu wenig Daten für Ply {PlyId+1} in Element {element.Id}")
                    continue
                
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