import KratosMultiphysics
from KratosMultiphysics.StructuralMechanicsApplication.structural_components.panel.panel_data import PanelCompositeResponse
from KratosMultiphysics.StructuralMechanicsApplication.structural_components.panel.panel_helper import NumberOfPlies 
import KratosMultiphysics.StructuralMechanicsApplication as SMA
import numpy as np

class PanelResponseExtractor:

    def ExtractCompositeResponse(self, sub_model_part, material):
               
        stresses = SMA.SHELL_ORTHOTROPIC_STRESS_THROUGH_THICKNESS
        force = SMA.SHELL_FORCE_GLOBAL
        moment = SMA.SHELL_MOMENT_GLOBAL
        process_info = sub_model_part.ProcessInfo      
        
        total_plies = 0
        for element in sub_model_part.sub_model_part.Elements:
            total_plies += NumberOfPlies(element.Properties) # Gesamtzahl der Plies über alle Elemente hinweg, um die Größe der Arrays zu bestimmen        
        
        if total_plies == 0:
            raise RuntimeError("No plies defined in any element of the submodelpart")
        
        # Leere Arrays anlegen in der zuvor ermittelten Länge. Doppelte Länge da je Ply (Top und Bottom)
        sigma_1 = np.zeros(2*total_plies)
        sigma_2 = np.zeros(2*total_plies)
        tau_21  = np.zeros(2*total_plies)

        num_elements = sub_model_part.sub_model_part.NumberOfElements()

        force_x = np.zeros(num_elements)  
        force_y = np.zeros(num_elements)
        force_xy  = np.zeros(num_elements)

        moment_x = np.zeros(num_elements)
        moment_y = np.zeros(num_elements)
        moment_xy  = np.zeros(num_elements)


        plies_element = total_plies/num_elements # Speichern der Gesamtzahl der Plies als Attribut der Klasse 
        num_elements = num_elements  # Speichern der Gesamtzahl der Elemente als Attribut der Klasse


        
        idx_element = 0 #globaler index zählt die Anzahl der Elemente, um die Kräfte in den richtigen Vektor zu schreiben

        idx = 0 #globaler Index für die Zuordnung der Werte zu den jeweiligen Plies über alle Elemente hinweg
        for element in sub_model_part.Elements:

            #Anzahl der Gausspunkte je Element bestimmen
            integration_points = element.GetIntegrationPoints()
            number_integrationpoints = len(integration_points)
            #print("Number of Integration Points for Element", element.Id, ":", number_integrationpoints)

            # Summe initialisieren (None = "noch nichts addiert")
            stresses_sum = None
            force_sum = None
            moment_sum = None         

#Avergagen Forces
            # Schleife über die Gausspunkte: 0, 1, 2, ... nacheinander einlesen und aufsummieren 
            for G_punkt in range(number_integrationpoints):
                force_gp = element.CalculateOnIntegrationPoints(force, process_info)[G_punkt]  # Wert am aktuellen Gausspunkt
                
                if force_sum is None:
                    force_sum = force_gp   # Gausspunkt 0 als Startwert 
                else:
                    force_sum = force_sum + force_gp   # Die Restlichen Gausspunkte darauf aufsummieren

            
            #Mittelwert bilden
            force_gp_elementaxis = force_sum / number_integrationpoints           

            #Auslesen der Dehnungen aus dem Dehnungstensor und Vektoren mit allen Elementen erstellen

            force_x[idx_element] = force_gp_elementaxis[0, 0]
            force_y[idx_element] = force_gp_elementaxis[1, 1]
            force_xy[idx_element]  = force_gp_elementaxis[0, 1] 
            
#Avergagen Moments
            # Schleife über die Gausspunkte: 0, 1, 2, ... nacheinander einlesen und aufsummieren
            for G_punkt in range(number_integrationpoints):
                moment_gp = element.CalculateOnIntegrationPoints(moment, process_info)[G_punkt]  # Wert am aktuellen Gausspunkt

                if moment_sum is None:
                    moment_sum = moment_gp   # Gausspunkt 0 als Startwert 
                else:
                    moment_sum = moment_sum + moment_gp   # Die Restlichen Gausspunkte darauf aufsummieren

            #Mittelwert bilden
            moment_gp_elementaxis = moment_sum / number_integrationpoints           

            #Auslesen der Dehnungen aus dem Dehnungstensor und Vektoren mit allen Eleeḿenten erstellen
            moment_x[idx_element] = moment_gp_elementaxis[0, 0]
            moment_y[idx_element] = moment_gp_elementaxis[1, 1]
            moment_xy[idx_element]  = moment_gp_elementaxis[0, 1] 

            idx_element += 1
 
#Averagen Stresses
            # Schleife über die Gausspunkte: 0, 1, 2, ... nacheinander einlesen und aufsummieren DER STRESSES
            for G_punkt in range(number_integrationpoints):
                stresses_gp = element.CalculateOnIntegrationPoints(stresses, process_info)[G_punkt]  # Wert am aktuellen Gausspunkt

                if stresses_sum is None:
                    stresses_sum = stresses_gp   # Gausspunkt 0 als Startwert 
                else:
                    stresses_sum = stresses_sum + stresses_gp   # Die Restlichen Gausspunkte darauf aufsummieren

            #Mittelwert bilden
            stresses_gp_elementaxis = stresses_sum / number_integrationpoints

            # print("STRESSES", stresses_gp_elementaxis)
            rows = stresses_gp_elementaxis.Size1()             # Schleife über die Gausspunkte: 0, 1, 2, ... nacheinander einlesen und aufsummieren DER STRESSES
            cols = stresses_gp_elementaxis.Size2() 

            for PlyId in range(rows): #durch die Plies je Element iterieren transformieren
                            #doppelt so lange durchitterieren um top und bottom aller abzugreifen

                # Zeile -> NumPy-Vektor
                current_stress_data = np.zeros(cols)
               
                for col in range(cols):
                    current_stress_data[col] = stresses_gp_elementaxis[PlyId, col]

                    current_stress_data_2D = current_stress_data[:3] #Verkürzen auf 3 Einträge wegen ShellThick liefert z.B 5 einträge wo die letzten 2 egal sind bzw unbekannt und quasi 0
#
#
#
                ply_angle_rad = np.radians(material.ply_angle[idx])
#
#
#
                cos = np.cos(ply_angle_rad)
                sin = np.sin(ply_angle_rad)

                T = np.array([
                    [cos**2, sin**2,  2*sin*cos],
                    [sin**2, cos**2, -2*sin*cos],
                    [-sin*cos, sin*cos, cos**2 - sin**2]
                ])
               
                stress_state_plyaxis = T @ current_stress_data_2D # Schritt 2 nach 3 transformation der Spannungen von Elementachsen in Plyachsen

                sigma_1[idx] = stress_state_plyaxis[0] #Schritt 4 zusammenführen und aufteilen
                sigma_2[idx] = stress_state_plyaxis[1]
                tau_21[idx] = stress_state_plyaxis[2]

                idx += 1 #Am Ende ein langer Vektor mit allen Ply-Spannungen über allen Elemente  


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