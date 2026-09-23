#  ____             _       |
# |  _ \ _   _  ___| | __   |
# | |_) | | | |/ __| |/ /   |
# |  __/| |_| | (__|   <    |
# |_|    \__,_|\___|_|\_\   ANALYSIS
#
#  Main authors:    Lucas Rimpl
#  Co-authors:      Tobias Siemer
#
import numpy as np
from KratosMultiphysics.StructuralMechanicsApplication.handbook_methods.analysis_result import AnalysisResult
from KratosMultiphysics.StructuralMechanicsApplication.handbook_methods.method_base import HandbookMethod
import KratosMultiphysics.StructuralMechanicsApplication as SMA

class PuckAnalysis:

    def IsApplicable(self, structural_component):
        return True

    def Evaluate(self, structural_component):
        self.PuckAnalysis(structural_component)

    @staticmethod
    def PuckFF(sigma_1: np.ndarray, R_pa_t: np.ndarray, R_pa_c: np.ndarray):
        """Parameter Description

        Returns:
            t:          tensile
            c:          compressive
            pa:         parallel
            tr:         transverse

            sigma_1:    normal stress x1-axis

            R_pa_t:     fiber parallel tensile strength for uniaxial sigma||t-stressing
            R_pa_c:     fiber parallel compressive strength for uniaxial sigma||c-stressing
            f_E_FF:     stress exposure Fiber Fracturen (FF)
            RF_FF:      Reserve Factor Fiber Fracture (FF)
        """
    
        R_pa_ct = np.where(sigma_1 < 0, R_pa_c, R_pa_t)

        failure_mode_FF = np.where(sigma_1 < 0, "FF Compression","FF Tension")

        f_E_FF = np.abs(sigma_1/R_pa_ct)

        RF_FF  = 1/f_E_FF

        return RF_FF,f_E_FF,failure_mode_FF
   
    
    @staticmethod
    def PuckIFF(sigma_2: np.ndarray, tau_21: np.ndarray, R_tr_t: np.ndarray, R_tr_c: np.ndarray, R_trpa: np.ndarray, p_trtr_c: np.ndarray, p_trpa_t: np.ndarray, p_trpa_c: np.ndarray):
           
        """Parameter Description

        Returns:
            t:          tensile
            c:          compressive
            tr:         transverse
            pa:         parallel

            sigma_1:    normal stress x1-axis
            sigma_2:    normal stress x2-axis
            tau_21:     sher stress
            tau_21_c:   shear stress transition Point Mode B to Mode C
            
            R_tr_t:     transverse tensile strength for uniaxial sigma⊥t-stressing
            R_tr_c:     transverse compressive strength for uniaxial sigma⊥c-stressing
            R_trtr:     transverse shear strength for pure tau⊥⊥-stressing
            R_trtr_A:   transverse shear strength for pure tau⊥⊥-stressing ---fracture resistance of the action plane
            R_trpa:     longitudinal shear strength for pure tau⊥||-stressing
            p_trtr_c:   inclination parameter
            p_trpa_t:   inclination parameter
            p_trpa_c:   inclination parameter
            theta_fp:   angle fracture plane
            f_E_IFF:    stress exposure Inter Fiber Fracturen (IFF)
            RF_IFF:     Reserve Factor Inter Fiber Fracture (IFF)
        """
        R_trtr_A = R_tr_c/(2*(1+p_trtr_c))
        tau_21_c = R_trpa*np.sqrt(1+2*p_trtr_c)

        # Bedingung für failure mode
        cond_A = sigma_2 >= 0
        cond_B = (sigma_2 < 0) & (np.abs(sigma_2 / tau_21) <= (R_trtr_A / np.abs(tau_21_c)))
        cond_C = (sigma_2 < 0) & (np.abs(tau_21 / sigma_2) <= (np.abs(tau_21_c) / R_trtr_A))
        
        # Stress exposure ratio ratio Mode A,B,C
        
        # Mode A
        f_E_A = np.sqrt((tau_21/R_trpa)**2+(1-p_trpa_t*R_tr_t/R_trpa)**2*(sigma_2/R_tr_t)**2)+p_trpa_t*sigma_2/R_trpa
        
        # Mode B
        f_E_B = 1/R_trpa*(np.sqrt(tau_21**2+(p_trpa_c*sigma_2)**2)+p_trpa_c*sigma_2)
        
        # Mode C
        f_E_C = ((tau_21/(2*(1+p_trtr_c)*R_trpa))**2+(sigma_2/R_tr_c)**2)*R_tr_c/-sigma_2

        # Auswahl failure Mode
        f_E_IFF = np.where(cond_A, f_E_A, np.where(cond_B, f_E_B, np.where(cond_C, f_E_C, np.nan)))

        # Vector failure mode 
        failure_mode_IFF = np.where(cond_A,"IFF Mode A",np.where(cond_B, "IFF Mode B", np.where(cond_C, "IFF Mode C", "NO IFF Mode")))

        # Bruchwinkel berechnen (nur für Mode C) np.clip als Abnsicherung falls durch Rundung Wurzel negativ oder arccos größer 1 werte bekommt
        theta_fp = np.degrees(np.arccos(np.sqrt(np.clip((1/(2*(1+p_trtr_c)))*(((R_trtr_A*tau_21)/(R_trpa*sigma_2))**2+1),0.0,1.0))))

        #  Zuweisung Bruchwinkel
        theta_fp_IFF = np.where(cond_C, theta_fp, np.degrees(0.0))

        RF_IFF = 1/f_E_IFF

        return RF_IFF, f_E_IFF, failure_mode_IFF, theta_fp_IFF

    


    @staticmethod
    def PuckDegradation (ElementID: np.array, PlyID: np.array, topbot: np.array, RF_FF_deg: np.ndarray, failure_mode_FF: np.ndarray, RF_IFF_deg: np.ndarray, failure_mode_IFF: np.ndarray, force_x: np.ndarray, force_y: np.ndarray, force_xy: np.ndarray, moment_x: np.ndarray, moment_y: np.ndarray, moment_xy: np.ndarray, E_pa: np.ndarray, E_tr: np.ndarray, nu_12: np.ndarray, G_trpa: np.ndarray, PlyAngle: np.ndarray, PlyThickness: np.ndarray, plies_element: int, num_elements: int,
                         R_pa_t: np.ndarray, R_pa_c: np.ndarray, R_tr_t: np.ndarray, R_tr_c: np.ndarray, R_trpa: np.ndarray, p_trtr_c: np.ndarray, p_trpa_t: np.ndarray, p_trpa_c: np.ndarray, E_tr_A: int, G_patr_A: int, E_tr_B: int, G_patr_B:int, DegSteps_max: int):
        """Parameter Description

        Returns:
            t:          tensile                    generic.degradationfactor_E_tr_A,
                    generic.degradationfactor_G_patr_A,
                    generic.degradationfactor_E_tr_B,
                    generic.degradationfactor_G_patr_A,
            c:          compressive
            pa:         parallel
            tr:         transverse    

            deg:        Symbol für Degradationsrechnung
            E_tr:     Young's modulus in transverse direction
            G_trpa:     Shear modulus in transverse direction
        """
        PuckDegradationIFF = np.zeros((plies_element*2, num_elements), dtype=object) #Ergebnisvektor der RF_IFF für Degradation
        element_mit_degradation = np.zeros(num_elements) #Ergebnisvektor der Elemente die Degradation erfahren haben
        Puck_degradations_analyse = np.zeros((plies_element*2, 4), dtype=object)
        

       #RF nach Element sortieren
        RF_FF_deg = np.reshape(RF_FF_deg.copy(), (plies_element*2, num_elements), order='F') #ACHTUNG Umwandlung von Vektor zu Matrix
        failure_mode_FF_deg= np.reshape(failure_mode_FF.copy(), (plies_element*2, num_elements), order='F')

        RF_IFF_deg = np.reshape(RF_IFF_deg.copy(), (plies_element*2, num_elements), order='F')
        failure_mode_IFF_deg = np.reshape(failure_mode_IFF.copy(), (plies_element*2, num_elements), order='F')

        #Richtige Dimension der Elementparameter etc füt Puck
        ElementID = np.reshape(ElementID, (plies_element*2, num_elements), order='F')
        PlyID     = np.reshape(PlyID, (plies_element*2, num_elements), order='F')
        topbot    = np.reshape(topbot, (plies_element*2, num_elements), order='F')

        PlyAngle = np.reshape(PlyAngle, (plies_element*2, num_elements), order='F')
        PlyThickness = np.reshape(PlyThickness, (plies_element, num_elements), order='F')

        R_pa_t = np.reshape(R_pa_t, (plies_element*2, num_elements), order='F')
        R_pa_c = np.reshape(R_pa_c, (plies_element*2, num_elements), order='F')
        R_tr_t = np.reshape(R_tr_t, (plies_element*2, num_elements), order='F')
        R_tr_c = np.reshape(R_tr_c, (plies_element*2, num_elements), order='F')
        R_trpa = np.reshape(R_trpa, (plies_element*2, num_elements), order='F')

        p_trtr_c = np.reshape(p_trtr_c, (plies_element*2, num_elements), order='F')
        p_trpa_t = np.reshape(p_trpa_t, (plies_element*2, num_elements), order='F')
        p_trpa_c = np.reshape(p_trpa_c, (plies_element*2, num_elements), order='F')

        E_pa = np.reshape(E_pa, (plies_element*2, num_elements), order='F')
        E_tr = np.reshape(E_tr, (plies_element*2, num_elements), order='F')
        nu_12 = np.reshape(nu_12, (plies_element*2, num_elements), order='F')
        G_trpa = np.reshape(G_trpa, (plies_element*2, num_elements), order='F')

#Implementiert nach structure sizing Last Ply Failure LPF            
        for element in range(num_elements): #Iterieren über alle Elemente
            Degradation_abbruch = False
            Degradation_number = 1

            while True:
                
                print("Degradationsanalyse Nummer",Degradation_number,"Startet für Element ", ElementID[0,element]) 

                if np.max(( RF_IFF_deg[ :, element])) < 1 : #Gesamtes Ply bricht
                    print ("Element ", ElementID[:,element], "Abruchkriterium erreicht gesamtes Element RF IFF < 1")
                    Degradation_abbruch = True
                    break
                            
                elif np.min((RF_IFF_deg[ :, element])) >= 1 : #Einzelnes Ply bricht nicht alles safe
                    print("Element ", ElementID[:,element], " sicher, alle RF IFF > 1, keine Degradation notwendig" )
                    Degradation_abbruch = True
                    break


                Q_vektor = np.zeros((plies_element, 3, 3))
                A = np.zeros((3, 3))
                B = np.zeros((3, 3))
                D = np.zeros((3, 3))
#Von Top Bottom betrachtung auf Ply betrachtung umstellen um ABD richtig zu berechnen
                RF_top = RF_IFF_deg[0::2, element]
                RF_bottom = RF_IFF_deg[1::2, element]

                RF_IFF_ABD = np.minimum(RF_top, RF_bottom) #RF IFF min der Lage von top und Bottom ist relevant

                #print("RF_ABD", RF_IFF_ABD)

                failure_mode_IFF_ABD = np.where(
                    RF_top <= RF_bottom,
                    failure_mode_IFF_deg[0::2, element],
                    failure_mode_IFF_deg[1::2, element])

                for ply in range(plies_element): 

                #FFF Check lässt sich hinzufügen !!
                    # if RF_FF_deg[ply, element]<1 : #Optional hinzuschalten von FFF- Die FF dieser Ply führen direkt zum Abbruch
                    #     print("In Element ", ElementID[ply,element], "Ply ID: ", PlyID[ply,element], "Top=0/Bottom=1:" ,topbot[ply,element], "FF Failure detected and not allowed")
                    #     print("RF_FF: ", RF_FF_deg[ply, element], "Failure Mode: ", failure_mode_FF_deg[ply, element])
                    #     Degradation_abbruch = True
                    #     break   

                    if RF_IFF_ABD[ply] >= 1:
                    # Berechnung Q-Matrix - UHRSPRUNGSWETE FÜR E UND G - O = Original

                        E_tr_degO = E_tr[2*ply, element]
                        G_trpa_degO = G_trpa[2*ply, element]

                        nu_12_degO = nu_12[2*ply, element]
                        E_pa_degO = E_pa[2*ply, element]

                        nu_21 = nu_12_degO * E_tr_degO/E_pa_degO
                        denom = 1 - nu_12_degO* nu_21     

                        Q11O = E_pa_degO/denom
                        Q22O = E_tr_degO/denom
                        Q12O = nu_12_degO * E_tr_degO/ denom
                        Q66O = G_trpa_degO
                        Q_O =  np.array([
                                [Q11O, Q12O, 0],
                                [Q12O, Q22O, 0],
                                [0, 0, Q66O]
                                ])

                        ply_angle_rad = np.radians(PlyAngle[2*ply,element])

                        cos = np.cos(ply_angle_rad)
                        sin = np.sin(ply_angle_rad)

                        Q11, Q12, Q22, Q66 = Q_O[0,0], Q_O[0,1], Q_O[1,1], Q_O[2,2]

                        Qb11O = Q11*cos**4 + 2*(Q12 + 2*Q66)*sin**2*cos**2 + Q22*sin**4
                        Qb22O = Q11*sin**4 + 2*(Q12 + 2*Q66)*sin**2*cos**2 + Q22*cos**4
                        Qb12O = (Q11 + Q22 - 4*Q66)*sin**2*cos**2 + Q12*(cos**4 + sin**4)
                        Qb66O = (Q11 + Q22 - 2*Q12 - 2*Q66)*sin**2*cos**2 + Q66*(sin**4 + cos**4)


                        Qb_O =  np.array([
                                        [Qb11O, Qb12O, 0],
                                        [Qb12O, Qb22O, 0],
                                        [0, 0, Qb66O]
                                    ])

                        Q_vektor[ply] = Qb_O
                        
                    elif failure_mode_IFF_ABD[ply] == "IFF Mode A":
                    # Berechnung Q-Matrix - WERTE BEI DEGRADATION FAILURE MODE A FÜR E UND G - A = MODE A
                        E_tr_degA = E_tr_A * E_tr[2*ply, element]
                        G_trpa_degA = G_patr_A * G_trpa[2*ply, element]

                        nu_12_degA = nu_12[2*ply, element]
                        E_pa_degA = E_pa[2*ply, element]
                        nu_21 = nu_12_degA * E_tr_degA/E_pa_degA
                        denom = 1 - nu_12_degA* nu_21     
                        
                        Q11A = E_pa_degA/denom
                        Q22A = E_tr_degA/denom
                        Q12A = nu_12_degA * E_tr_degA/ denom
                        Q66A = G_trpa_degA

                        Q_A =  np.array([
                                [Q11A, Q12A, 0],
                                [Q12A, Q22A, 0],
                                [0, 0, Q66A]
                                ])

                        ply_angle_rad = np.radians(PlyAngle[2*ply,element])

                        cos = np.cos(ply_angle_rad)
                        sin = np.sin(ply_angle_rad)

                        Q11, Q12, Q22, Q66 = Q_A[0,0], Q_A[0,1], Q_A[1,1], Q_A[2,2]

                        Qb11A = Q11*cos**4 + 2*(Q12 + 2*Q66)*sin**2*cos**2 + Q22*sin**4
                        Qb22A = Q11*sin**4 + 2*(Q12 + 2*Q66)*sin**2*cos**2 + Q22*cos**4
                        Qb12A = (Q11 + Q22 - 4*Q66)*sin**2*cos**2 + Q12*(cos**4 + sin**4)
                        Qb66A = (Q11 + Q22 - 2*Q12 - 2*Q66)*sin**2*cos**2 + Q66*(sin**4 + cos**4)


                        Qb_A = np.array([
                                        [Qb11A, Qb12A, 0],
                                        [Qb12A, Qb22A, 0],
                                        [0, 0, Qb66A]
                                    ])

                        Q_vektor[ply] = Qb_A                        


                    elif failure_mode_IFF_ABD[ply] == "IFF Mode B":
                    # Berechnung Q-Matrix - WERTE BEI DEGRADATION FAILURE MODE B FÜR E UND G - B = MODE B
                        E_tr_degB = E_tr_B * E_tr[2*ply, element]
                        G_trpa_degB = G_patr_B * G_trpa[2*ply, element]

                        nu_12_degB = nu_12[2*ply, element]
                        E_pa_degB = E_pa[2*ply, element]
                        nu_21 = nu_12_degB * E_tr_degB/E_pa_degB
                        denom = 1 - nu_12_degB* nu_21     
                        
                        Q11B = E_pa_degB/denom
                        Q22B = E_tr_degB/denom
                        Q12B = nu_12_degB * E_tr_degB/ denom
                        Q66B = G_trpa_degB
                        
                        Q_B =  np.array([
                                [Q11B, Q12B, 0],
                                [Q12B, Q22B, 0],
                                [0, 0, Q66B]
                                ])

                        ply_angle_rad = np.radians(PlyAngle[2*ply,element])

                        cos = np.cos(ply_angle_rad)
                        sin = np.sin(ply_angle_rad)

                        Q11, Q12, Q22, Q66 = Q_B[0,0], Q_B[0,1], Q_B[1,1], Q_B[2,2]

                        Qb11B = Q11*cos**4 + 2*(Q12 + 2*Q66)*sin**2*cos**2 + Q22*sin**4
                        Qb22B = Q11*sin**4 + 2*(Q12 + 2*Q66)*sin**2*cos**2 + Q22*cos**4
                        Qb12B = (Q11 + Q22 - 4*Q66)*sin**2*cos**2 + Q12*(cos**4 + sin**4)
                        Qb66B = (Q11 + Q22 - 2*Q12 - 2*Q66)*sin**2*cos**2 + Q66*(sin**4 + cos**4)

                        Qb_B = np.array([
                                        [Qb11B, Qb12B, 0],
                                        [Qb12B, Qb22B, 0],
                                        [0, 0, Qb66B]
                                    ])

                        Q_vektor[ply] = Qb_B                        



                    elif failure_mode_IFF_ABD[ply] == "IFF Mode C":
                        #ABBRUCH FALL NICHT ERLAUBT
                        print("In Element ", ElementID[ply,element], "Ply ID: ", PlyID[ply,element],"IFF Mode C detected and not allowed")
                        Degradation_abbruch = True
                        break

                if Degradation_abbruch:
                    break

                thickness_gesamt = np.sum(PlyThickness[:, element])  # Summe der PlyThickness für jedes Element = Laminat Dicke
                z_interfaces = -thickness_gesamt / 2 + np.concatenate(([0], np.cumsum(PlyThickness[:, element])))  # Berechnung der z-Positionen der Ply-Schnittstellen und anschließend befüllen von der mitte aus, dafür sorgt np.concatenate[0] mit extra 0 im array als Platz für z=0
                z = np.column_stack((z_interfaces[:-1], z_interfaces[1:])).flatten() #Sortieren, dass Top Bottom Top Bottom, .. reihenfolge ohne überschneidene werte nur einfach 

                z_ABD = z_interfaces

                # print("z", z)
                # print("z Schnittstellen", z_ABD)


                for k in range(plies_element):

                    A += Q_vektor[k] * (z_ABD[k+1] - z_ABD[k])
                    B += 0.5 * Q_vektor[k] * (z_ABD[k+1]**2 - z_ABD[k]**2)
                    D += (1/3) * Q_vektor[k] * (z_ABD[k+1]**3 - z_ABD[k]**3)

                ABD = np.block([
                                [A, B],
                                [B, D]
                            ])    
                
#3x1 Vektor der Dehnungen am element im globalen KS
                Belastung = np.zeros(6)
                Belastung[0] = force_x[element]
                Belastung[1] = force_y[element]
                Belastung[2] = force_xy[element]
                Belastung[3] = moment_x[element]
                Belastung[4] = moment_y[element]
                Belastung[5] = moment_xy[element]

                Dehnung_Kruemmung = np.linalg.solve(ABD, Belastung)

                Dehnung_element_global = Dehnung_Kruemmung[:3]
                Kruemmung_element_global = Dehnung_Kruemmung[3:]

            

                Q_vektordoppelt = np.repeat(Q_vektor, 2, axis=0)

                sigma_1_deg = np.zeros(plies_element*2)
                sigma_2_deg = np.zeros(plies_element*2)
                tau_21_deg = np.zeros(plies_element*2)

                for i in range(plies_element*2): 
                    Dehnungen_element_ply = Dehnung_element_global + z[i] * Kruemmung_element_global 

                    ply_angle_rad = np.radians(PlyAngle[i,element])

                    #Drehmatrix Dehnunghen
                    cos = np.cos(ply_angle_rad)
                    sin = np.sin(ply_angle_rad)    

                    T = np.array([
                                                [cos**2, sin**2,  sin*cos],
                                                [sin**2, cos**2, -1*sin*cos],
                                                [-2*sin*cos, 2*sin*cos, cos**2 - sin**2]]
                                                )
                    Dehnung_element_plyaxis = T @ Dehnungen_element_ply

                    Spannungen_deg = Q_vektordoppelt[i] @ Dehnung_element_plyaxis

                    #Das Spannungstrio an der Ply des Elements wo Kritischer IFF aufgetreten ist
                    sigma_1_deg[i] = Spannungen_deg[0] 
                    sigma_2_deg[i] = Spannungen_deg[1]
                    tau_21_deg[i] = Spannungen_deg[2]

                #print(sigma_1_deg)
                    #Irgendwie in Puck und hin und  wieder Zurück
                R_pa_t_deg = R_pa_t[:, element]
                R_pa_c_deg = R_pa_c[:, element]
                R_tr_t_deg = R_tr_t[:, element]
                R_tr_c_deg = R_tr_c[:, element]
                R_trpa_deg = R_trpa [:, element]

                p_trtr_c_deg = p_trtr_c[:, element]
                p_trpa_t_deg = p_trpa_t[:, element]
                p_trpa_c_deg = p_trpa_c[:, element]

                


                RF_FF_deg_new,f_E_FF,failure_mode_FF_new = PuckAnalysis.PuckFF(sigma_1_deg, R_pa_t_deg, R_pa_c_deg)
        
                RF_IFF_deg_new,f_E_IFF,failure_mode_IFF_deg_new,theta_fp_IFF_new = PuckAnalysis.PuckIFF(sigma_2_deg, tau_21_deg, R_tr_t_deg, R_tr_c_deg, R_trpa_deg, p_trtr_c_deg, p_trpa_t_deg, p_trpa_c_deg)

                RF_FF_deg[:,element] = RF_FF_deg_new
                failure_mode_FF_deg[:,element] = failure_mode_FF_new

                RF_IFF_deg[:,element] = RF_IFF_deg_new
                failure_mode_IFF_deg[:,element] = failure_mode_IFF_deg_new

#Zur Ausgabe während der Analyse                        
                Puck_degradations_analyse[:,0] = RF_IFF_deg_new
                Puck_degradations_analyse[:,1] = f_E_IFF
                Puck_degradations_analyse[:,2] = failure_mode_IFF_deg_new
                Puck_degradations_analyse[:,3] = theta_fp_IFF_new

                print ("-RF_IFF-", "-f_E_IFF-", "-failure_mode_IFF-", "-theta_fp_IFF-")
                PuckAnalysis.PrintMatrix(Puck_degradations_analyse)

#Zur Ausgabe in gesamter Matrix
                PuckDegradationIFF[:,element] = RF_IFF_deg_new
                element_mit_degradation[element] = ElementID[ply,element]
                                        
                Degradation_number += 1

#Abbruch falls Anbruchkiriterium nicht erreicht wird und Degradationsanalyse zu lange dauert anzahl 10 sind beliebig gewählt
                if Degradation_number > DegSteps_max:
                    print("Maximale Anzahl von", DegSteps_max, "Schritten in der Degradationsanalyse erreicht")
                    break
                                
        return PuckDegradationIFF, element_mit_degradation


        
    #Methode zur strukturierten Ausgabe der Puck_Analyse Matrix
    @staticmethod
    def PrintMatrix(matrix):
        n_cols = matrix.shape[1]
        col_widths = [max(len(str(val)) for val in matrix[:, c]) for c in range(n_cols)]

        separator = "+" + "+".join("-" * (w + 2) for w in col_widths) + "+"

        print(separator)
        for row in matrix:
            line = "|"
            for val, w in zip(row, col_widths):
                line += f" {str(val):<{w}} |"
            print(line)
            print(separator)
        

    @staticmethod
    def PuckAnalysis(structural_component):

        # ensure numpy arrays
        ElementID = np.array(structural_component.material.element_id)
        PlyID = np.array(structural_component.material.ply_id)
        topbot = np.array(structural_component.material.ply_side) # 0 für Top und 1 für Bottom Seite der Lage

        PlyAngle = np.array(structural_component.material.ply_angle)
        PlyThickness = np.array(structural_component.material.ply_thickness)

        sigma_1 = np.array(structural_component.response.sigma_1)
        sigma_2 = np.array(structural_component.response.sigma_2)
        tau_21 = np.array(structural_component.response.tau_21)

        force_x = np.array(structural_component.response.force_x)
        force_y = np.array(structural_component.response.force_y)
        force_xy = np.array(structural_component.response.force_xy)

        moment_x = np.array(structural_component.response.moment_x)
        moment_y = np.array(structural_component.response.moment_y)
        moment_xy = np.array(structural_component.response.moment_xy)

        R_pa_t = np.array(structural_component.material.strength_R_pa_t)
        R_pa_c = np.array(structural_component.material.strength_R_pa_c)
        R_tr_t = np.array(structural_component.material.strength_R_tr_t)
        R_tr_c = np.array(structural_component.material.strength_R_tr_c)
        R_trpa = np.array(structural_component.material.strength_R_trpa)

        p_trtr_c = np.array(structural_component.material.inclination_p_trtr_c)
        p_trpa_t = np.array(structural_component.material.inclination_p_trpa_t)
        p_trpa_c = np.array(structural_component.material.inclination_p_trpa_c)

        E_pa   = np.array(structural_component.material.youngs_modul_E_pa)
        E_tr   = np.array(structural_component.material.youngs_modul_E_tr)
        nu_12  = np.array(structural_component.material.poissons_ratio_nu_12)
        G_trpa = np.array(structural_component.material.shear_modulus_G_trpa)

        plies_element = int(structural_component.response.plies_element)
        num_elements  = int(structural_component.response.num_elements)

        E_tr_A   = structural_component.material.degradationfactor_E_tr_A
        G_patr_A = structural_component.material.degradationfactor_G_patr_A
        E_tr_B   = structural_component.material.degradationfactor_E_tr_B
        G_patr_B = structural_component.material.degradationfactor_G_patr_B

        DegSteps_max = structural_component.material.max_degradation_steps

        
        RF_FF,f_E_FF,failure_mode_FF = PuckAnalysis.PuckFF(sigma_1, R_pa_t, R_pa_c)
        

        RF_IFF,f_E_IFF,failure_mode_IFF,theta_fp_IFF = PuckAnalysis.PuckIFF(sigma_2,tau_21, R_tr_t, R_tr_c, R_trpa, p_trtr_c, p_trpa_t, p_trpa_c)

        RF_FF_deg = RF_FF
        RF_IFF_deg = RF_IFF

        PuckDegradationIFF, element_mit_degradation = PuckAnalysis.PuckDegradation(ElementID, PlyID, topbot, RF_FF_deg, failure_mode_FF, RF_IFF_deg, failure_mode_IFF, force_x, force_y, force_xy, moment_x, moment_y, moment_xy, E_pa, E_tr, nu_12, G_trpa, PlyAngle, PlyThickness, plies_element, num_elements,
                                                               R_pa_t, R_pa_c, R_tr_t, R_tr_c, R_trpa, p_trtr_c, p_trpa_t, p_trpa_c, E_tr_A, G_patr_A, E_tr_B, G_patr_B, DegSteps_max)

        Puck_analyse = np.zeros((PlyID.size, 12), dtype=object) #PlyID ist schon ein vector mit doppelter größe wegen TOP und BOTTOM
 

        Puck_DegradationIFF = PuckDegradationIFF.flatten(order="F")
    
        # Das Array befüllen um alle Variablen des Puck-Kriteriums gebündelt auszugeben 
        Puck_analyse[:,0] = ElementID
        Puck_analyse[:,1] = PlyID
        Puck_analyse[:,2] = PlyAngle
        Puck_analyse[:,3] = topbot

        Puck_analyse[:,4] = RF_FF
        Puck_analyse[:,5] = f_E_FF
        Puck_analyse[:,6] = failure_mode_FF  

        Puck_analyse[:,7]  = RF_IFF
        Puck_analyse[:,8]  = f_E_IFF
        Puck_analyse[:,9]  = failure_mode_IFF
        Puck_analyse[:,10] = theta_fp_IFF 
        Puck_analyse[:,11] = Puck_DegradationIFF
   

        print ("FOLGENDE REIHENFOLGE DER SPALTEN:")

        print ("-ElementID-", "-PlyID-", "-PlyAngle-", "-topbot-", "-RF_FF-", "-f_E_FF-", "-failure_mode_FF-", "-RF_IFF-", "-f_E_IFF-", "-failure_mode_IFF-", "-theta_fp_IFF-", "-PuckDegradationIFF-")
    
        PuckAnalysis.PrintMatrix(Puck_analyse)   


