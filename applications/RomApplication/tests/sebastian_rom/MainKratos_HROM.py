import KratosMultiphysics
import KratosMultiphysics.RomApplication as romapp
from KratosMultiphysics.RomApplication.randomized_singular_value_decomposition import RandomizedSingularValueDecomposition
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis
from KratosMultiphysics.RomApplication.structural_mechanics_analysis_rom import StructuralMechanicsAnalysisROM
from matplotlib import pyplot as plt
import numpy as np
import json
from scipy.sparse import csr_matrix, load_npz
import pdb
from postprocessing import RunPostProcess

def rSVD(X,r,q,p):
    # Step 1: Sample column space of X with P matrix
    ny = X.shape[1]
    P = np.random.randn(ny,r+p)
    Z = X @ P
    for k in range(q):
        Z = X @ (X.T @ Z)

    Q, R = np.linalg.qr(Z,mode='reduced')

    # Step 2: Compute SVD on projected Y = Q.T @ X
    Y = Q.T @ X
    UY, S, VT = np.linalg.svd(Y,full_matrices=0)
    U = Q @ UY

    return U, S, VT

def ArrangeSnapshotMatrix(ResidualSnapshots):
    ### Building the Snapshot matrix ####
    for i in range (len(ResidualSnapshots)):
        if i == 0:
            SnapshotMatrix = ResidualSnapshots[i]
        else:
            SnapshotMatrix = np.c_[SnapshotMatrix,ResidualSnapshots[i]]
    return SnapshotMatrix

def FindDuplicateIndex(List_1,List_2):
    Index_list = []
    for i in range(len(List_1)):
        for j in range(len(List_2)):
            if List_1[i]==List_2[j]:
                Index_list.append(i)
                continue
    return Index_list

class RunHROM(StructuralMechanicsAnalysisROM):

    def ModifyInitialGeometry(self):
        """Here is the place where the HROM_WEIGHTS are assigned to the selected elements and conditions"""
        super().ModifyInitialGeometry()
        computing_model_part = self._solver.GetComputingModelPart()
        self.HR_elements = []
        self.HR_conditions = []
        self.HR_residuals = []
        ## Adding the weights to the corresponding elements
        with open('ElementsAndWeights.json') as f:
            HR_data = json.load(f)
            for key in HR_data["Elements"].keys():
                computing_model_part.GetElement(int(key)+1).SetValue(romapp.HROM_WEIGHT, HR_data["Elements"][key])
                self.HR_elements.append(int(key))
            for key in HR_data["Conditions"].keys():
                computing_model_part.GetCondition(int(key)+1).SetValue(romapp.HROM_WEIGHT, HR_data["Conditions"][key])
                self.HR_conditions.append(int(key))
            self.HR_elements.sort()
            self.HR_conditions.sort()
            restricted_residual_parameters = KratosMultiphysics.Parameters("""
        {
            "nodal_unknowns" : ["REACTION_X","REACTION_Y","REACTION_Z"],
            "number_of_dofs" : 1
        }""" )
        self.ResidualUtilityObject_Restricted_Residuals = romapp.RomModelPartUtility(self._GetSolver().GetComputingModelPart().GetSubModelPart("COMPUTE_HROM").GetSubModelPart("GENERIC_Interface_frame_1"), restricted_residual_parameters, KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme())
        self.ResidualUtilityObject_Main_Part = romapp.RomModelPartUtility(computing_model_part, restricted_residual_parameters, KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme())
        self.HR_restricted_residual_elements = self.ResidualUtilityObject_Restricted_Residuals.GetElementListFromNode(computing_model_part)
        with open('RomParameters_Residuals.json') as s:
            HR_data_residuals = json.load(s) # Load SVD of Residuals
            HR_modes_size = HR_data_residuals["rom_settings"]["number_of_rom_dofs"]
            HR_dofs_per_element = HR_data_residuals["rom_settings"]["number_of_dofs_per_element"]
            HR_elements_size = len(HR_data_residuals["elemental_modes"])
            total_restricted_dofs = int(HR_dofs_per_element*HR_elements_size)
            self.U_RESIDUAL = np.zeros((total_restricted_dofs,HR_modes_size))
            restricted_residual_elements = []
            counter_in = 0
            for key in HR_data_residuals["elemental_modes"].keys():
                counter_fin = counter_in+HR_dofs_per_element
                self.U_RESIDUAL[counter_in:counter_fin,:] = np.array(HR_data_residuals["elemental_modes"][key])
                restricted_residual_elements.append(int(key))
                counter_in = counter_fin
            self.restricted_residual_elements_index = FindDuplicateIndex(restricted_residual_elements,self.HR_restricted_residual_elements) #No need to be saved "self"
            HR_total_restricted_dofs = int(HR_dofs_per_element*len(self.HR_restricted_residual_elements))
            self.HR_U_RESIDUAL = np.zeros((HR_total_restricted_dofs,HR_modes_size))
            counter_in = 0
            for i in self.restricted_residual_elements_index:
                counter_fin = counter_in+HR_dofs_per_element
                aux_index_in = int(i*HR_dofs_per_element)
                aux_index_fin = aux_index_in + HR_dofs_per_element
                self.HR_U_RESIDUAL[counter_in:counter_fin,:] = self.U_RESIDUAL[aux_index_in:aux_index_fin,:]
                


        with open('RomParameters_Stresses.json') as s:
            HR_data_stresses = json.load(s) # Load SVD of Stresses
            HR_dofs_size = HR_data_stresses["rom_settings"]["number_of_rom_dofs"] # To initialize left singular values of Stresses
            HR_nodal_modes_size = len(HR_data_stresses["nodal_modes"]) # To initialize left singular values of Stresses
            self.U_STRESS = np.zeros((HR_nodal_modes_size,HR_dofs_size)) # Initialize left singular values of Stresses
            for key in HR_data_stresses["nodal_modes"].keys(): 
                self.U_STRESS[int(key)-1,:] = np.array(HR_data_stresses["nodal_modes"][key]) # Create numpy array of left singular values
            i = 0
        with open('test.npy', 'rb') as f:
            SnapshotMatrix_stresses = np.load(f)
            number_of_hrom_elements = len(self.HR_elements)
            self.S_hrom = np.zeros((number_of_hrom_elements,60))
            self.HR_U_STRESS = np.zeros((number_of_hrom_elements,HR_dofs_size))
            for j in self.HR_elements:
                self.HR_U_STRESS[i,:] = self.U_STRESS[j,:]
                self.S_hrom[i,:] = SnapshotMatrix_stresses[j,:]
                i += 1
            self.HR_stresses = [] # Initialize vector of stresses
            # HR_restricted_residual_parameters = KratosMultiphysics.Parameters("""
            # {
            #     "nodal_unknowns" : ["DISPLACEMENT_X","DISPLACEMENT_Y","DISPLACEMENT_Z"],
            #     "number_of_dofs" : 1
            # }""" )
            # self.ResidualUtilityObject_Restricted_Residuals = romapp.RomModelPartUtility(self._GetSolver().GetComputingModelPart().GetSubModelPart("COMPUTE_HROM"), HR_restricted_residual_parameters, KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme())
            ##########################################3  
            self.time_step=0
            ##########################################3  
        # with open('RomParameters_Residuals.json') as s:
        #     HR_data_restricted_residual = json.load(s) # Load 
        #     HR_residual_dofs_size = HR_data_restricted_residual["rom_settings"]["number_of_rom_dofs"]
        #     number_of_dofs_per_condition = HR_data_restricted_residual["rom_settings"]["number_of_dofs_per_condition"]
        #     restricted_residuals_conditions = []
        #     for key in HR_data_restricted_residual["nodal_modes"].keys():
        #         restricted_residuals_conditions.append(int(key))
        #     self.HR_restricted_residual_conditions = set(restricted_residuals_conditions).intersection(self.HR_conditions)
        #     self.HR_restricted_residual_conditions.sort()
        #     self.U_RESTRICTED_RESIDUAL = np.zeros((len(self.HR_restricted_residual_conditions)*number_of_dofs_per_condition,HR_residual_dofs_size))
        #     counter = 0
        #     for i in self.HR_restricted_residual_conditions:
        #         self.U_RESTRICTED_RESIDUAL[counter:counter+number_of_dofs_per_condition,:] = HR_data_restricted_residual["nodal_modes"][i]
        #         counter+=1



        
    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        # Added by me ####--------------------------------------------------------------------------------------------------------------------------------------
        # Save the column of von mises stresses for each solution step
        # extrapolation_parameters = KratosMultiphysics.Parameters("""
        # {
        #     "model_part_name"            : "",
        #     "echo_level"                 : 0,
        #     "average_variable"           : "NODAL_AREA",
        #     "area_average"               : true,
        #     "list_of_variables"          : ["VON_MISES_STRESS"],
        #     "extrapolate_non_historical" : true
        # }
        # """)
        # self.integration_values_extrapolation_to_nodes_process = KratosMultiphysics.IntegrationValuesExtrapolationToNodesProcess(self._GetSolver().GetComputingModelPart(), extrapolation_parameters)
        # self.integration_values_extrapolation_to_nodes_process.Execute()
        
        new_stress = []
        counter=0
        model_part = self._GetSolver().GetComputingModelPart()
        for elem in model_part.Elements:
            new_stress.append(elem.CalculateOnIntegrationPoints(KratosMultiphysics.StructuralMechanicsApplication.VON_MISES_STRESS,model_part.ProcessInfo)[0])
            # new_stress.append(self.S_hrom[counter,self.time_step].tolist())
            counter+=1
        self.HR_stresses.append(new_stress) 
        self.time_step+=1  
        HR_restricted_residual = self.ResidualUtilityObject_Main_Part.GetNonProjectedResidualsFromElementList(self.HR_restricted_residual_elements) 
        self.HR_residuals.append(HR_restricted_residual)
        ##########################################3  
        # ############################## EN REVISION #############################
        # with open('SnapshotMatrix_residuals.npy', 'rb') as f:
        #     SnapshotMatrix_residuals = np.load(f) 
        # #U_RESIDUAL,_,_,_= RandomizedSingularValueDecomposition().Calculate(SnapshotMatrix_residuals)#####
        # U_RESIDUAL,_,_ = rSVD(SnapshotMatrix_residuals,10,1,0)
        # i=0
        # HR_model_part = self._GetSolver().GetComputingModelPart().GetSubModelPart("VISUALIZE_HROM")
        # number_of_conditions = HR_model_part.NumberOfConditions()
        # number_of_fom_elements = np.shape(U_RESIDUAL)[0]-number_of_conditions
        # number_of_hrom_elements_and_conditions = len(self.HR_elements)+len(self.HR_conditions)
        # HR_U_RESIDUAL = np.zeros((number_of_hrom_elements_and_conditions,np.shape(U_RESIDUAL)[1]))
        # Residual_hrom = np.zeros((number_of_hrom_elements_and_conditions,np.shape(SnapshotMatrix_residuals)[1]))
        # for j in self.HR_elements:
        #     HR_U_RESIDUAL[i,:] = U_RESIDUAL[j,:]
        #     Residual_hrom[i,:] = SnapshotMatrix_residuals[j,:]
        #     i += 1
        # for j in self.HR_conditions:
        #     HR_U_RESIDUAL[i,:] = U_RESIDUAL[j+number_of_fom_elements,:]
        #     Residual_hrom[i,:] = SnapshotMatrix_residuals[j+number_of_fom_elements,:]
        #     i+=1
        # B = np.linalg.pinv(HR_U_RESIDUAL)@HR_residual
        # reaction_comparisson = U_RESIDUAL@B
        # ResidualUtilityObject_ReactionsSkin = romapp.RomResidualsUtility(HR_model_part, self.project_parameters["solver_settings"]["rom_settings"], KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme())
        # # Reactions = ResidualUtilityObject_ReactionsSkin.GetReactions(KratosMultiphysics.Matrix(SnapshotMatrix_residuals[number_of_fom_elements:,0:10]),1)
        # Reactions = ResidualUtilityObject_ReactionsSkin.GetReactionsNew(KratosMultiphysics.Matrix(reaction_comparisson[number_of_fom_elements:,:]),1)
        # counter_1 = 0
        # for cond in HR_model_part.Conditions:
        #     counter_2 = 0
        #     for node in cond.GetNodes():
        #         node.SetSolutionStepValue(KratosMultiphysics.REACTION_X,Reactions[counter_1,counter_2])
        #         counter_2 += 1
        #         node.SetSolutionStepValue(KratosMultiphysics.REACTION_Y,Reactions[counter_1,counter_2])
        #         counter_2 += 1
        #         node.SetSolutionStepValue(KratosMultiphysics.REACTION_Z,Reactions[counter_1,counter_2])
        #         counter_2 += 1
        # #################################################
            # counter_1 += 1
        ##########################################3    
        # self.HR_residuals.append(HR_restricted_residual)
        
   
    def Finalize(self):
        super().Finalize()
        ##########################################3  
        # Build the hrom projection for the residuals and compare the norm of the fom residuals
        self.HR_residuals = ArrangeSnapshotMatrix(self.HR_residuals)
        B = np.linalg.pinv(self.HR_U_RESIDUAL)@self.HR_residuals
        residual_comparisson = self.U_RESIDUAL@B
        with open('test_restricted_residuals.npy', 'rb') as f:
            SnapshotMatrix_residuals = np.load(f) 
        ## Comparisson purposes
        Residual_hrom = np.zeros((np.shape(self.HR_residuals)))
        counter = 0
        for i in self.restricted_residual_elements_index:
            aux_in = int(i*12)
            aux_fin = aux_in+12
            Residual_hrom[counter:counter+12,:] =  SnapshotMatrix_residuals[aux_in:aux_fin,:]
            counter += 12
        ################## Eliminar por completo, es solo para verificar la reconstruccion con los valors fom de los residuales
        B = np.linalg.pinv(self.HR_U_RESIDUAL)@Residual_hrom
        residual_comparisson = self.U_RESIDUAL@B
        ##################
        frobenius_norm = np.linalg.norm(self.HR_residuals-Residual_hrom)
        self.relative_error_5 = np.divide(frobenius_norm,np.linalg.norm(Residual_hrom))
        print("HROM reconstruction error of residuals on selected elements:")
        print(self.relative_error_5)
        frobenius_norm = np.linalg.norm(residual_comparisson-SnapshotMatrix_residuals)
        self.relative_error_4 = np.divide(frobenius_norm,np.linalg.norm(SnapshotMatrix_residuals)) 
        print("HROM reconstruction error of residuals:")
        print(self.relative_error_4)
        # #U_RESIDUAL,_,_,_= RandomizedSingularValueDecomposition().Calculate(SnapshotMatrix_residuals)#####
        # U_RESIDUAL,_,_ = rSVD(SnapshotMatrix_residuals,10,1,0)
        # i=0
        # HR_model_part = self._GetSolver().GetComputingModelPart().GetSubModelPart("VISUALIZE_HROM")
        # number_of_conditions = HR_model_part.NumberOfConditions()
        # number_of_fom_elements = np.shape(U_RESIDUAL)[0]-number_of_conditions
        # number_of_hrom_elements_and_conditions = len(self.HR_elements)+len(self.HR_conditions)
        # HR_U_RESIDUAL = np.zeros((number_of_hrom_elements_and_conditions,np.shape(U_RESIDUAL)[1]))
        # Residual_hrom = np.zeros((number_of_hrom_elements_and_conditions,np.shape(SnapshotMatrix_residuals)[1]))
        # for j in self.HR_elements:
        #     HR_U_RESIDUAL[i,:] = U_RESIDUAL[j,:]
        #     Residual_hrom[i,:] = SnapshotMatrix_residuals[j,:]
        #     i += 1
        # for j in self.HR_conditions:
        #     HR_U_RESIDUAL[i,:] = U_RESIDUAL[j+number_of_fom_elements,:]
        #     Residual_hrom[i,:] = SnapshotMatrix_residuals[j+number_of_fom_elements,:]
        #     i+=1
        # B = np.linalg.pinv(HR_U_RESIDUAL)@self.HR_residuals
        # reaction_comparisson = U_RESIDUAL@B
        # frobenius_norm = np.linalg.norm(self.HR_residuals-Residual_hrom)
        # self.relative_error_5 = np.divide(frobenius_norm,np.linalg.norm(Residual_hrom))
        # print("HROM reconstruction error of residuals on selected elements and conditions:")
        # print(self.relative_error_5)
        # # scenarios = 60
        # # ResidualUtilityObject_ReactionsSkin = romapp.RomResidualsUtility(HR_model_part, self.project_parameters["solver_settings"]["rom_settings"], KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme())
        # # Reactions = ResidualUtilityObject_ReactionsSkin.GetReactions(KratosMultiphysics.Matrix(reaction_comparisson[number_of_fom_elements:,:]),scenarios)
        # # Dof_array = ResidualUtilityObject_ReactionsSkin.GetDofArray()
        # # dof_list=[]
        # # for cond in HR_model_part.Conditions:
        # #     for node in cond.GetNodes():
        # #         dof_list.append(node.Id)
        # # count_repeated_nodes = KratosMultiphysics.Vector(HR_model_part.NumberOfNodes(),0)
        # # for i in dof_list:
        # #     count_repeated_nodes[i]+=1
        # frobenius_norm = np.linalg.norm(reaction_comparisson-SnapshotMatrix_residuals)
        # self.relative_error_4 = np.divide(frobenius_norm,np.linalg.norm(SnapshotMatrix_residuals)) 
        # print("HROM reconstruction error of residuals:")
        # print(self.relative_error_4)
        ##########################################3 
        # Build the hrom projection for the von mises stresses and compare the norm of the fom stresses
        self.HR_stresses = np.array(self.HR_stresses).T
        q = np.linalg.pinv(self.HR_U_STRESS)@self.HR_stresses
        stress_comparisson = self.U_STRESS@q
        with open('test.npy', 'rb') as f:
            SnapshotMatrix_stresses = np.load(f)
        frobenius_norm = np.linalg.norm(self.HR_stresses-self.S_hrom)
        self.relative_error_1 = np.divide(frobenius_norm,np.linalg.norm(self.S_hrom))
        print("HROM reconstruction error on selected elements' gauss points:")
        print(self.relative_error_1)
        frobenius_norm = np.linalg.norm(stress_comparisson-SnapshotMatrix_stresses)
        self.relative_error_2 = np.divide(frobenius_norm,np.linalg.norm(SnapshotMatrix_stresses))
        print("HROM reconstruction error on all gauss points:" )
        print(self.relative_error_2)
        P = load_npz("shape.npz")
        stress = P@stress_comparisson
        with open('test_stress_nodes.npy', 'rb') as f:
            stress_nodes = np.load(f)
        frobenius_norm = np.linalg.norm(stress-stress_nodes)
        self.relative_error_3 = np.divide(frobenius_norm,np.linalg.norm(stress_nodes))
        print("HROM reconstruction error on nodes:")
        print(self.relative_error_3)
        

        with open('test_stress.npy', 'wb') as g:
            np.save(g, stress_comparisson)
        with open('test_stress_snap.npy', 'wb') as g:
            np.save(g, self.HR_stresses)
        with open('new_stress.npy', 'wb') as g:
            np.save(g, stress)
        
        

def Run_HROM_KRATOS():
    with open("ProjectParameters_HROM.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())
    model = KratosMultiphysics.Model()
    simulation = RunHROM(model,parameters)
    simulation.Run()
    return simulation.relative_error_1, simulation.relative_error_2, simulation.relative_error_3
