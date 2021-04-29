import KratosMultiphysics
import KratosMultiphysics.RomApplication as romapp
from KratosMultiphysics.RomApplication.randomized_singular_value_decomposition import RandomizedSingularValueDecomposition
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis
from KratosMultiphysics.RomApplication.structural_mechanics_analysis_rom import StructuralMechanicsAnalysisROM
from matplotlib import pyplot as plt
import numpy as np
import json
from scipy.sparse import csr_matrix, save_npz, lil_matrix
from MainKratos_HROM import Run_HROM_KRATOS
from postprocessing import RunPostProcess

import pdb



"""
For user-scripting it is intended that a new class is derived
from StructuralMechanicsAnalysis to do modifications
"""

class StructuralMechanicsAnalysisSavingData(StructuralMechanicsAnalysis):

    def __init__(self,model,project_parameters):
        super().__init__(model,project_parameters)
        self.time_step_solution_container = []
        #Added by me
        self.time_step_solution_container_stresses = []
        self.time_step_solution_container_stresses_nodes = []
        self.time_step_solution_container_restricted_residuals = []

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        ArrayOfDisplacements = []
        
        # Initialize array of stresses and extrapolate von mises to nodes.
        ArrayOfStresses = []
        ArrayOfStressesNodes = []
        self.time_step_solution_container_restricted_residuals.append(self.ResidualUtilityObject_Main_Part.GetNonProjectedResidualsFromElementList(self.restricted_residual_elements))
        model_part = self._GetSolver().GetComputingModelPart()
        for elem in model_part.Elements:
            PK2_STRESS = elem.CalculateOnIntegrationPoints(KratosMultiphysics.PK2_STRESS_VECTOR,model_part.ProcessInfo)
            for i in range(PK2_STRESS[0].Size()):
                ArrayOfStresses.append(PK2_STRESS[0][i])
        for node in model_part.Nodes:
            ArrayOfDisplacements.append(node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X, 0))
            ArrayOfDisplacements.append(node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y, 0))
            ArrayOfDisplacements.append(node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z, 0))
            ArrayOfStressesNodes.append(node.GetValue(KratosMultiphysics.PK2_STRESS_VECTOR))
        self.time_step_solution_container.append(ArrayOfDisplacements)
        self.time_step_solution_container_stresses.append(ArrayOfStresses)
        self.time_step_solution_container_stresses_nodes.append(ArrayOfStressesNodes)

    def ModifyInitialGeometry(self):
        super().ModifyInitialGeometry()
        model_part = self._GetSolver().GetComputingModelPart()
        # P = np.zeros((model_part.NumberOfNodes(),model_part.NumberOfElements()))

        # Build a sparse extrapolation operator for stresses (from GP to nodes)
        #############
        P = lil_matrix((model_part.NumberOfNodes(),model_part.NumberOfElements()))
        num_of_nodes = P.shape[0]
        P_diag = np.zeros(num_of_nodes)
        for elem in model_part.Elements:
            elem_geom = elem.GetGeometry()
            shape_functions_evaluated_in_gp = np.array(elem_geom.ShapeFunctionsValues())
            elem_area = elem_geom.Area()
            sum_shape_functions_evaluated_in_gp = shape_functions_evaluated_in_gp.sum(axis=0)
            nodes_in_element = elem.GetNodes()
            for node in nodes_in_element:
                P_diag[node.Id-1] += elem_area
            for i in range(shape_functions_evaluated_in_gp.shape[0]):
                for j in range(shape_functions_evaluated_in_gp.shape[1]):
                    P[nodes_in_element[j].Id-1,elem.Id-1] = shape_functions_evaluated_in_gp[i,j]*elem_area/sum_shape_functions_evaluated_in_gp[j]
        # P = csr_matrix(P)
        P_diag = np.diag(1/P_diag)
        P_diag = csr_matrix(P_diag)
        P = P_diag@P
        P = csr_matrix(P)
        save_npz("shape.npz", P)
        #############

        #### Reactions
        #############
        restricted_residual_parameters = KratosMultiphysics.Parameters("""
        {
            "nodal_unknowns" : ["DISPLACEMENT_X","DISPLACEMENT_Y","DISPLACEMENT_Z"],
            "number_of_dofs" : 1
        }""" )
        self.ResidualUtilityObject_Restricted_Part = romapp.RomModelPartUtility(model_part.GetSubModelPart("GENERIC_Interface_frame_1"), restricted_residual_parameters, KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme())
        self.ResidualUtilityObject_Main_Part = romapp.RomModelPartUtility(model_part, restricted_residual_parameters, KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme())
        self.restricted_residual_elements = self.ResidualUtilityObject_Restricted_Part.GetElementListFromNode(model_part)
        restricted_residual_nodes = self.ResidualUtilityObject_Restricted_Part.GetNodeList()
        # Build list of lists to assemble the residuals to restricted surface's nodes
        aux_list_nodes = [[] for i in range(len(restricted_residual_nodes))] # List of lists which contains the neighbour elements per node
        aux_list_elem = [[] for i in range(len(self.restricted_residual_elements))] # List of lists which contains the neighbour nodes per element
        counter=0
        for i in self.restricted_residual_elements:
            i_elem = model_part.GetElement(i)
            for node in i_elem.GetNodes():
                i_node = node.Id
                aux_list_elem[counter].append(i_node)
                for j in range(len(restricted_residual_nodes)):
                    if i_node==restricted_residual_nodes[j]:
                        aux_list_nodes[j].append(counter)
            counter+=1
        reaction_assembling={"Nodal_Unknowns":{},"Element_Dofs":{},"List_Nodes":{},"List_Nodal_Elements":{},"List_Elemental_Nodes":{}}
        reaction_assembling["Nodal_Unknowns"] = ["REACTION_X","REACTION_Y","REACTION_Z"]
        reaction_assembling["Nodes_per_element"] = model_part.GetElement(1).GetGeometry().PointsNumber()
        reaction_assembling["List_Nodes"] = restricted_residual_nodes
        reaction_assembling["List_Nodal_Elements"] = aux_list_nodes
        reaction_assembling["List_Elemental_Nodes"] = aux_list_elem
        
        with open('Reaction_parameters.json', 'w') as f:
            json.dump(reaction_assembling,f, indent=2)
        #############


    def EvaluateQuantityOfInterest(self):
       ##############################################################################################
       #     Functions evaluating the QoI of the problem: Snapshot matrix of every time step        #
       #                                                                                            #
       ##############################################################################################
        SnapshotMatrix = np.zeros((len(self.time_step_solution_container[0]), len(self.time_step_solution_container)))
        SnapshotMatrix_stresses = np.zeros((len(self.time_step_solution_container_stresses[0]), len(self.time_step_solution_container_stresses)))
        SnapshotMatrix_stresses_nodes = np.zeros((len(self.time_step_solution_container_stresses_nodes[0]), len(self.time_step_solution_container_stresses_nodes)))
        SnapshotMatrix_restricted_residuals = np.zeros((len(self.time_step_solution_container_restricted_residuals[0]), len(self.time_step_solution_container_restricted_residuals)))
        for i in range(len(self.time_step_solution_container)):
            Snapshot_i= np.array(self.time_step_solution_container[i])
            SnapshotMatrix[:,i] = Snapshot_i.transpose()
            Snapshot_i_stresses = np.array(self.time_step_solution_container_stresses[i])
            SnapshotMatrix_stresses[:,i] = Snapshot_i_stresses.transpose()
            # Snapshot_i_stresses_nodes= np.array(self.time_step_solution_container_stresses_nodes[i])
            # SnapshotMatrix_stresses_nodes[:,i] = Snapshot_i_stresses_nodes.transpose()
            Snapshot_i_restricted_residuals= np.array(self.time_step_solution_container_restricted_residuals[i])
            SnapshotMatrix_restricted_residuals[:,i] = Snapshot_i_restricted_residuals.transpose()
        return SnapshotMatrix, SnapshotMatrix_stresses, SnapshotMatrix_stresses_nodes, SnapshotMatrix_restricted_residuals

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







if __name__ == "__main__":
    # tolerance_1 = [0,1e-9,1e-7]
    # tolerance_2 = [0,1e-6,1e-4]
    # tol_1 = [10,9,8,7,6,5,4,3,2]
    # tol_2 = [10,9,8,7,6,5,4,3,2]

    relative_error_1 = []
    relative_error_2 = []
    relative_error_3 = []
    index=1
    modes_displacement = []
    modes_stresses = []
    for i in range(1):
        ##############################################################################################
        #                                           TRAIN ROM                                        #
        ##############################################################################################
        with open("ProjectParameters.json",'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())


        model = KratosMultiphysics.Model()
        simulation = StructuralMechanicsAnalysisSavingData(model,parameters)
        simulation.Run()
        SnapshotMatrix, SnapshotMatrix_stresses, SnapshotMatrix_stresses_nodes, SnapshotMatrix_restricted_residuals = simulation.EvaluateQuantityOfInterest()
        
        # # Load npy file to avoid running the FOM again
        # with open('test.npy', 'rb') as f:
        #     SnapshotMatrix_stresses = np.load(f)
        # with open('test_disp.npy', 'rb') as f:
        #     SnapshotMatrix = np.load(f)
        # with open('test_restricted_residuals.npy', 'rb') as f:
        #     SnapshotMatrix_restricted_residuals = np.load(f)
        # # Save as npy file to check the norm
        with open('test.npy', 'wb') as f:
            np.save(f, SnapshotMatrix_stresses)
        with open('test_disp.npy', 'wb') as f:
            np.save(f, SnapshotMatrix)
        with open('test_stress_nodes.npy', 'wb') as f:
            np.save(f, SnapshotMatrix_stresses_nodes)
        with open('test_restricted_residuals.npy', 'wb') as f:
            np.save(f, SnapshotMatrix_restricted_residuals)
        

        # Taking the SVD ###
        u,s,_,_= RandomizedSingularValueDecomposition().Calculate(SnapshotMatrix,1e-9)#,tolerance_1[i])
        U,S,_,_= RandomizedSingularValueDecomposition().Calculate(SnapshotMatrix_stresses,1e-5)#,tolerance_2[i])
        U_residuals,S_residuals,_,_= RandomizedSingularValueDecomposition().Calculate(SnapshotMatrix_restricted_residuals,1e-9)
        
        # u,s,_ = rSVD(SnapshotMatrix,20,1,0)
        # U,S,_ = rSVD(SnapshotMatrix_stresses,10,1,0)
        # U_residuals,S_residuals, _ = rSVD(SnapshotMatrix_restricted_residuals,10,1,0)

        ### Plotting singular values  ###
        # plt.plot( range(0,len(s)), np.log(s), marker='o', markerfacecolor='blue', markersize=12, color='skyblue', linewidth=4)
        # plt.title('Singular Values')
        # plt.show()

        
        # # Plot results of singular values for stresses 
        # fig = plt.figure()
        # ax1 = fig.add_subplot(121)
        # ax1.semilogy(S,'-o',color='k')
        # ax2 = fig.add_subplot(122)
        # ax2.plot(np.cumsum(S)/np.sum(S),'-o',color='k')
        # plt.show()
        


        ### Saving the nodal basis ###  (Need to make this more robust, hard coded here)
        basis_POD={"rom_settings":{},"nodal_modes":{}}
        basis_POD["rom_settings"]["nodal_unknowns"] = ["DISPLACEMENT_X","DISPLACEMENT_Y","DISPLACEMENT_Z"]
        basis_POD["rom_settings"]["number_of_rom_dofs"] = np.shape(u)[1]
        Dimensions = len(basis_POD["rom_settings"]["nodal_unknowns"])
        N_nodes=np.shape(u)[0]/Dimensions
        N_nodes = int(N_nodes)
        node_Id=np.linspace(1,N_nodes,N_nodes)
        i = 0
        for j in range (0,N_nodes):
            basis_POD["nodal_modes"][int(node_Id[j])] = (u[i:i+Dimensions].tolist())
            i=i+Dimensions

        with open('RomParameters.json', 'w') as f:
            json.dump(basis_POD,f, indent=2)

        print('\n\nDisplacement basis printed in json format\n\n')

        
        # Create rom parameters for stresses
        basis_POD_stresses={"rom_settings":{},"gauss_points_modes":{}}
        basis_POD_stresses["rom_settings"]["gauss_points_unknowns"] = ["PK2_STRESS_VECTOR"]
        basis_POD_stresses["rom_settings"]["number_of_rom_modes"] = np.shape(U)[1]
        Dimensions = 6
        N_elem=np.shape(U)[0]/Dimensions
        N_elem = int(N_elem)
        elem_Id=np.linspace(1,N_elem,N_elem)
        i = 0
        for j in range (0,N_elem):
            basis_POD_stresses["gauss_points_modes"][int(elem_Id[j])] = (U[i:i+Dimensions].tolist())
            i=i+Dimensions

        with open('RomParameters_Stresses.json', 'w') as w:
            json.dump(basis_POD_stresses,w, indent=2)

        print('\n\nStress basis printed in json format\n\n')

        # Create rom parameters for restricted residuals
        basis_POD={"rom_settings":{},"elemental_modes":{}}
        basis_POD["rom_settings"]["nodal_unknowns"] = ["REACTION_X","REACTION_Y","REACTION_Z"]
        basis_POD["rom_settings"]["number_of_rom_dofs"] = np.shape(U_residuals)[1]
        basis_POD["rom_settings"]["restricted_elements"] = simulation.restricted_residual_elements
        N_cond= len(simulation.restricted_residual_elements)
        Dimensions = int(np.shape(U_residuals)[0]/N_cond)
        basis_POD["rom_settings"]["number_of_dofs_per_element"] = Dimensions
        i = 0
        for j in simulation.restricted_residual_elements:
            basis_POD["elemental_modes"][int(j)] = (U_residuals[i:i+Dimensions].tolist())
            i=i+Dimensions

        with open('RomParameters_Residuals.json', 'w') as f:
            json.dump(basis_POD,f, indent=2)

        print('\n\nReaction basis printed in json format\n\n')
        


        ##############################################################################################
        #                                          TRAIN HROM                                        #
        ##############################################################################################
        with open("ProjectParameters_ROM.json",'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())
        model = KratosMultiphysics.Model()
        simulation = StructuralMechanicsAnalysisROM(model,parameters,"EmpiricalCubature")

        simulation.Run()

        mode_disp = np.shape(u)[1]
        mode_stress = np.shape(U)[1]

        modes_displacement.append(mode_disp)
        modes_stresses.append(mode_stress)
        rel_error_1, rel_error_2 = Run_HROM_KRATOS()
        relative_error_1.append(rel_error_1)
        relative_error_2.append(rel_error_2)
        # relative_error_3.append(rel_error_3)
        
        # RunPostProcess(index)


        index += 1

    plt.plot( modes_stresses, relative_error_1, marker='o', markerfacecolor='blue', markersize=12, color='skyblue', linewidth=4)
    plt.title("HROM reconstruction error on selected elements' gauss points")
    # plt.xlim(max(modes_stresses), min(modes_stresses))
    plt.show()

    plt.plot( modes_stresses, relative_error_2, marker='o', markerfacecolor='blue', markersize=12, color='skyblue', linewidth=4)
    plt.title('HROM reconstruction error on all gauss points')
    # plt.xlim(max(modes_stresses), min(modes_stresses))
    plt.show()

    # plt.plot( modes_stresses, relative_error_3, marker='o', markerfacecolor='blue', markersize=12, color='skyblue', linewidth=4)
    # plt.title('HROM reconstruction error on nodes')
    # # plt.xlim(max(modes_stresses), min(modes_stresses))
    # plt.show()

    with open('relative_error_1.npy', 'wb') as g:
            np.save(g, np.array(relative_error_1))
    with open('relative_error_2.npy', 'wb') as g:
            np.save(g, np.array(relative_error_2))
    # with open('relative_error_3.npy', 'wb') as g:
    #         np.save(g, np.array(relative_error_3))