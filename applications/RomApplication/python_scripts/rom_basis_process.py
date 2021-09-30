import KratosMultiphysics
import KratosMultiphysics.RomApplication as romapp
from KratosMultiphysics.RomApplication.randomized_singular_value_decomposition import RandomizedSingularValueDecomposition
import numpy as np
from scipy.sparse import csr_matrix, save_npz, lil_matrix
import json

def Factory(parameters, model):
    if not isinstance(parameters, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return RomBasisProcess(model, parameters["Parameters"])

## All the processes python should be derived from "Process"

class RomBasisProcess(KratosMultiphysics.Process):
    def __init__(self, Model, parameters):
        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters("""{
            "variables_to_create_basis"     : ["DISPLACEMENT","PK2_STRESS_VECTOR","REACTION"],
            "displacement_basis_tolerance" : 1e-9,
            "stress_basis_tolerance" : 1e-6,
            "reaction_basis_tolerance" : 1e-9,
            "restricted_model_part" : "Name_of_restricted_model_part",
            "plasticity" : false
        }""")

        # Add missing settings that the user did not provide but that
        # are necessary for this process
        parameters.ValidateAndAssignDefaults(default_settings)
        # Alternative:
        # settings.RecursivelyValidateAndAssignDefaults(default_settings)
        self.parameters = parameters
        variables_to_create_basis = self.parameters["variables_to_create_basis"].GetStringArray()
        self.stress_flag = False
        self.reaction_flag = False
        self.plasticity_flag = False
        for i in variables_to_create_basis:
            if i=="DISPLACEMENT":
                pass
            elif i=="PK2_STRESS_VECTOR":
                self.stress_flag = True
            elif i=="REACTION":
                self.reaction_flag = True
            else: 
                raise Exception("At least one variable to reconstruct is not yet implemented")
        
        if self.parameters["plasticity"].GetBool():
            self.plasticity_flag = True
        
        self.model_part = Model.GetModelPart("Structure")
        self.time_step_solution_container = []
        self.time_step_solution_container_stresses = []
        self.time_step_solution_container_restricted_residuals = []
        

        
    def ExecuteInitialize(self):
        """ This method is executed at the begining to initialize the process

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        if  self.stress_flag:
            self.BuildSimpleExtrapolationOperator() # Build a sparse extrapolation operator for stresses (from GP to nodes)
        
        if self.reaction_flag:
            restricted_model_part_name = self.parameters["restricted_model_part"].GetString()
            #### Reactions
            restricted_residual_parameters = KratosMultiphysics.Parameters("""
            {
                "nodal_unknowns" : ["REACTION_X","REACTION_Y","REACTION_Z"],
                "number_of_dofs" : 1
            }""")
            self.ResidualUtilityObject_Restricted_Part = romapp.RomModelPartUtility(self.model_part.GetSubModelPart(restricted_model_part_name), restricted_residual_parameters, KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme())
            self.ResidualUtilityObject_Main_Part = romapp.RomModelPartUtility(self.model_part, restricted_residual_parameters, KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme())
            self.restricted_residual_elements = self.ResidualUtilityObject_Restricted_Part.GetElementListFromNode(self.model_part)
            

    def ExecuteFinalizeSolutionStep(self):
        """ This method is executed in order to finalize the current step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        if self.reaction_flag:
            self.time_step_solution_container_restricted_residuals.append(self.ResidualUtilityObject_Main_Part.GetNonProjectedResidualsFromElementList(self.restricted_residual_elements))
        
        if self.stress_flag:    
            ArrayOfStresses = []
            for elem in self.model_part.Elements:
                PK2_STRESS = elem.CalculateOnIntegrationPoints(KratosMultiphysics.PK2_STRESS_VECTOR,self.model_part.ProcessInfo)
                for i in range(PK2_STRESS[0].Size()):
                    ArrayOfStresses.append(PK2_STRESS[0][i])
            self.time_step_solution_container_stresses.append(ArrayOfStresses)

        ArrayOfDisplacements = []    
        for node in self.model_part.Nodes:
            ArrayOfDisplacements.append(node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X, 0))
            ArrayOfDisplacements.append(node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y, 0))
            ArrayOfDisplacements.append(node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z, 0))
        self.time_step_solution_container.append(ArrayOfDisplacements)


    def ExecuteFinalize(self):
        """ This method is executed after the computations, at the end of the solution-loop

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        SnapshotMatrix = self.EvaluateQuantityOfInterest(self.time_step_solution_container)

        # Taking the SVD ###
        displacement_tolerance = self.parameters["displacement_basis_tolerance"].GetDouble()
        if self.plasticity_flag:
            u = self.BuildElastoPlasticBasis(SnapshotMatrix,displacement_tolerance)
        else:
            u,_,_,_= RandomizedSingularValueDecomposition().Calculate(SnapshotMatrix,displacement_tolerance)

        ### Saving the nodal basis ###  
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
        
        # Create rom parameters for stresses
        if self.stress_flag:
            SnapshotMatrix_stresses = self.EvaluateQuantityOfInterest(self.time_step_solution_container_stresses)
            stress_tolerance = self.parameters["stress_basis_tolerance"].GetDouble()
            if self.plasticity_flag:
                U = self.BuildElastoPlasticBasis(SnapshotMatrix_stresses,stress_tolerance)
            else:
                U,_,_,_= RandomizedSingularValueDecomposition().Calculate(SnapshotMatrix_stresses,stress_tolerance)# Taking the SVD ###
            basis_POD["stress_parameters"]={"rom_settings":{},"gauss_points_modes":{}}
            basis_POD["stress_parameters"]["rom_settings"]["gauss_points_unknowns"] = ["PK2_STRESS_VECTOR"]
            basis_POD["stress_parameters"]["rom_settings"]["number_of_rom_modes"] = np.shape(U)[1]
            Dimensions = 6
            N_elem=np.shape(U)[0]/Dimensions
            N_elem = int(N_elem)
            elem_Id=np.linspace(1,N_elem,N_elem)
            i = 0
            for j in range (0,N_elem):
                basis_POD["stress_parameters"]["gauss_points_modes"][int(elem_Id[j])] = (U[i:i+Dimensions].tolist())
                i=i+Dimensions

        # Create rom parameters for restricted residuals and reactions
        if self.reaction_flag:
            SnapshotMatrix_restricted_residuals = self.EvaluateQuantityOfInterest(self.time_step_solution_container_restricted_residuals)
            residual_tolerance = self.parameters["reaction_basis_tolerance"].GetDouble()
            if self.plasticity_flag:
                U_residuals = self.BuildElastoPlasticBasis(SnapshotMatrix_restricted_residuals,residual_tolerance)
            else:
                U_residuals,_,_,_= RandomizedSingularValueDecomposition().Calculate(SnapshotMatrix_restricted_residuals,residual_tolerance)# Taking the SVD ###
            basis_POD["residual_parameters"]={"rom_settings":{},"elemental_modes":{}}
            basis_POD["residual_parameters"]["rom_settings"]["nodal_unknowns"] = ["REACTION_X","REACTION_Y","REACTION_Z"]
            basis_POD["residual_parameters"]["rom_settings"]["number_of_rom_dofs"] = np.shape(U_residuals)[1]
            basis_POD["residual_parameters"]["rom_settings"]["restricted_elements"] = self.restricted_residual_elements
            N_cond= len(self.restricted_residual_elements)
            Dimensions = int(np.shape(U_residuals)[0]/N_cond)
            basis_POD["residual_parameters"]["rom_settings"]["number_of_dofs_per_element"] = Dimensions
            i = 0
            for j in self.restricted_residual_elements:
                basis_POD["residual_parameters"]["elemental_modes"][int(j)] = (U_residuals[i:i+Dimensions].tolist())
                i=i+Dimensions

            # Reaction Assembling paraameters
            restricted_residual_nodes = self.ResidualUtilityObject_Restricted_Part.GetNodeList()
            aux_list_nodes, aux_list_elem = self.CreateAuxiliarListForReactions(restricted_residual_nodes) # Builds list of lists to assemble the residuals to restricted surface's nodes
            basis_POD["reaction_parameters"]={"Nodal_Unknowns":{},"Element_Dofs":{},"List_Nodes":{},"List_Nodal_Elements":{},"List_Elemental_Nodes":{}}
            basis_POD["reaction_parameters"]["Nodal_Unknowns"] = ["REACTION_X","REACTION_Y","REACTION_Z"]
            basis_POD["reaction_parameters"]["Nodes_per_element"] = self.model_part.GetElement(1).GetGeometry().PointsNumber()
            basis_POD["reaction_parameters"]["List_Nodes"] = restricted_residual_nodes
            basis_POD["reaction_parameters"]["List_Nodal_Elements"] = aux_list_nodes
            basis_POD["reaction_parameters"]["List_Elemental_Nodes"] = aux_list_elem

        with open('RomParameters.json', 'w') as f:
            json.dump(basis_POD,f, indent=2)

        print('\n\nROM basis printed in json format\n\n')

    def CreateAuxiliarListForReactions(self, restricted_residual_nodes):
        # Builds list of lists to assemble the residuals to restricted surface's nodes
        aux_list_nodes = [[] for i in range(len(restricted_residual_nodes))] # List of lists which contains the neighbour elements per node
        aux_list_elem = [[] for i in range(len(self.restricted_residual_elements))] # List of lists which contains the neighbour nodes per element
        counter=0
        for i in self.restricted_residual_elements:
            i_elem = self.model_part.GetElement(i)
            for node in i_elem.GetNodes():
                i_node = node.Id
                aux_list_elem[counter].append(i_node)
                for j in range(len(restricted_residual_nodes)):
                    if i_node==restricted_residual_nodes[j]:
                        aux_list_nodes[j].append(counter)
            counter+=1
        return aux_list_nodes, aux_list_elem
    
    def BuildExtrapolationOperator(self):
        # Build a sparse extrapolation operator for stresses (from GP to nodes)
        #############
        P = lil_matrix((self.model_part.NumberOfNodes(),self.model_part.NumberOfElements()))
        num_of_nodes = P.shape[0]
        P_diag = np.zeros(num_of_nodes)
        for elem in self.model_part.Elements:
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
        P_diag = np.diag(1/P_diag)
        P_diag = lil_matrix(P_diag)
        P = P_diag@P
        save_npz("shape.npz", P)
        #############

    def BuildSimpleExtrapolationOperator(self):
        P=lil_matrix((self.model_part.NumberOfNodes(),self.model_part.NumberOfElements()))
        num_of_nodes = P.shape[0]
        P_diag=np.zeros(num_of_nodes)
        for elem in self.model_part.Elements:
            elem_geom = elem.GetGeometry()
            elem_area = elem_geom.Area()
            nodes_in_element = elem.GetNodes()
            for node in nodes_in_element:
                P_diag[node.Id-1] += elem_area
                P[node.Id-1,elem.Id-1] = elem_area
        P_diag = np.diag(1/P_diag)
        P_diag = lil_matrix(P_diag)
        P = P_diag@P
        save_npz("shape.npz", P)
    
    def EvaluateQuantityOfInterest(self,time_step_solution_container):
       ##############################################################################################
       #     Functions evaluating the QoI of the problem: Snapshot matrix of every time step        #
       #                                                                                            #
       ##############################################################################################
        SnapshotMatrix = np.zeros((len(time_step_solution_container[0]), len(time_step_solution_container)))
        for i in range(len(time_step_solution_container)):
            Snapshot_i= np.array(time_step_solution_container[i])
            SnapshotMatrix[:,i] = Snapshot_i.transpose()
        return SnapshotMatrix

    def rSVD(self,X,r,q,p):
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
    
    def BuildElastoPlasticBasis(self,SnapshotMatrix,basis_tolerance):
        training_size = np.shape(SnapshotMatrix)[1]
        plasticity_basis_tolerance = 0.5
        u_shape = []
        for i in range(1,training_size):
            u,_,_,_= RandomizedSingularValueDecomposition().Calculate(SnapshotMatrix[:,0:i],basis_tolerance)
            u_shape.append(np.shape(u)[1])
            try:
                gradient = np.gradient(u_shape)
                if gradient[-1]>plasticity_basis_tolerance:
                    elastic_step_range = i
                    break
            except:
                pass
        u,_,_,_= RandomizedSingularValueDecomposition().Calculate(SnapshotMatrix[:,0:elastic_step_range],basis_tolerance)
        for j in range(3):
            aux = u.T@SnapshotMatrix
            aux = u@aux
            SnapshotMatrix = SnapshotMatrix-aux
        u_p,_,_,_= RandomizedSingularValueDecomposition().Calculate(SnapshotMatrix[:,:],basis_tolerance)
        u = np.concatenate((u,u_p),axis=1)
        Q, R = np.linalg.qr(u,mode='reduced')
        return Q

    def DMD(self,X,Xprime,r):
        U,Sigma,VT = np.linalg.svd(X,full_matrices=0) # Step 1
        Ur = U[:,:r]
        Sigmar = np.diag(Sigma[:r])
        VTr = VT[:r,:]
        Atilde = np.linalg.solve(Sigmar.T,(Ur.T @ Xprime @ VTr.T).T).T # Step 2
        Lambda, W = np.linalg.eig(Atilde) # Step 3
        Lambda = np.diag(Lambda)
        
        Phi = Xprime @ np.linalg.solve(Sigmar.T,VTr).T @ W # Step 4
        alpha1 = Sigmar @ VTr[:,0]
        b = np.linalg.solve(W @ Lambda,alpha1)
        return np.real(Phi)