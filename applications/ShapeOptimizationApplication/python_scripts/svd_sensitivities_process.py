import KratosMultiphysics
import KratosMultiphysics.ShapeOptimizationApplication as ShapeOptimizationApplication

import numpy as np

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return SVDSensitivitiesProcess(Model, settings["Parameters"])

class SVDSensitivitiesProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings):
        KratosMultiphysics.Process.__init__(self)
        
        self.model_part = Model[settings["model_part_name"].GetString()]
        self.svd_file_name = settings["file_name"].GetString()
        #ShapeOptimizationApplication.SVDSensitivitiesProcess.__init__(self, self.model_part, settings)
        
        
    # --------------------------------------------------------------------------
    def ExecuteInitialize(self):
        self.svd_file = open(self.svd_file_name, 'w')
        self.svd_file.write("#Output svd results\n")
        #self.__getAdjointDisplacements()
        
    # --------------------------------------------------------------------------
    def ExecuteFinalizeSolutionStep(self):
        #self.__test_np_svd()
        #self.__svdAdjointDisplacements()
        self.__svdElementalSensitivities()
        
    # --------------------------------------------------------------------------
    # svd Elemental Sensitivities
    def __svdElementalSensitivities(self):
        
        numElement = len(self.model_part.Elements)
        
        elemsensi_1 = np.empty(numElement)
        elemsensi_2 = np.empty(numElement)
        
        for elem in self.model_part.Elements:
            elemsensi_1[elem.Id - 1] = elem.GetValue(ShapeOptimizationApplication.IY_SENSITIVITY_1)
            elemsensi_2[elem.Id - 1] = elem.GetValue(ShapeOptimizationApplication.IY_SENSITIVITY_2)
            
        #print(elemsensi_1)
        #print(elemsensi_2)
        
        # SVD for IY_SENSITIVITY_1
        sensi_matrix = np.stack((elemsensi_1,elemsensi_2))
        
        U,s,V = np.linalg.svd(sensi_matrix, full_matrices = True)
        
        print("singular values are:")
        print(s)
        print("objective modes are:")
        print(U)
        print("input modes are:")
        print(V)
        
        # set input modes to the elements
        for elem in self.model_part.Elements:
            elem.SetValue(ShapeOptimizationApplication.IY_SENSITIVITY_1, V[(elem.Id - 1),0])
            elem.SetValue(ShapeOptimizationApplication.IY_SENSITIVITY_2, V[(elem.Id - 1),1])
            #print(elem)
    # --------------------------------------------------------------------------
    # svd adjoint displacements 
    def __svdAdjointDisplacements(self):
        
        numNode = len(self.model_part.Nodes)
        
        adjoint_disp_1 = np.empty(numNode*3)
        adjoint_disp_2 = np.empty(numNode*3)
        
        for node in self.model_part.Nodes:
            adjoint_disp_1[(node.Id - 1)*3] = node.GetSolutionStepValue(ShapeOptimizationApplication.ADJOINT_DISPLACEMENT_1_X,0)
            adjoint_disp_1[(node.Id - 1)*3 + 1] = node.GetSolutionStepValue(ShapeOptimizationApplication.ADJOINT_DISPLACEMENT_1_Y,0)
            adjoint_disp_1[(node.Id - 1)*3 + 2] = node.GetSolutionStepValue(ShapeOptimizationApplication.ADJOINT_DISPLACEMENT_1_Z,0)
            adjoint_disp_2[(node.Id - 1)*3] = node.GetSolutionStepValue(ShapeOptimizationApplication.ADJOINT_DISPLACEMENT_2_X,0)
            adjoint_disp_2[(node.Id - 1)*3 + 1] = node.GetSolutionStepValue(ShapeOptimizationApplication.ADJOINT_DISPLACEMENT_2_Y,0)
            adjoint_disp_2[(node.Id - 1)*3 + 2] = node.GetSolutionStepValue(ShapeOptimizationApplication.ADJOINT_DISPLACEMENT_2_Z,0)
        
        print(adjoint_disp_1)
        # SVD for displacement_x
        adjoint_disp_matrix = np.stack((adjoint_disp_1,adjoint_disp_2))
        
        U,s,V = np.linalg.svd(adjoint_disp_matrix, full_matrices = True)
        S = np.zeros((2,2))
        S[:2,:2] = np.diag(s)
        
        print("singular values are:")
        print(s)
        print("objective modes are:")
        print(U)
        
        # set adjoint displacement to the nodes
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(ShapeOptimizationApplication.ADJOINT_DISPLACEMENT_1_X,0,V[(node.Id - 1)*3,0])
            node.SetSolutionStepValue(ShapeOptimizationApplication.ADJOINT_DISPLACEMENT_1_Y,0,V[(node.Id - 1)*3 + 1,0])
            node.SetSolutionStepValue(ShapeOptimizationApplication.ADJOINT_DISPLACEMENT_1_Z,0,V[(node.Id - 1)*3 + 2,0])
            node.SetSolutionStepValue(ShapeOptimizationApplication.ADJOINT_DISPLACEMENT_2_X,0,V[(node.Id - 1)*3,1])
            node.SetSolutionStepValue(ShapeOptimizationApplication.ADJOINT_DISPLACEMENT_2_Y,0,V[(node.Id - 1)*3 + 1,1])
            node.SetSolutionStepValue(ShapeOptimizationApplication.ADJOINT_DISPLACEMENT_2_Z,0,V[(node.Id - 1)*3 + 2,1])      
            #node.SetSolutionStepValue(ShapeOptimizationApplication.ADJOINT_DISPLACEMENT_2_X,0,1.0)
            #node.SetSolutionStepValue(ShapeOptimizationApplication.ADJOINT_DISPLACEMENT_2_Y,0,2.0)
            #node.SetSolutionStepValue(ShapeOptimizationApplication.ADJOINT_DISPLACEMENT_2_Z,0,3.0)                  
        
    # --------------------------------------------------------------------------
    # get adjoint displacements and store them in matrix
    def __getAdjointDisplacements_seperated(self):
        #numNodes = self.model_part.NumberOfNodes()
        #print(self.model_part.Nodes)
        adjoint_disp_1_x = np.empty(794)
        adjoint_disp_1_y = np.empty(794)
        adjoint_disp_1_z = np.empty(794)
        adjoint_disp_2_x = np.empty(794)
        adjoint_disp_2_y = np.empty(794)
        adjoint_disp_2_z = np.empty(794)
        for node in self.model_part.Nodes:
            adjoint_disp_1_x[node.Id-1] = node.GetSolutionStepValue(ShapeOptimizationApplication.ADJOINT_DISPLACEMENT_1_X,0)
            adjoint_disp_1_y[node.Id-1] = node.GetSolutionStepValue(ShapeOptimizationApplication.ADJOINT_DISPLACEMENT_1_Y,0)
            adjoint_disp_1_z[node.Id-1] = node.GetSolutionStepValue(ShapeOptimizationApplication.ADJOINT_DISPLACEMENT_1_Z,0)
            adjoint_disp_2_x[node.Id-1] = node.GetSolutionStepValue(ShapeOptimizationApplication.ADJOINT_DISPLACEMENT_2_X,0)
            adjoint_disp_2_y[node.Id-1] = node.GetSolutionStepValue(ShapeOptimizationApplication.ADJOINT_DISPLACEMENT_2_Y,0)
            adjoint_disp_2_z[node.Id-1] = node.GetSolutionStepValue(ShapeOptimizationApplication.ADJOINT_DISPLACEMENT_2_Z,0)
            
        print(adjoint_disp_1_z)
        
        #adjoint_disp_matrix = np.zeros([794,6])
        #adjoint_disp_matrix = np.stack((adjoint_disp_1_x,adjoint_disp_1_y,adjoint_disp_1_z,adjoint_disp_2_x,adjoint_disp_2_y,adjoint_disp_2_z))

        # SVD for displacement_x
        adjoint_disp_x_matrix = np.stack((adjoint_disp_1_x,adjoint_disp_2_x))
        
        U_x,s_x,V_x = np.linalg.svd(adjoint_disp_x_matrix, full_matrices = True)
        S_x = np.zeros((2,2))
        S_x[:2,:2] = np.diag(s_x)
        
        # SVD for displacement_y
        adjoint_disp_y_matrix = np.stack((adjoint_disp_1_y,adjoint_disp_2_y))
        
        U_y,s_y,V_y = np.linalg.svd(adjoint_disp_y_matrix, full_matrices = True)
        S_y = np.zeros((2,2))
        S_y[:2,:2] = np.diag(s_y)
        
        # SVD for displacement_z
        adjoint_disp_z_matrix = np.stack((adjoint_disp_1_z,adjoint_disp_2_z))
        
        U_z,s_z,V_z = np.linalg.svd(adjoint_disp_z_matrix, full_matrices = True)
        S_z = np.zeros((2,2))
        S_z[:2,:2] = np.diag(s_z)
        
        print(S_x)
        print(S_y)
        print(S_z)
        print("Objective modes:\n")
        print(U_x)
        print(U_y)
        print(U_z)
        #print(S[0,0])
        #print(S[1,1])
        #print(S[2,2])
        #print(S[3,3])
        #print(S[4,4])
        #print(S[5,5])
        #print(s)
     
        # Extract Input vectors from V matrices
        
        
    
    # --------------------------------------------------------------------------    
    #def __storeResultOfSensitivityAnalysisOnNodes(self):
    
    def __column(matrix, i):
        return [row[i] for row in matrix]
    
    
    
    
    
    
    # --------------------------------------------------------------------------    
    def __test_np_svd(self):
        a = np.random.randn(3, 5) #+ 1j*np.random.randn(9, 6)
        U, s, V = np.linalg.svd(a, full_matrices=True)
        #U.shape, V.shape, s.shape((9, 9), (6, 6), (6,))
        S = np.zeros((3, 5))
        S[:3, :3] = np.diag(s)
        np.allclose(a, np.dot(U, np.dot(S, V)))
        
        print(U)
        print(S)
        print(V)
        #self.svd_file.write('{0:12.5e} {1:22.15e} {2:22.15e}\n'.format(S[1],S[2],S[3]))
        #self.svd_file.write.flush()
        
    # --------------------------------------------------------------------------
    def ExecuteFinalize(self):
        self.svd_file.close()

