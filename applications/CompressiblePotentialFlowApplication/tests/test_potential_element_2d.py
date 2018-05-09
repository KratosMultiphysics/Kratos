from KratosMultiphysics import *
from KratosMultiphysics.CompressiblePotentialFlowApplication import *
import KratosMultiphysics.KratosUnittest as KratosUnittest
import random
from numpy.linalg import inv
from numpy.linalg import det
#from sympy import *

class TestPotentialElement2D(KratosUnittest.TestCase):

    def setUp(self):
        # create test model part
        self.model_part = ModelPart("test")
        self.model_part.AddNodalSolutionStepVariable(POSITIVE_FACE_PRESSURE)
        self.model_part.AddNodalSolutionStepVariable(NEGATIVE_FACE_PRESSURE)
        self.model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        self.model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
        self.model_part.CreateNewNode(3, 1.0, 1.0, 0.0)
        prop = self.model_part.GetProperties()[0]
        self.model_part.CreateNewElement("CompressiblePotentialFlowElement2D3N", 1, [1, 2, 3], prop)

        self.potential_element = self.model_part.GetElement(1)

        self._assign_solution_step_data1(0)
        self._assign_solution_step_data2(1)

    def _assign_solution_step_data1(self, step=0):
        # generate nodal solution step test data
        random.seed(1.0)
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(POSITIVE_FACE_PRESSURE,step,random.random())
            node.SetSolutionStepValue(NEGATIVE_FACE_PRESSURE,step,random.random())

    def _assign_solution_step_data2(self, step=0):
        # generate nodal solution step test data
        random.seed(2.0)
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(POSITIVE_FACE_PRESSURE,step,random.random())
            node.SetSolutionStepValue(NEGATIVE_FACE_PRESSURE,step,random.random())

    def _zero_vector(self,size):
        v = Vector(size)
        for i in range(size):
            v[i] = 0.0
        return v

    def _transpose(self, m):
        tmp = Matrix(m.Size1(), m.Size2())
        for i in range(m.Size1()):
            for j in range(m.Size2()):
                tmp[i,j] = m[j,i]
        return tmp

    
    def _assert_vector_almost_equal(self, vector1, vector2, prec=7):
        self.assertEqual(vector1.Size(), vector2.Size())
        for i in range(vector1.Size()):
            self.assertAlmostEqual(vector1[i], vector2[i], prec)
    
    def _assert_matrix_almost_equal(self, matrix1, matrix2, prec=7):
        self.assertEqual(matrix1.Size1(), matrix2.Size1())
        self.assertEqual(matrix1.Size2(), matrix2.Size2())
        for i in range(matrix1.Size1()):
            for j in range(matrix1.Size2()):
                self.assertAlmostEqual(matrix1[i,j], matrix2[i,j], prec)

    def test_CalculateRHS(self):
        LHS = Matrix(3,3)
        RHS = self._zero_vector(3)
        self.potential_element.CalculateLocalSystem(LHS,RHS,self.model_part.ProcessInfo)
        
        potential = self._zero_vector(3)
        counter = 0
        for node in self.model_part.Nodes:
            potential[counter] = node.GetSolutionStepValue(POSITIVE_FACE_PRESSURE)
            counter +=1
        
        rhs2 = self._zero_vector(3)
        rhs2 -= LHS*potential
        self._assert_vector_almost_equal(RHS,rhs2)
        
    def test_CalculateLocalSystem(self):
        LHS = Matrix(3,3)
        RHS = self._zero_vector(3)
        self.potential_element.CalculateLocalSystem(LHS,RHS,self.model_part.ProcessInfo)
        
        dN_dxi_deta = Matrix(2,3)
        self.calculate_shape_functions_derivatives(dN_dxi_deta)
        
        dN_dxi_deta_trans = Matrix(3,2)
        self.calculate_shape_functions_derivatives_trans(dN_dxi_deta_trans)
        
        coordinates = Matrix(3,2)
        self.get_coordinates(coordinates)
        
        J = Matrix(2,2)
        J = dN_dxi_deta*coordinates
        
        J_inv = Matrix(2,2)
        self.calculate_inv_J(J,J_inv)
        
        J_inv_trans = Matrix(2,2)
        self.calculate_inv_J_trans(J,J_inv_trans)
        
        dN_dx_dy = Matrix(2,3)
        dN_dx_dy = J_inv*dN_dxi_deta
        
        dN_dx_dy_trans = Matrix(3,2)
        dN_dx_dy_trans = dN_dxi_deta_trans*J_inv_trans
        
        lhs2 = Matrix(3,3)
        lhs2 = dN_dx_dy_trans*dN_dx_dy/2
        self._assert_matrix_almost_equal(LHS,lhs2)
        
        
    def calculate_shape_functions_derivatives(self, dN_dxi_deta):
        dN_dxi_deta[0,0] = -1
        dN_dxi_deta[0,1] =  1
        dN_dxi_deta[0,2] =  0
        dN_dxi_deta[1,0] = -1
        dN_dxi_deta[1,1] =  0
        dN_dxi_deta[1,2] =  1
    
    def calculate_shape_functions_derivatives_trans(self, dN_dxi_deta_trans):
        dN_dxi_deta_trans[0,0] = -1
        dN_dxi_deta_trans[1,0] =  1
        dN_dxi_deta_trans[2,0] =  0
        dN_dxi_deta_trans[0,1] = -1
        dN_dxi_deta_trans[1,1] =  0
        dN_dxi_deta_trans[2,1] =  1
    
    
    def get_coordinates(self,coordinates):
        counter = 0
        for node in self.model_part.Nodes: #for counter, node in enumerate(self.model_part.Nodes)
            coordinates[counter,0] = node.X 
            coordinates[counter,1] = node.Y
            counter +=1
        
    def calculate_inv_J(self,J,J_inv):
        det_J = J[0,0]*J[1,1] - J[0,1]*J[1,0]
        J_inv[0,0] =  J[1,1]/det_J
        J_inv[0,1] = -J[0,1]/det_J
        J_inv[1,0] = -J[1,0]/det_J
        J_inv[1,1] =  J[0,0]/det_J
    
    def calculate_inv_J_trans(self,J,J_inv_trans):
        det_J = J[0,0]*J[1,1] - J[0,1]*J[1,0]
        J_inv_trans[0,0] =  J[1,1]/det_J
        J_inv_trans[0,1] = -J[1,0]/det_J
        J_inv_trans[1,0] = -J[0,1]/det_J
        J_inv_trans[1,1] =  J[0,0]/det_J

if __name__ == '__main__':
    KratosUnittest.main()
