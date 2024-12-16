#import python class test
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utilities

#import python packages
import numpy as np
from KratosMultiphysics.RomApplication.empirical_cubature_method import EmpiricalCubatureMethod
from KratosMultiphysics.RomApplication.randomized_singular_value_decomposition import RandomizedSingularValueDecomposition

def synthetic_matrix(degree, rows = 100):
    TestMatrix = np.zeros((rows,degree+1))
    x = np.linspace(0,1,rows)
    for i in range(degree+1):
        TestMatrix[:,i] = np.power(x,i)
    return TestMatrix

def calculate_basis(TestMatrix):
    # Calculate the randomized and truncated SVD of the matrix
    u,_,_,_ = RandomizedSingularValueDecomposition(COMPUTE_V=False).Calculate(TestMatrix)
    return u

class TestEmpiricalCubatureMethod(KratosUnittest.TestCase):

    def test_empirical_cubature_method(self):

        for degree in range(5,10):
            TestMatrix = synthetic_matrix(degree) #Get a synthetic matrix (During the training of a ROM model, this is a matrix of residuals projected onto a basis)
            TestMatrixBasis = calculate_basis(TestMatrix)
            ElementSelector = EmpiricalCubatureMethod()
            ElementSelector.SetUp(TestMatrixBasis, InitialCandidatesSet = np.array([0]), constrain_sum_of_weights=False)
            ElementSelector.Run()
            if ElementSelector.success is not True:
                ElementSelector.SetUp(TestMatrixBasis, InitialCandidatesSet = None, constrain_sum_of_weights=False)
                ElementSelector.Run()
            self.assertEqual( np.size(ElementSelector.z) , degree + 1) #for a polynomial of degree n, n+1 points are to be selected
    # Cleaning
    kratos_utilities.DeleteDirectoryIfExisting("__pycache__")

if __name__=='__main__':
    KratosUnittest.main()