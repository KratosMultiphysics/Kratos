#import python class test
import KratosMultiphysics.KratosUnittest as KratosUnittest

#import python packages
from KratosMultiphysics.RomApplication.randomized_singular_value_decomposition import RandomizedSingularValueDecomposition
import numpy as np
import numpy.matlib

def synthetic_matrix(degree, rows = 100,repetitions=20):
    TestMatrix = np.zeros((rows,degree))
    x = np.linspace(0,1,rows)
    for i in range(degree):
        TestMatrix[:,i] = np.power(x,i)          
    return np.matlib.repmat(TestMatrix, 1, repetitions)

class TestRandomizedSVD(KratosUnittest.TestCase):

    def test_radomized_svd(self):

        svd_truncation_tolerance = 1e-5
        for rank in range(5,10):        
            TestMatrix = synthetic_matrix(rank) #create a matrix of known rank using polynomials
            U,S,V,error = RandomizedSingularValueDecomposition().Calculate(TestMatrix, svd_truncation_tolerance) #calculate randomized svd
            Randomized_Reconstruction = U@np.diag(S)@V.T #reconstruct matrix

            #check that the difference of the reconstruction is below tolerance
            self.assertLess( np.linalg.norm(Randomized_Reconstruction - TestMatrix), svd_truncation_tolerance*np.linalg.norm(TestMatrix))

if __name__=='__main__':
    KratosUnittest.main()