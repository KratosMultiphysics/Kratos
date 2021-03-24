#import python class test
import KratosMultiphysics.KratosUnittest as KratosUnittest

#import python packages
try:
    import numpy as np
    from KratosMultiphysics.RomApplication.empirical_cubature_method import EmpiricalCubatureMethod
    numpy_available = True
except:
    numpy_available = False

def synthetic_matrix(degree, rows = 100):
    TestMatrix = np.zeros((rows,degree+1))
    x = np.linspace(0,1,rows)
    for i in range(degree+1):
        TestMatrix[:,i] = np.power(x,i)
    return TestMatrix

class TestEmpiricalCubatureMethod(KratosUnittest.TestCase):

    @KratosUnittest.skipIf(numpy_available == False, "numpy is required for RomApplication")
    def test_empirical_cubature_method(self):

        for degree in range(5,10):
            TestMatrix = synthetic_matrix(degree) #Get a synthetic matrix (During the training of a ROM model, this is a matrix of residuals projected onto a basis)

            #Pass the matrix to the ECM and obtain a set of elements(points) and weights (these steps are contained in the Run method of the ElementSelector base class)
            ElementSelector = EmpiricalCubatureMethod(SVD_tolerance = 0,  ECM_tolerance = 0)
            ElementSelector.SetUp(TestMatrix, 'test_number_of_elements', 'test_model_part_name')
            ElementSelector.Initialize()
            ElementSelector.Calculate()

            self.assertEqual( len(ElementSelector.z) , degree + 1 ) #for a polynomial of degree n, n+1 points are to be selected

if __name__=='__main__':
    KratosUnittest.main()