import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils

import numpy as np 

try:
    import scipy
    import scipy.sparse
    import KratosMultiphysics.scipy_conversion_tools
    scipy_available = True
except ImportError:
    scipy_available = False

class TestSparseMatrixInterface(KratosUnittest.TestCase):

    def test_scalar_assembly(self):
        G = KratosMultiphysics.SparseContiguousRowGraph(3)
        G.AddEntry(0,0)
        G.AddEntry(2,2)
        G.AddEntry(1,2)
        G.Finalize()

        A = KratosMultiphysics.CsrMatrix(G)
        #5 0 0
        #0 0 7
        #0 0 6
        A.BeginAssemble()
        A.AssembleEntry(5.0,0,0)
        A.AssembleEntry(6.0,2,2)
        A.AssembleEntry(7.0,1,2)
        A.FinalizeAssemble()

        x = KratosMultiphysics.SystemVector(3)
        x.SetValue(1.0)
        y = KratosMultiphysics.SystemVector(3)
        y.SetValue(0.0)
        A.SpMV(x,y)
        self.assertEqual(y[0],5.0)
        self.assertEqual(y[1],7.0)
        self.assertEqual(y[2],6.0)

        y = KratosMultiphysics.SystemVector(3)
        y.SetValue(0.0)
        A.TransposeSpMV(x,y)
        self.assertEqual(y[0],5.0)
        self.assertEqual(y[1],0.0)
        self.assertEqual(y[2],13.0)

    def test_matrix_assembly(self):
        #data to be assembled
        values = np.array([[1.0,-1.0],[-1.0,1.0]])
        
        x = KratosMultiphysics.SystemVector(3)
        x[0] = 1.0
        x[1] = 2.0
        x[2] = 3.0
        y = KratosMultiphysics.SystemVector(3)
        y.SetValue(0.0)

        #construction of matrix graph (by vector input)
        G = KratosMultiphysics.SparseContiguousRowGraph(3)
        G.AddEntries(np.array([0,1])) #for square matrices
        G.AddEntries([1,2],[1,2]) #thos would allow non square input
        G.Finalize()

        #assembling matrix
        A = KratosMultiphysics.CsrMatrix(G)
        A.BeginAssemble()
        A.Assemble(values,[0,1]) #input by list
        A.Assemble(values,np.array([1,2])) #input by numpy array
        A.FinalizeAssemble()

        A.SpMV(x,y)

        validation_data = [ 2., -3.,  1., -3.,  6., -3.,  1., -3.,  2.]                                                                                                                                 
        validation_index2 = [0, 1, 2, 0, 1, 2, 0, 1, 2]                                                                                                                                                 
        validation_index1 = [0, 3, 6, 9]    

        B = A.SpMM(A)

        for i in range(len(validation_data)):
            self.assertEqual(B.value_data()[i], validation_data[i])
            self.assertEqual(B.index2_data()[i], validation_index2[i])

        for i in range(len(validation_index1)):
            self.assertEqual(B.index1_data()[i], validation_index1[i])

        # the following should be added back in case scipy support is enabled in testing
        if scipy_available:
             #test conversion to scipy matrix
            Ascipy = KratosMultiphysics.scipy_conversion_tools.to_csr(A)

            #verify compatibility with numpy arrays
            xscipy = np.array(x).T
            yscipy = Ascipy@x 

            self.assertEqual(yscipy[0],-1.0)
            self.assertEqual(yscipy[1],0.0)
            self.assertEqual(yscipy[2],1.0)

            B_scipy = Ascipy@Ascipy
            B_scipy.sort_indices() #it is crucial that the index are sorted to do the following comparison

            for i in range(len(validation_data)):
                self.assertEqual(B_scipy.data[i], validation_data[i])
                self.assertEqual(B_scipy.indices[i], validation_index2[i])

            for i in range(len(validation_index1)):
                self.assertEqual(B_scipy.indptr[i], validation_index1[i])


if __name__ == '__main__':
    KratosUnittest.main()