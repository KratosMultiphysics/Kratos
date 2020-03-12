from __future__ import print_function, absolute_import, division

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest

try:
    import numpy as np
    missing_numpy = False
except ImportError as e:
    missing_numpy = True

class TestNumpyExportDenseMatrix(KratosUnittest.TestCase):

    @KratosUnittest.skipIf(missing_numpy,"Missing python libraries (numpy)")
    def test_numpy_export_dense_matrix_no_copying(self):
        # Create a Kratos matrix
        KratosMatrix = KratosMultiphysics.Matrix(3,3)
        KratosMatrix.fill(1.0)

        # Export it to numpy array (No copying is performed)
        NumpyMatrix = np.array(KratosMatrix,copy=False)

        # Test correct creation
        for i in range(3):
            for j in range(3):
                self.assertEqual(KratosMatrix[i,j], NumpyMatrix[i,j])

        # Change an entry in Kratos Matrix
        KratosMatrix[0,0] = 123

        # Test change in Numpy matrix
        self.assertEqual(KratosMatrix[0,0], NumpyMatrix[0,0])

        # Change an entry in Numpy matrix
        NumpyMatrix[2,2] = 555

        # Test change in Kratos matrix
        self.assertEqual(KratosMatrix[2,2], NumpyMatrix[2,2])
        
    @KratosUnittest.skipIf(missing_numpy,"Missing python libraries (numpy)")
    def test_numpy_export_dense_matrix_copying(self):
        # Create a Kratos matrix
        KratosMatrix = KratosMultiphysics.Matrix(3,3)
        KratosMatrix.fill(1.0)

        # Export to numpy (Create a copy)
        CopyNumpyMatrix = np.array(KratosMatrix)

        # Test correct creation
        for i in range(3):
            for j in range(3):
                self.assertAlmostEqual(KratosMatrix[i,j], CopyNumpyMatrix[i,j], delta=1e-12)

        # Modify values in Kratos Matrix
        KratosMatrix.fill(555)

        # Test that values in CopyNumpyMatrix remain unchanged
        for i in range(3):
            for j in range(3):
                self.assertAlmostEqual(1.0, CopyNumpyMatrix[i,j], delta=1e-12)


if __name__ == '__main__':
    KratosUnittest.main()
