import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import numpy as np 

class TestNumpyExportDenseMatrix(KratosUnittest.TestCase):

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

    def test_numpy_export_complex_dense_matrix_no_copying(self):
        # Create a Kratos matrix
        KratosMatrix = KratosMultiphysics.ComplexMatrix(3,3)
        KratosMatrix.fill(1.0+1.0j)

        # Export it to numpy array (No copying is performed)
        NumpyMatrix = np.array(KratosMatrix,copy=False)

        # Test correct creation
        for i in range(3):
            for j in range(3):
                self.assertEqual(KratosMatrix[i,j], NumpyMatrix[i,j])

        # Change an entry in Kratos Matrix
        KratosMatrix[0,0] = 123-123j

        # Test change in Numpy matrix
        self.assertEqual(KratosMatrix[0,0], NumpyMatrix[0,0])

        # Change an entry in Numpy matrix
        NumpyMatrix[2,2] = 555-555j

        # Test change in Kratos matrix
        self.assertEqual(KratosMatrix[2,2], NumpyMatrix[2,2])

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

    
    def test_numpy_import_dense_matrix(self):
        # Create numpy array
        np_array = np.ones((3,2))
        np_array[0,0] = 12
        np_array[2,1] = 42

        # Import (i.e. copy) to Kratos
        kratos_matrix = KratosMultiphysics.Matrix(np_array)

        # Test
        for i in range(3):
            for j in range(2):
                self.assertAlmostEqual(kratos_matrix[i,j], np_array[i,j], delta=1e-12)

    def test_numpy_import_complex_dense_matrix(self):
        # Create numpy array
        np_array = np.ones((3,2), dtype=complex)
        np_array[0,0] = 12.-3.j
        np_array[2,1] = 42.-7.j

        # Assert no double array can be created with complex data
        with self.assertRaises(RuntimeError):
            KratosMultiphysics.Matrix(np_array)

        # Import (i.e. copy) to Kratos
        kratos_matrix = KratosMultiphysics.ComplexMatrix(np_array)

        # Test
        for i in range(3):
            for j in range(2):
                self.assertAlmostEqual(kratos_matrix[i,j], np_array[i,j], delta=1e-12)

    def test_numpy_import_real_to_complex_dense_matrix(self):
        # Create numpy array
        np_array = np.ones((3,2))
        np_array[0,0] = 12.
        np_array[2,1] = 42.

        # Import (i.e. copy) to Kratos
        kratos_matrix = KratosMultiphysics.ComplexMatrix(np_array)

        # Test
        for i in range(3):
            for j in range(2):
                self.assertAlmostEqual(kratos_matrix[i,j].real, np_array[i,j], delta=1e-12)
                self.assertAlmostEqual(kratos_matrix[i,j].imag, 0.0, delta=1e-12)

    def test_numpy_import_dense_matrix_format(self):
        # Create three-dimensional numpy array
        np_array = np.ones((3,2,4))
        with self.assertRaises(RuntimeError):
            KratosMultiphysics.Matrix(np_array)

if __name__ == '__main__':
    KratosUnittest.main()
