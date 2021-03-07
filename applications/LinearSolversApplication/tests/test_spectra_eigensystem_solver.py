
from __future__ import print_function, absolute_import, division
import os
import KratosMultiphysics
import KratosMultiphysics.LinearSolversApplication as LinearSolversApplication

import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics import eigen_solver_factory

class TestSpectraEigensystemSolver(KratosUnittest.TestCase):
    def test_mass_normalization(self):
        # check if solver is available
        if (not hasattr(LinearSolversApplication, "SpectraSymGEigsShiftSolver")):
            self.skipTest("SpectraSymGEigsShiftSolver is not included in the compilation of the LinearSolversApplication")

        space = KratosMultiphysics.UblasSparseSpace()

        settings = KratosMultiphysics.Parameters('''{
            "solver_type": "spectra_sym_g_eigs_shift",
            "number_of_eigenvalues": 3,
            "normalize_eigenvectors": true,
            "shift": 0.0,
            "echo_level": 1
        }''')

        K = KratosMultiphysics.CompressedMatrix()

        this_file_dir = os.path.dirname(os.path.realpath(__file__))
        base_dir = os.path.dirname(os.path.dirname(os.path.dirname(this_file_dir)))
        matrix_file_path = os.path.join(base_dir, "kratos", "tests", "auxiliar_files_for_python_unittest", "sparse_matrix_files", "A.mm")

        file_read = KratosMultiphysics.ReadMatrixMarketMatrix(matrix_file_path, K) # symmetric test matrix
        self.assertTrue(file_read, msg="The MatrixFile could not be read")

        n = K.Size1()
        self.assertEqual(n, 900)

        M = KratosMultiphysics.CompressedMatrix(n, n)

        for i in range(n):
            for j in range(n):
                if (i == j):
                    M[i, j] = 1.0

        # create result containers (they will be resized inside the solver)
        eigenvalues = KratosMultiphysics.Vector(n)
        eigenvectors = KratosMultiphysics.Matrix(n, 1)

        # Construct the solver
        eigen_solver = eigen_solver_factory.ConstructSolver(settings)

        # Solve
        eigen_solver.Solve(K, M, eigenvalues, eigenvectors)

        self.assertAlmostEqual(eigenvalues[0], 0.061463, 4)
        self.assertAlmostEqual(eigenvalues[1], 0.15318431112733275, 7)
        self.assertAlmostEqual(eigenvalues[2], 0.153184311127333, 7)


        # test mass normalization of eigenvectors and if they are correct
        for i in range(eigenvectors.Size1()):
            eigenvector = KratosMultiphysics.Vector(n)
            for j in range(n):
                eigenvector[j] = eigenvectors[i,j]

            _aux = KratosMultiphysics.Vector(n)
            space.Mult(M, eigenvector, _aux)

            value = 0.0
            for j in range(n):
                value += eigenvector[j] * _aux[j]

            self.assertAlmostEqual(value, 1.0, 7)


if __name__ == '__main__':
    KratosUnittest.main()
