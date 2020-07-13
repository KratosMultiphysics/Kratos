
from __future__ import print_function, absolute_import, division

import KratosMultiphysics

import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics import eigen_solver_factory

class TestDenseEigenvalueSolver(KratosUnittest.TestCase):
    def test_eigen_dense_eigenvalue_solver(self):
        eigensolver_settings = KratosMultiphysics.Parameters('''{
            "solver_type"       : "dense_eigensolver",
            "ascending_order"   : true,
            "echo_level"        : 0
         }''')

        # Construct matrix
        n = 3
        A = KratosMultiphysics.Matrix(n,n)
        A[0,0] = 3
        A[0,1] = 2
        A[0,2] = -1
        A[1,0] = 2
        A[1,1] = 10
        A[1,2] = 2
        A[2,0] = -1
        A[2,1] = 2
        A[2,2] = 3

        # Create dummy matrix
        Dummy = KratosMultiphysics.Matrix(n, n)

        # Create result containers (they will be resized inside the solver)
        eigenvalues = KratosMultiphysics.Vector(n)
        eigenvectors = KratosMultiphysics.Matrix(n,n)

        # Construct the solver
        eigen_solver = eigen_solver_factory.ConstructSolver(eigensolver_settings)
        # Solve
        eigen_solver.Solve(A, Dummy, eigenvalues, eigenvectors)

        # Check eigenvalues
        self.assertAlmostEqual(eigenvalues[0], 1.101020514433644, 8)
        self.assertAlmostEqual(eigenvalues[1], 4.0, 8)
        self.assertAlmostEqual(eigenvalues[2], 10.898979485566358, 8)
        # Check eigenvectors
        self.assertAlmostEqual(abs(eigenvectors[0,0]), 0.673887338679049, 8)
        self.assertAlmostEqual(abs(eigenvectors[1,0]), 0.302905446527686, 8)
        self.assertAlmostEqual(abs(eigenvectors[2,0]), 0.673887338679049, 8)

        self.assertAlmostEqual(abs(eigenvectors[0,1]), 0.707106781186548, 8)
        self.assertAlmostEqual(abs(eigenvectors[1,1]), 0.0, 8)
        self.assertAlmostEqual(abs(eigenvectors[2,1]), 0.707106781186547, 8)

        self.assertAlmostEqual(abs(eigenvectors[0,2]), 0.214186495298066, 8)
        self.assertAlmostEqual(abs(eigenvectors[1,2]), 0.953020613871422, 8)
        self.assertAlmostEqual(abs(eigenvectors[2,2]), 0.214186495298066, 8)

if __name__ == '__main__':
    KratosUnittest.main()