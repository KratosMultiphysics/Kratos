
from __future__ import print_function, absolute_import, division
import KratosMultiphysics

import KratosMultiphysics.EigenSolversApplication as EigenSolversApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest
from new_linear_solver_factory import ConstructSolver

class TestSparseEigensystemSolver(KratosUnittest.TestCase):
    def test_sparse_eigensystem_solver(self):
        # check if solver is available
        if (not hasattr(EigenSolversApplication, 'SparseEigensystemSolver')):
            self.skipTest("SparseEigensystemSolver is not included in the compilation of the EigenSolversApplication")

        space = KratosMultiphysics.UblasSparseSpace()

        settings = KratosMultiphysics.Parameters('''{
            "solver_type": "eigen_sparse_eigensystem",
            "number_of_eigenvalues" : 3,
            "tolerance" : 1e-6,
            "max_iteration" : 10,
            "echo_level": 2
        }''')
        
        solver = EigenSolversApplication.SparseEigensystemSolver(settings)

        a = KratosMultiphysics.CompressedMatrix()
        b = KratosMultiphysics.CompressedMatrix()
        
        KratosMultiphysics.ReadMatrixMarketMatrix(r'C:\Users\oberbichler\source\matrices\hook_01\hook_01k.mtx', a)
        KratosMultiphysics.ReadMatrixMarketMatrix(r'C:\Users\oberbichler\source\matrices\hook_01\hook_01m.mtx', b)
        
        dimension = a.Size1()

        eigenvalues = KratosMultiphysics.Vector(dimension)
        eigenvectors = KratosMultiphysics.Matrix(dimension, 3)

        solver.Solve(a, b, eigenvalues, eigenvectors)

        # lam = eigenvalues[0]
        # x = eigenvectors[0]

        # k_x = KratosMultiphysics.Vector(dimension)
        # k_x = KratosMultiphysics.Vector(dimension)
        # space.Mult(k, x, k_x)
        # space.Mult(m, x, k_x)

        # for i in range(dimension):
        #     self.assertAlmostEqual(b_act[i], b_exp[i], 7)

if __name__ == '__main__':
    KratosUnittest.main()