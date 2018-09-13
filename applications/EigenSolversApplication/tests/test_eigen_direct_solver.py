
from __future__ import print_function, absolute_import, division
import os
import KratosMultiphysics

import KratosMultiphysics.EigenSolversApplication as EigenSolversApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest
from new_linear_solver_factory import ConstructSolver

class TestEigenDirectSolver(KratosUnittest.TestCase):
    def _execute_eigen_direct_solver_test(self, class_name, solver_type):
        # check if solver is available
        if (not hasattr(EigenSolversApplication, class_name)):
            self.skipTest(class_name + " is not included in the compilation of the EigenSolversApplication")

        space = KratosMultiphysics.UblasSparseSpace()

        settings = KratosMultiphysics.Parameters('{ "solver_type" : "' + solver_type + '" }')

        solver = ConstructSolver(settings)

        a = KratosMultiphysics.CompressedMatrix()

        this_file_dir = os.path.dirname(os.path.realpath(__file__))
        base_dir = os.path.dirname(os.path.dirname(os.path.dirname(this_file_dir)))
        matrix_file_path = os.path.join(base_dir, "kratos", "tests", "A.mm")

        KratosMultiphysics.ReadMatrixMarketMatrix(matrix_file_path, a) # symmetric test matrix

        dimension = a.Size1()

        self.assertEqual(dimension, 900)

        b_exp = KratosMultiphysics.Vector(dimension) # [1, 2, ..., dimension-1, dimension]

        for i in range(dimension):
            b_exp[i] = i + 1

        x = KratosMultiphysics.Vector(dimension)

        solver.Solve(a, x, b_exp)

        b_act = KratosMultiphysics.Vector(dimension)
        space.Mult(a, x, b_act)

        for i in range(dimension):
            self.assertAlmostEqual(b_act[i], b_exp[i], 7)

    def test_eigen_sparse_lu(self):
        self._execute_eigen_direct_solver_test('SparseLUSolver', 'eigen_sparse_lu')

    def test_eigen_pardiso_lu(self):
        self._execute_eigen_direct_solver_test('PardisoLUSolver', 'eigen_pardiso_lu')

    def test_eigen_pardiso_ldlt(self):
        self._execute_eigen_direct_solver_test('PardisoLDLTSolver', 'eigen_pardiso_ldlt')

    def test_eigen_pardiso_llt(self):
        self._execute_eigen_direct_solver_test('PardisoLLTSolver', 'eigen_pardiso_llt')

if __name__ == '__main__':
    KratosUnittest.main()