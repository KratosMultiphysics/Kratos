
from __future__ import print_function, absolute_import, division
import os
import KratosMultiphysics

import KratosMultiphysics.EigenSolversApplication as EigenSolversApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.python_linear_solver_factory import ConstructSolver

class TestEigenDenseDirectSolver(KratosUnittest.TestCase):
    def _execute_eigen_dense_direct_solver_test(self, class_name, solver_type):
        # check if solver is available
        if (not hasattr(EigenSolversApplication, class_name)):
            self.skipTest(class_name + " is not included in the compilation of the EigenSolversApplication")

        settings = KratosMultiphysics.Parameters('{ "solver_type" : "EigenSolversApplication.' + solver_type + '" }')

        solver = ConstructSolver(settings)

        A = KratosMultiphysics.Matrix(3,3)
        val = 0
        for i in range(3):
            for j in range(3):
                val += 1
                A[i,j] = val
        A[2,2] = 10

        b_act = KratosMultiphysics.Vector([3., 3., 4.])
        x = KratosMultiphysics.Vector(3)

        solver.Solve(A, x, b_act)

        b_exp = A * x

        for i in range(3):
            self.assertAlmostEqual(b_act[i], b_exp[i], 7)

    def _execute_eigen_dense_direct_complex_solver_test(self, class_name, solver_type):
        # check if solver is available
        if (not hasattr(EigenSolversApplication, class_name)):
            self.skipTest(class_name + " is not included in the compilation of the EigenSolversApplication")
            
        settings = KratosMultiphysics.Parameters('{ "solver_type" : "EigenSolversApplication.' + solver_type + '" }')

        solver = ConstructSolver(settings)

        A = KratosMultiphysics.ComplexMatrix(3,3)
        val = 0
        for i in range(3):
            for j in range(3):
                val += 1+2j
                A[i,j] = val
        A[2,2] = 10

        b_act = KratosMultiphysics.ComplexVector([3., 3., 4.-2j])
        x = KratosMultiphysics.ComplexVector(3)

        solver.Solve(A, x, b_act)

        b_exp = A * x

        for i in range(3):
            self.assertAlmostEqual(b_act[i], b_exp[i], 7)

    def test_eigen_dense_colpivhouseholderqr(self):
        self._execute_eigen_dense_direct_solver_test('DenseColPivHouseholderQRSolver', 'dense_colpivhouseholderqr')

    def test_eigen_dense_householderqr(self):
        self._execute_eigen_dense_direct_solver_test('DenseHouseholderQRSolver', 'dense_householderqr')

    def test_eigen_dense_llt(self):
        self._execute_eigen_dense_direct_solver_test('DenseLLTSolver', 'dense_llt')

    def test_eigen_dense_colpivhouseholderqr_complex(self):
        self._execute_eigen_dense_direct_complex_solver_test('ComplexDenseColPivHouseholderQRSolver', 'complex_dense_colpivhouseholderqr')

    def test_eigen_dense_householderqr_complex(self):
        self._execute_eigen_dense_direct_complex_solver_test('ComplexDenseHouseholderQRSolver', 'complex_dense_householderqr')

    def test_eigen_dense_llt_complex(self):
        self._execute_eigen_dense_direct_complex_solver_test('ComplexDenseLLTSolver', 'complex_dense_llt')

if __name__ == '__main__':
    KratosUnittest.main()