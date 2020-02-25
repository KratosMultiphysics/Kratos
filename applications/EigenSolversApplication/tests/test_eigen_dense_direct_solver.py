
from __future__ import print_function, absolute_import, division
import KratosMultiphysics

import KratosMultiphysics.EigenSolversApplication as EigenSolversApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.EigenSolversApplication import dense_linear_solver_factory

class TestEigenDenseDirectSolver(KratosUnittest.TestCase):
    def _real_eq_system(self):
        A = KratosMultiphysics.Matrix(3,3)
        val = 0
        for i in range(3):
            for j in range(3):
                val += 1
                A[i,j] = val
        A[2,2] = 10

        b = KratosMultiphysics.Vector([3., 3., 4.])
        x = KratosMultiphysics.Vector(3)

        return A, b, x

    def _real_posdef_eq_system(self):
        A = KratosMultiphysics.Matrix(3,3,0.)
        A[0,0] = 4
        A[0,1] = -1
        A[1,0] = -1
        A[2,0] = 2
        A[0,2] = 2
        A[1,1] = 6
        A[2,2] = 5

        b = KratosMultiphysics.Vector([3., 3., 4.])
        x = KratosMultiphysics.Vector(3)

        return A, b, x

    def _cplx_eq_system(self):
        A = KratosMultiphysics.ComplexMatrix(3,3)
        val = 0
        for i in range(3):
            for j in range(3):
                val += 1+2j
                A[i,j] = val
        A[2,2] = 10

        b = KratosMultiphysics.ComplexVector([3., 3., 4.-2.j])
        x = KratosMultiphysics.ComplexVector(3)
        
        return A, b, x

    def _execute_eigen_dense_direct_solver_test(self, class_name, solver_type, eq_system_type):
        # check if solver is available
        if (not hasattr(EigenSolversApplication, class_name)):
            self.skipTest(class_name + " is not included in the compilation of the EigenSolversApplication")

        settings = KratosMultiphysics.Parameters('{ "solver_type" : "EigenSolversApplication.' + solver_type + '" }')

        solver = dense_linear_solver_factory.ConstructSolver(settings)
        
        A, b_act, x = eq_system_type()

        solver.Solve(A, x, b_act)

        b_exp = A * x

        for i in range(3):
            self.assertAlmostEqual(b_act[i], b_exp[i], 7)

    def test_eigen_dense_colpivhouseholderqr(self):
        self._execute_eigen_dense_direct_solver_test('DenseColPivHouseholderQRSolver', 'dense_col_piv_householder_qr',self._real_eq_system)

    def test_eigen_dense_householderqr(self):
        self._execute_eigen_dense_direct_solver_test('DenseHouseholderQRSolver', 'dense_householder_qr',self._real_eq_system)

    def test_eigen_dense_llt(self):
        self._execute_eigen_dense_direct_solver_test('DenseLLTSolver', 'dense_llt', self._real_posdef_eq_system)

    def test_eigen_dense_partialpivlu(self):
        self._execute_eigen_dense_direct_solver_test('DensePartialPivLUSolver', 'dense_partial_piv_lu', self._real_eq_system)

    def test_eigen_dense_colpivhouseholderqr_complex(self):
        self._execute_eigen_dense_direct_solver_test('ComplexDenseColPivHouseholderQRSolver', 'complex_dense_col_piv_householder_qr', self._cplx_eq_system)

    def test_eigen_dense_householderqr_complex(self):
        self._execute_eigen_dense_direct_solver_test('ComplexDenseHouseholderQRSolver', 'complex_dense_householder_qr', self._cplx_eq_system)

    def test_eigen_dense_partialpivlu_complex(self):
        self._execute_eigen_dense_direct_solver_test('ComplexDensePartialPivLUSolver', 'complex_dense_partial_piv_lu', self._cplx_eq_system)

if __name__ == '__main__':
    KratosUnittest.main()