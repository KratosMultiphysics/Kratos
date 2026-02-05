
import os
import KratosMultiphysics

import KratosMultiphysics.LinearSolversApplication as LinearSolversApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.python_linear_solver_factory import ConstructSolver

class TestEigenDirectSolver(KratosUnittest.TestCase):
    def __ExecuteEigenDirectSolverTest(self,
                                       class_name: str,
                                       solver_type: str) -> None:
        space = KratosMultiphysics.UblasSparseSpace()

        settings = KratosMultiphysics.Parameters('{ "solver_type" : "LinearSolversApplication.' + solver_type + '" }')

        solver = ConstructSolver(settings)

        a = KratosMultiphysics.CompressedMatrix()

        this_file_dir = os.path.dirname(os.path.realpath(__file__))
        base_dir = os.path.dirname(os.path.dirname(os.path.dirname(this_file_dir)))
        matrix_file_path = os.path.join(base_dir, "kratos", "tests", "auxiliar_files_for_python_unittest", "sparse_matrix_files", "A.mm")

        file_read = KratosMultiphysics.ReadMatrixMarketMatrix(matrix_file_path, a) # symmetric test matrix
        self.assertTrue(file_read, msg="The MatrixFile could not be read")

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

    def __ExecuteEigenDirectComplexSolverTest(self,
                                              class_name: str,
                                              solver_type: str) -> None:
        space = KratosMultiphysics.UblasComplexSparseSpace()

        settings = KratosMultiphysics.Parameters('{ "solver_type" : "LinearSolversApplication.' + solver_type + '" }')

        solver = ConstructSolver(settings)

        a = KratosMultiphysics.CompressedMatrix()

        this_file_dir = os.path.dirname(os.path.realpath(__file__))
        base_dir = os.path.dirname(os.path.dirname(os.path.dirname(this_file_dir)))
        matrix_file_path = os.path.join(base_dir, "kratos", "tests", "auxiliar_files_for_python_unittest", "sparse_matrix_files", "A.mm")

        file_read = KratosMultiphysics.ReadMatrixMarketMatrix(matrix_file_path, a) # symmetric test matrix
        self.assertTrue(file_read, msg="The MatrixFile could not be read")

        a = KratosMultiphysics.ComplexCompressedMatrix(a)
        dimension = a.Size1()

        self.assertEqual(dimension, 900)

        b_exp = KratosMultiphysics.ComplexVector(dimension)

        for i in range(dimension):
            b_exp[i] = complex(i+1,i-1)

        x = KratosMultiphysics.ComplexVector(dimension)

        solver.Solve(a, x, b_exp)

        b_act = KratosMultiphysics.ComplexVector(dimension)
        space.Mult(a, x, b_act)

        for i in range(dimension):
            self.assertAlmostEqual(b_act[i], b_exp[i], 7)

    def test_EigenSparseLU(self):
        self.__ExecuteEigenDirectSolverTest('SparseLUSolver', 'sparse_lu')

    def test_EigenSparseCG(self):
        self.__ExecuteEigenDirectSolverTest('SparseCGSolver', 'sparse_cg')

    def test_EigenSparseQR(self):
        self.__ExecuteEigenDirectSolverTest('SparseQRSolver', 'sparse_qr')

    @KratosUnittest.skipIf(not LinearSolversApplication.HasMKL(), "Kratos was compiled without MKL support.")
    def test_EigenPardisoLU(self):
        self.__ExecuteEigenDirectSolverTest('PardisoLUSolver', 'pardiso_lu')

    @KratosUnittest.skipIf(not LinearSolversApplication.HasMKL(), "Kratos was compiled without MKL support.")
    def test_EigenPardisoLDLT(self):
        self.__ExecuteEigenDirectSolverTest('PardisoLDLTSolver', 'pardiso_ldlt')

    @KratosUnittest.skipIf(not LinearSolversApplication.HasMKL(), "Kratos was compiled without MKL support.")
    def test_EigenPardisoLLT(self):
        self.__ExecuteEigenDirectSolverTest('PardisoLLTSolver', 'pardiso_llt')

    def test_EigenComplexSparseLU(self):
        self.__ExecuteEigenDirectComplexSolverTest('ComplexSparseLUSolver', 'sparse_lu_complex')

    @KratosUnittest.skipIf(not LinearSolversApplication.HasMKL(), "Kratos was compiled without MKL support.")
    def test_EigenComplexPardisoLU(self):
        self.__ExecuteEigenDirectComplexSolverTest('ComplexPardisoLUSolver', 'pardiso_lu_complex')

    @KratosUnittest.skipIf(not LinearSolversApplication.HasMKL(), "Kratos was compiled without MKL support.")
    def test_EigenComplexPardisoLDLT(self):
        self.__ExecuteEigenDirectComplexSolverTest('ComplexPardisoLDLTSolver', 'pardiso_ldlt_complex')

    @KratosUnittest.skipIf(not LinearSolversApplication.HasMKL(), "Kratos was compiled without MKL support.")
    def test_EigenComplexPardisoLLT(self):
        self.__ExecuteEigenDirectComplexSolverTest('ComplexPardisoLLTSolver', 'pardiso_llt_complex')

    @KratosUnittest.skipIf(not LinearSolversApplication.HasSuiteSparse(), "Kratos was compiled without SuiteSparse support.")
    def test_EigenCholmod(self) -> None:
        self.__ExecuteEigenDirectSolverTest("CholmodSolver", "cholmod")

    @KratosUnittest.skipIf(not LinearSolversApplication.HasSuiteSparse(), "Kratos was compiled without SuiteSparse support.")
    def test_EigenUmfPack(self) -> None:
        self.__ExecuteEigenDirectSolverTest("UmfPackSolver", "umfpack")

    @KratosUnittest.skipIf(not LinearSolversApplication.HasSuiteSparse(), "Kratos was compiled without SuiteSparse support.")
    def test_EigenSPQR(self) -> None:
        self.__ExecuteEigenDirectSolverTest("SPQRSolver", "spqr")

    @KratosUnittest.skipIf(not LinearSolversApplication.HasSuiteSparse(), "Kratos was compiled without SuiteSparse support.")
    def test_ComplexEigenUmfPack(self) -> None:
        self.__ExecuteEigenDirectComplexSolverTest("ComplexUmfPackSolver", "umfpack_complex")

    @KratosUnittest.skipIf(not LinearSolversApplication.HasSuiteSparse(), "Kratos was compiled without SuiteSparse support.")
    def test_ComplexEigenSPQR(self) -> None:
        self.__ExecuteEigenDirectComplexSolverTest("ComplexSPQRSolver", "spqr_complex")

if __name__ == '__main__':
    KratosUnittest.main()