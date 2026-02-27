
# --- Kratos Imports ---
import KratosMultiphysics
import KratosMultiphysics.LinearSolversApplication as LinearSolversApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.python_linear_solver_factory import ConstructSolver

# --- STD Imports ---
import pathlib


class TestMKLSmoothers(KratosUnittest.TestCase):

    def __Run(self, solver_type: str):
        settings: KratosMultiphysics.Parameters = KratosMultiphysics.Parameters('{"solver_type" : "LinearSolversApplication.' + solver_type + '" }')
        solver: KratosMultiphysics.LinearSolver = ConstructSolver(settings)

        this_file_dir: pathlib.Path = pathlib.Path(__file__).absolute().parent
        base_dir: pathlib.Path = this_file_dir.parent.parent.parent
        matrix_file_path: pathlib.Path = base_dir / "kratos" / "tests" / "auxiliar_files_for_python_unittest" / "sparse_matrix_files" / "A.mm"

        lhs = KratosMultiphysics.CompressedMatrix()
        file_read = KratosMultiphysics.ReadMatrixMarketMatrix(str(matrix_file_path), lhs)
        self.assertTrue(file_read, msg="The MatrixFile could not be read")

        dimension = lhs.Size1()
        rhs = KratosMultiphysics.Vector(dimension)

        for i in range(dimension):
            rhs[i] = i + 1

        # Iteratively call the smoother on the residual, and check whether it decreases its norm.
        residual_norm: float = KratosMultiphysics.UblasSparseSpace().TwoNorm(rhs)
        for _ in range(10):
            solution = KratosMultiphysics.Vector(dimension)
            solver.Solve(lhs, solution, rhs)
            residual = KratosMultiphysics.Vector(dimension)
            KratosMultiphysics.UblasSparseSpace().Mult(lhs, solution, residual)
            rhs -= residual
            new_residual_norm: float = KratosMultiphysics.UblasSparseSpace().TwoNorm(rhs)
            self.assertLess(new_residual_norm, residual_norm)
            residual_norm = new_residual_norm


    @KratosUnittest.skipIf(not LinearSolversApplication.HasMKL(), "Kratos was compiled without MKL support.")
    def test_MKLILU0(self):
        self.__Run('mkl_ilu0')


    @KratosUnittest.skipIf(not LinearSolversApplication.HasMKL(), "Kratos was compiled without MKL support.")
    def test_MKLILUT(self):
        self.__Run('mkl_ilut')


if __name__ == '__main__':
    KratosUnittest.main()
