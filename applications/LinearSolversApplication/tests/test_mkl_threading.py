
# --- Kratos Imports ---
import KratosMultiphysics
import KratosMultiphysics.LinearSolversApplication as LinearSolversApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.python_linear_solver_factory import ConstructSolver

# --- STD Imports ---
from pathlib import Path

# Lazily evaluated: whether MKLUtilities.SetNumThreads is actually reflected by
# GetNumThreads. Environment variables such as MKL_NUM_THREADS or MKL_DYNAMIC
# can override the runtime setting, in which case thread-count assertions are
# meaningless and the affected tests skip.
_mkl_thread_control = None

def _HasMKLThreadControl() -> bool:
    global _mkl_thread_control
    if _mkl_thread_control is None:
        saved_threads = LinearSolversApplication.MKLUtilities.GetNumThreads()
        LinearSolversApplication.MKLUtilities.SetNumThreads(1)
        _mkl_thread_control = LinearSolversApplication.MKLUtilities.GetNumThreads() == 1
        LinearSolversApplication.MKLUtilities.SetNumThreads(saved_threads)
    return _mkl_thread_control


@KratosUnittest.skipIf(not LinearSolversApplication.HasMKL(), "Kratos was compiled without MKL support.")
class TestMKLUtilities(KratosUnittest.TestCase):
    """Unit tests for the MKLUtilities threading-control layer.

    Conventions under test (see custom_utilities/mkl_utilities.cpp):
    - CheckThreadNumber returns True when the current MKL thread count already
      satisfies the requested setting, False when an adjustment would be needed.
    - ComputeMKLThreadCount returns None when no adjustment is needed, or the
      thread count the caller should apply.
    """

    def setUp(self):
        self.initial_mkl_threads = LinearSolversApplication.MKLUtilities.GetNumThreads()
        self.initial_omp_threads = KratosMultiphysics.ParallelUtilities.GetNumThreads()

    def tearDown(self):
        KratosMultiphysics.ParallelUtilities.SetNumThreads(self.initial_omp_threads)
        LinearSolversApplication.MKLUtilities.SetNumThreads(self.initial_mkl_threads)

    def _SkipWithoutThreadControl(self):
        if not _HasMKLThreadControl():
            self.skipTest("MKL thread count is overridden by the environment (e.g. MKL_NUM_THREADS).")

    def test_GetSetNumThreads(self):
        self._SkipWithoutThreadControl()
        LinearSolversApplication.MKLUtilities.SetNumThreads(1)
        self.assertEqual(LinearSolversApplication.MKLUtilities.GetNumThreads(), 1)
        if KratosMultiphysics.ParallelUtilities.GetNumProcs() > 1:
            LinearSolversApplication.MKLUtilities.SetNumThreads(2)
            self.assertEqual(LinearSolversApplication.MKLUtilities.GetNumThreads(), 2)

    def test_CheckThreadNumber_DoNothing(self):
        # Do_nothing never requires an adjustment, regardless of the current state.
        self.assertTrue(LinearSolversApplication.MKLUtilities.CheckThreadNumber(0))
        self.assertTrue(LinearSolversApplication.MKLUtilities.CheckThreadNumber(
            LinearSolversApplication.MKLThreadSetting.Do_nothing))

    def test_CheckThreadNumber_Manual(self):
        self._SkipWithoutThreadControl()
        LinearSolversApplication.MKLUtilities.SetNumThreads(1)
        self.assertTrue(LinearSolversApplication.MKLUtilities.CheckThreadNumber(1))
        if KratosMultiphysics.ParallelUtilities.GetNumProcs() > 1:
            self.assertFalse(LinearSolversApplication.MKLUtilities.CheckThreadNumber(2))

    def test_CheckThreadNumber_Consistent(self):
        self._SkipWithoutThreadControl()
        # MKL = 1 <= OMP threads: always consistent.
        LinearSolversApplication.MKLUtilities.SetNumThreads(1)
        self.assertTrue(LinearSolversApplication.MKLUtilities.CheckThreadNumber(
            LinearSolversApplication.MKLThreadSetting.Consistent))
        if KratosMultiphysics.ParallelUtilities.GetNumProcs() > 1:
            # MKL = 2 > OMP = 1: inconsistent.
            KratosMultiphysics.ParallelUtilities.SetNumThreads(1)
            LinearSolversApplication.MKLUtilities.SetNumThreads(2)
            self.assertFalse(LinearSolversApplication.MKLUtilities.CheckThreadNumber(
                LinearSolversApplication.MKLThreadSetting.Consistent))

    def test_CheckThreadNumber_Minimal(self):
        self._SkipWithoutThreadControl()
        # MKL = 1 == min(1, OMP): already minimal.
        LinearSolversApplication.MKLUtilities.SetNumThreads(1)
        self.assertTrue(LinearSolversApplication.MKLUtilities.CheckThreadNumber(
            LinearSolversApplication.MKLThreadSetting.Minimal))
        if KratosMultiphysics.ParallelUtilities.GetNumProcs() > 1:
            # MKL = 2 != min(2, 1): needs adjustment.
            KratosMultiphysics.ParallelUtilities.SetNumThreads(1)
            LinearSolversApplication.MKLUtilities.SetNumThreads(2)
            self.assertFalse(LinearSolversApplication.MKLUtilities.CheckThreadNumber(
                LinearSolversApplication.MKLThreadSetting.Minimal))

    def test_CheckThreadNumber_Invalid(self):
        with self.assertRaises(RuntimeError):
            LinearSolversApplication.MKLUtilities.CheckThreadNumber(-3)

    def test_ComputeMKLThreadCount_DoNothing(self):
        self.assertIsNone(LinearSolversApplication.MKLUtilities.ComputeMKLThreadCount(0))
        self.assertIsNone(LinearSolversApplication.MKLUtilities.ComputeMKLThreadCount(
            LinearSolversApplication.MKLThreadSetting.Do_nothing))

    def test_ComputeMKLThreadCount_Manual(self):
        self._SkipWithoutThreadControl()
        LinearSolversApplication.MKLUtilities.SetNumThreads(1)
        # Already at the requested count: no adjustment.
        self.assertIsNone(LinearSolversApplication.MKLUtilities.ComputeMKLThreadCount(1))
        if KratosMultiphysics.ParallelUtilities.GetNumProcs() > 1:
            # Mismatch: the requested manual count is returned.
            self.assertEqual(LinearSolversApplication.MKLUtilities.ComputeMKLThreadCount(2), 2)

    def test_ComputeMKLThreadCount_Consistent(self):
        self._SkipWithoutThreadControl()
        # MKL = 1 <= OMP: consistent, no adjustment.
        LinearSolversApplication.MKLUtilities.SetNumThreads(1)
        self.assertIsNone(LinearSolversApplication.MKLUtilities.ComputeMKLThreadCount(
            LinearSolversApplication.MKLThreadSetting.Consistent))
        if KratosMultiphysics.ParallelUtilities.GetNumProcs() > 1:
            # MKL = 2 > OMP = 1: the OMP thread count is returned.
            KratosMultiphysics.ParallelUtilities.SetNumThreads(1)
            LinearSolversApplication.MKLUtilities.SetNumThreads(2)
            self.assertEqual(LinearSolversApplication.MKLUtilities.ComputeMKLThreadCount(
                LinearSolversApplication.MKLThreadSetting.Consistent), 1)

    def test_ComputeMKLThreadCount_Minimal(self):
        self._SkipWithoutThreadControl()
        # MKL = 1 == min(1, OMP): already minimal, no adjustment.
        LinearSolversApplication.MKLUtilities.SetNumThreads(1)
        self.assertIsNone(LinearSolversApplication.MKLUtilities.ComputeMKLThreadCount(
            LinearSolversApplication.MKLThreadSetting.Minimal))
        if KratosMultiphysics.ParallelUtilities.GetNumProcs() > 1:
            # MKL = 2, OMP = 1: min(2, 1) = 1 is returned.
            KratosMultiphysics.ParallelUtilities.SetNumThreads(1)
            LinearSolversApplication.MKLUtilities.SetNumThreads(2)
            self.assertEqual(LinearSolversApplication.MKLUtilities.ComputeMKLThreadCount(
                LinearSolversApplication.MKLThreadSetting.Minimal), 1)

    def test_ComputeMKLThreadCount_Parameters(self):
        # Explicit and defaulted do_nothing never adjust.
        self.assertIsNone(LinearSolversApplication.MKLUtilities.ComputeMKLThreadCount(
            KratosMultiphysics.Parameters('{"num_threads_mkl": "do_nothing"}')))
        self.assertIsNone(LinearSolversApplication.MKLUtilities.ComputeMKLThreadCount(
            KratosMultiphysics.Parameters('{}')))
        if _HasMKLThreadControl() and KratosMultiphysics.ParallelUtilities.GetNumProcs() > 1:
            LinearSolversApplication.MKLUtilities.SetNumThreads(1)
            self.assertEqual(LinearSolversApplication.MKLUtilities.ComputeMKLThreadCount(
                KratosMultiphysics.Parameters('{"num_threads_mkl": 2}')), 2)
            KratosMultiphysics.ParallelUtilities.SetNumThreads(1)
            LinearSolversApplication.MKLUtilities.SetNumThreads(2)
            self.assertEqual(LinearSolversApplication.MKLUtilities.ComputeMKLThreadCount(
                KratosMultiphysics.Parameters('{"num_threads_mkl": "consistent"}')), 1)
            KratosMultiphysics.ParallelUtilities.SetNumThreads(self.initial_omp_threads)
            LinearSolversApplication.MKLUtilities.SetNumThreads(self.initial_mkl_threads)
        with self.assertRaises(RuntimeError):
            LinearSolversApplication.MKLUtilities.ComputeMKLThreadCount(
                KratosMultiphysics.Parameters('{"num_threads_mkl": "invalid_policy"}'))

    def test_ComputeMKLThreadCount_Invalid(self):
        with self.assertRaises(RuntimeError):
            LinearSolversApplication.MKLUtilities.ComputeMKLThreadCount(-3)


@KratosUnittest.skipIf(not LinearSolversApplication.HasMKL(), "Kratos was compiled without MKL support.")
class TestPardisoMKLThreading(KratosUnittest.TestCase):
    """End-to-end Pardiso solves under explicit MKL/OMP threading configurations.

    These tests guard the scenario of MKL solvers running with default
    multi-threaded MKL on the active (backend-selected) sparse space, and
    exercise the num_threads_mkl solver setting, which no other test covers.
    Note that the num_threads_mkl policy is resolved in the solver's
    Initialize, i.e. at construction time, so the thread state must be
    prepared BEFORE calling ConstructSolver.
    """

    def setUp(self):
        self.initial_mkl_threads = LinearSolversApplication.MKLUtilities.GetNumThreads()
        self.initial_omp_threads = KratosMultiphysics.ParallelUtilities.GetNumThreads()

    def tearDown(self):
        KratosMultiphysics.ParallelUtilities.SetNumThreads(self.initial_omp_threads)
        LinearSolversApplication.MKLUtilities.SetNumThreads(self.initial_mkl_threads)

    def __SolveAndCheck(self, solver_type: str, num_threads_mkl=None) -> None:
        # See the comment in test_eigen_direct_solver.py about the
        # backend-agnostic SparseSpace/SparseMatrix/SparseVector aliases.
        space = KratosMultiphysics.SparseSpace()

        settings_string = '{"solver_type": "LinearSolversApplication.' + solver_type + '"'
        if num_threads_mkl is not None:
            if isinstance(num_threads_mkl, str):
                settings_string += ', "num_threads_mkl": "' + num_threads_mkl + '"'
            else:
                settings_string += ', "num_threads_mkl": ' + str(num_threads_mkl)
        settings_string += '}'
        solver = ConstructSolver(KratosMultiphysics.Parameters(settings_string))

        a = KratosMultiphysics.SparseMatrix()

        base_dir = Path(__file__).resolve().parents[3]
        matrix_file_path = (
            base_dir
            / "kratos"
            / "tests"
            / "auxiliar_files_for_python_unittest"
            / "sparse_matrix_files"
            / "A.mm"
        )

        file_read = KratosMultiphysics.ReadMatrixMarketMatrix(str(matrix_file_path), a)
        self.assertTrue(file_read, msg="The MatrixFile could not be read")

        dimension = a.Size1()
        self.assertEqual(dimension, 900)

        b_exp = KratosMultiphysics.SparseVector(dimension)
        for i in range(dimension):
            b_exp[i] = i + 1

        x = KratosMultiphysics.SparseVector(dimension)

        mkl_threads_before_solve = LinearSolversApplication.MKLUtilities.GetNumThreads()
        solver.Solve(a, x, b_exp)
        self.assertEqual(
            LinearSolversApplication.MKLUtilities.GetNumThreads(),
            mkl_threads_before_solve,
            msg="Solve must restore the MKL thread count it found")

        b_act = KratosMultiphysics.SparseVector(dimension)
        space.Mult(a, x, b_act)
        for i in range(dimension):
            self.assertAlmostEqual(b_act[i], b_exp[i], 7)

    def test_PardisoLU_DefaultThreading(self):
        # Default multi-threaded MKL with the default OMP thread count: the
        # configuration in which MKL solvers must work out of the box.
        self.__SolveAndCheck('pardiso_lu')

    def test_PardisoLU_ManualThreads_One(self):
        self.__SolveAndCheck('pardiso_lu', num_threads_mkl=1)

    def test_PardisoLU_ManualThreads_Max(self):
        self.__SolveAndCheck('pardiso_lu', num_threads_mkl=KratosMultiphysics.ParallelUtilities.GetNumThreads())

    def test_PardisoLU_ConsistentPolicy(self):
        self.__SolveAndCheck('pardiso_lu', num_threads_mkl="consistent")

    def test_PardisoLU_MinimalPolicy(self):
        # With MKL threads <= OMP threads the minimal policy is a no-op.
        if _HasMKLThreadControl():
            LinearSolversApplication.MKLUtilities.SetNumThreads(1)
        self.__SolveAndCheck('pardiso_lu', num_threads_mkl="minimal")

    def test_PardisoLU_MinimalPolicy_MKLThreadsMismatch(self):
        # MKL threads > OMP threads: the minimal policy must clamp the solve
        # to the OMP thread count and restore the MKL count afterwards.
        if not _HasMKLThreadControl():
            self.skipTest("MKL thread count is overridden by the environment (e.g. MKL_NUM_THREADS).")
        if KratosMultiphysics.ParallelUtilities.GetNumProcs() < 2:
            self.skipTest("Requires more than one processor.")
        KratosMultiphysics.ParallelUtilities.SetNumThreads(1)
        LinearSolversApplication.MKLUtilities.SetNumThreads(2)
        self.__SolveAndCheck('pardiso_lu', num_threads_mkl="minimal")

    def test_PardisoLU_OMPThreadVariation(self):
        # Default MKL threading while the Kratos OMP thread count varies.
        for num_threads in (1, KratosMultiphysics.ParallelUtilities.GetNumProcs()):
            with self.subTest(num_threads=num_threads):
                KratosMultiphysics.ParallelUtilities.SetNumThreads(num_threads)
                self.__SolveAndCheck('pardiso_lu')

    def test_PardisoLDLT_ManualThreads(self):
        self.__SolveAndCheck('pardiso_ldlt', num_threads_mkl=1)

    def test_PardisoLLT_ManualThreads(self):
        self.__SolveAndCheck('pardiso_llt', num_threads_mkl=1)


if __name__ == '__main__':
    KratosUnittest.main()
