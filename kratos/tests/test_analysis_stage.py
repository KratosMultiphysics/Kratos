import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics

from KratosMultiphysics.analysis_stage import AnalysisStage
from KratosMultiphysics.python_null_solver import PYTHON_NULL_SOLVER


class TestAnalysisStage(KratosUnittest.TestCase):
    """
    Tests for the integration of the Python null solver in AnalysisStage.

    This verifies that:
    - A null solver is used by default
    - A real solver overrides the null solver
    - The null solver behaves as a no-op
    - Solver creation is lazy and cached
    """

    class _AnalysisStageWithDummySolver(AnalysisStage):
        """Stage that provides a dummy solver.

        This is used to test the solver-based behavior.
        """

        class DummySolver:
            """A dummy solver that can be used to test the solver-based behavior."""
            def __init__(self):
                self.called = False

            def AddVariables(self):
                pass

            def SolveSolutionStep(self):
                self.called = True

        def _CreateSolver(self):
            return self.DummySolver()

    def testNullSolver(self):
        """Check that null solver methods do not raise and do nothing."""
        minimal_parameters = KratosMultiphysics.Parameters("""{
            "problem_data": {
                "echo_level": 0,
                "parallel_type": "OpenMP"
            }
        }""")
        stage = AnalysisStage(KratosMultiphysics.Model(), minimal_parameters)
        solver = stage._GetSolver()

        # These should not raise any error as
        solver.Initialize()
        solver.SolveSolutionStep()
        solver.SomeFutureMethod()

        # Check that the solver is the null solver
        self.assertIs(solver, PYTHON_NULL_SOLVER)
        self.assertFalse(stage._HasSolver())

    def testRealSolver(self):
        """Check that the real (dummy) solver is actually executed."""
        minimal_parameters = KratosMultiphysics.Parameters("""{
            "problem_data": {
                "echo_level": 0,
                "parallel_type": "OpenMP"
            }
        }""")
        stage = self._AnalysisStageWithDummySolver(KratosMultiphysics.Model(), minimal_parameters)
        solver = stage._GetSolver()

        # Test dummy solver execution
        solver.SolveSolutionStep()
        self.assertTrue(solver.called)

        # Test that the solver is not the null solver
        self.assertIsNot(solver, PYTHON_NULL_SOLVER)
        self.assertTrue(stage._HasSolver())

    def testSolverIsCached(self):
        """Check that solver is created only once (lazy initialization)."""
        minimal_parameters = KratosMultiphysics.Parameters("""{
            "problem_data": {
                "echo_level": 0,
                "parallel_type": "OpenMP"
            }
        }""")
        stage = self._AnalysisStageWithDummySolver(KratosMultiphysics.Model(), minimal_parameters)

        solver1 = stage._GetSolver()
        solver2 = stage._GetSolver()

        self.assertIs(solver1, solver2)

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()