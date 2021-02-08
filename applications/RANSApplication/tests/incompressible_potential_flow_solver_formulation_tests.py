import KratosMultiphysics as km
import KratosMultiphysics.KratosUnittest as UnitTest
import KratosMultiphysics.kratos_utilities as kratos_utilities
from test_utilities import RunParametricTestCase


class IncompressiblePotentialFlowSolverFormulationTest(UnitTest.TestCase):
    def setUp(self):
        # Set to true to get post-process files for the test
        self.print_output = False

    def testSteady(self):
        self._runTest()

    def _runTest(self):
        self.addCleanup(lambda: kratos_utilities.DeleteTimeFiles("."))

        self.parameters = {}
        if (km.IsDistributedRun()):
            self.parameters["<PARALLEL_TYPE>"] = "MPI"
        else:
            self.parameters["<PARALLEL_TYPE>"] = "OpenMP"
        RunParametricTestCase(
            "backward_facing_step_incompressible_potential_flow_parameters.json",
            "BackwardFacingStepTest", self.parameters, self.print_output)


if __name__ == '__main__':
    UnitTest.main()
