import KratosMultiphysics as km
import KratosMultiphysics.KratosUnittest as UnitTest
import KratosMultiphysics.kratos_utilities as kratos_utilities
from test_utilities import RunParametricTestCase


class FlowSolverTestCase(UnitTest.TestCase):
    @classmethod
    def setUpCase(cls, working_folder, parameters_file_name, materials_file_name, print_output):
        cls.working_folder = working_folder
        cls.parameters_file_name = parameters_file_name
        cls.materials_file_name = materials_file_name
        cls.print_output = print_output
        cls.parameters = {}

    @classmethod
    def tearDownClass(cls):
        with UnitTest.WorkFolderScope(cls.working_folder , __file__):
            kratos_utilities.DeleteTimeFiles(".")
            kratos_utilities.DeleteFileIfExisting(cls.materials_file_name)

    def testSteady(self):
        self.parameters["<TIME_SCHEME_TYPE>"] = "steady"
        self.parameters["<FLOW_SOLVER_FORMULATION>"] = "vms"

        self._runTest()

    def testBossak(self):
        self.parameters["<TIME_SCHEME_TYPE>"] = self.transient_scheme_type
        self.parameters["<FLOW_SOLVER_FORMULATION>"] = "vms"

        self._runTest()

    def _runTest(self):
        if (km.IsDistributedRun()):
            self.parameters["<PARALLEL_TYPE>"] = "MPI"
        else:
            self.parameters["<PARALLEL_TYPE>"] = "OpenMP"

        self.addCleanup(lambda: kratos_utilities.DeleteTimeFiles("."))
        self.addCleanup(lambda: kratos_utilities.DeleteFileIfExisting(self.materials_file_name))

        RunParametricTestCase(self.parameters_file_name, self.materials_file_name, self.working_folder,
                               self.parameters, self.print_output)

