from KratosMultiphysics import IsDistributedRun
import KratosMultiphysics.KratosUnittest as UnitTest
import KratosMultiphysics.kratos_utilities as kratos_utilities
from test_utilities import RunParametricTestCase

class PeriodicTurbulenceModellingTestCase(UnitTest.TestCase):
    @classmethod
    def setUpCase(cls, working_folder, parameters_file_name, print_output):
        cls.working_folder = working_folder
        cls.parameters_file_name = parameters_file_name
        cls.print_output = print_output
        cls.parameters = {}
        cls.parameters["<LINEAR_SOLVER_TYPE>"] = "skyline_lu_factorization"

    def setUp(self):
        if (IsDistributedRun()):
            self.skipTest("Skipping since Periodic tests are not designed to be run in MPI.")

    def testVMSAfcTkeSteady(self):
        self.parameters["<STABILIZATION_METHOD>"] = "algebraic_flux_corrected"
        self.parameters["<WALL_FRICTION_VELOCITY_CALCULATION_METHOD>"] = "turbulent_kinetic_energy_based"
        self.parameters["<FLOW_SOLVER_FORMULATION>"] = "vms"

        self._runTest()

    def testVMSAfcVelocitySteady(self):
        self.parameters["<STABILIZATION_METHOD>"] = "algebraic_flux_corrected"
        self.parameters["<WALL_FRICTION_VELOCITY_CALCULATION_METHOD>"] = "velocity_based"
        self.parameters["<FLOW_SOLVER_FORMULATION>"] = "vms"

        self._runTest()

    def testVMSRfcTkeSteady(self):
        self.parameters["<STABILIZATION_METHOD>"] = "residual_based_flux_corrected"
        self.parameters["<WALL_FRICTION_VELOCITY_CALCULATION_METHOD>"] = "turbulent_kinetic_energy_based"
        self.parameters["<FLOW_SOLVER_FORMULATION>"] = "vms"

        self._runTest()

    def testVMSRfcVelocitySteady(self):
        self.parameters["<STABILIZATION_METHOD>"] = "residual_based_flux_corrected"
        self.parameters["<WALL_FRICTION_VELOCITY_CALCULATION_METHOD>"] = "velocity_based"
        self.parameters["<FLOW_SOLVER_FORMULATION>"] = "vms"

        self._runTest()

    def testQSVMSRfcVelocitySteady(self):
        self.parameters["<STABILIZATION_METHOD>"] = "residual_based_flux_corrected"
        self.parameters["<WALL_FRICTION_VELOCITY_CALCULATION_METHOD>"] = "velocity_based"
        self.parameters["<FLOW_SOLVER_FORMULATION>"] = "qsvms"

        self._runTest()

    def _runTest(self):
        self.addCleanup(lambda: kratos_utilities.DeleteTimeFiles("."))

        RunParametricTestCase(self.parameters_file_name, self.working_folder,
                               self.parameters, self.print_output)

