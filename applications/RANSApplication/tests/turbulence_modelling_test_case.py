import KratosMultiphysics as km
import KratosMultiphysics.KratosUnittest as UnitTest
import KratosMultiphysics.kratos_utilities as kratos_utilities
from test_utilities import RunParametricTestCase


class TurbulenceModellingTestCase(UnitTest.TestCase):
    @classmethod
    def setUpCase(cls, working_folder, parameters_file_name, print_output):
        cls.working_folder = working_folder
        cls.parameters_file_name = parameters_file_name
        cls.print_output = print_output
        cls.parameters = {}

    def testVMSAfcTkeSteady(self):
        self.parameters["<STABILIZATION_METHOD>"] = "algebraic_flux_corrected"
        self.parameters["<WALL_FRICTION_VELOCITY_CALCULATION_METHOD>"] = "turbulent_kinetic_energy_based"
        self.parameters["<TIME_SCHEME_TYPE>"] = "steady"
        self.parameters["<FLOW_SOLVER_FORMULATION>"] = "vms"

        self._runTest()

    def testVMSAfcVelocitySteady(self):
        self.parameters["<STABILIZATION_METHOD>"] = "algebraic_flux_corrected"
        self.parameters["<WALL_FRICTION_VELOCITY_CALCULATION_METHOD>"] = "velocity_based"
        self.parameters["<TIME_SCHEME_TYPE>"] = "steady"
        self.parameters["<FLOW_SOLVER_FORMULATION>"] = "vms"

        self._runTest()

    def testVMSRfcTkeSteady(self):
        self.parameters["<STABILIZATION_METHOD>"] = "residual_based_flux_corrected"
        self.parameters["<WALL_FRICTION_VELOCITY_CALCULATION_METHOD>"] = "turbulent_kinetic_energy_based"
        self.parameters["<TIME_SCHEME_TYPE>"] = "steady"
        self.parameters["<FLOW_SOLVER_FORMULATION>"] = "vms"

        self._runTest()

    def testVMSRfcVelocitySteady(self):
        self.parameters["<STABILIZATION_METHOD>"] = "residual_based_flux_corrected"
        self.parameters["<WALL_FRICTION_VELOCITY_CALCULATION_METHOD>"] = "velocity_based"
        self.parameters["<TIME_SCHEME_TYPE>"] = "steady"
        self.parameters["<FLOW_SOLVER_FORMULATION>"] = "vms"

        self._runTest()

    def testVMSRfcTkeTransient(self):
        self.parameters["<STABILIZATION_METHOD>"] = "residual_based_flux_corrected"
        self.parameters["<WALL_FRICTION_VELOCITY_CALCULATION_METHOD>"] = "turbulent_kinetic_energy_based"
        self.parameters["<TIME_SCHEME_TYPE>"] = self.transient_scheme_type
        self.parameters["<FLOW_SOLVER_FORMULATION>"] = "vms"

        self._runTest()

    def testVMSRfcVelocityTransient(self):
        self.parameters["<STABILIZATION_METHOD>"] = "residual_based_flux_corrected"
        self.parameters["<WALL_FRICTION_VELOCITY_CALCULATION_METHOD>"] = "velocity_based"
        self.parameters["<TIME_SCHEME_TYPE>"] = self.transient_scheme_type
        self.parameters["<FLOW_SOLVER_FORMULATION>"] = "vms"

        self._runTest()

    def testQSVMSRfcVelocitySteady(self):
        if (self.transient_scheme_type == "bdf2"):
            self.skipTest("Skipping since QSVMS is not supported with FractionalStep")

        self.parameters["<STABILIZATION_METHOD>"] = "residual_based_flux_corrected"
        self.parameters["<WALL_FRICTION_VELOCITY_CALCULATION_METHOD>"] = "turbulent_kinetic_energy_based"
        self.parameters["<TIME_SCHEME_TYPE>"] = "steady"
        self.parameters["<FLOW_SOLVER_FORMULATION>"] = "qsvms"

        self._runTest()

    def testQSVMSRfcVelocityTransient(self):
        if (self.transient_scheme_type == "bdf2"):
            self.skipTest("Skipping since QSVMS is not supported with FractionalStep")

        self.parameters["<STABILIZATION_METHOD>"] = "residual_based_flux_corrected"
        self.parameters["<WALL_FRICTION_VELOCITY_CALCULATION_METHOD>"] = "turbulent_kinetic_energy_based"
        self.parameters["<TIME_SCHEME_TYPE>"] = self.transient_scheme_type
        self.parameters["<FLOW_SOLVER_FORMULATION>"] = "qsvms"

        self._runTest()

    def _runTest(self):
        if (km.IsDistributedRun()):
            self.parameters["<PARALLEL_TYPE>"] = "MPI"
        else:
            self.parameters["<PARALLEL_TYPE>"] = "OpenMP"

        self.addCleanup(lambda: kratos_utilities.DeleteTimeFiles("."))

        RunParametricTestCase(self.parameters_file_name, self.working_folder,
                               self.parameters, self.print_output)

