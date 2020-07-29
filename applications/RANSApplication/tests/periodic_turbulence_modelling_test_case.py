import KratosMultiphysics.KratosUnittest as UnitTest
import KratosMultiphysics.kratos_utilities as kratos_utilities
from KratosMultiphysics.RANSApplication.test_utilities import RunParametericTestCase


class PeriodicTurbulenceModellingTestCase(UnitTest.TestCase):
    @classmethod
    def setUpCase(cls, working_folder, parameters_file_name, print_output):
        cls.working_folder = working_folder
        cls.parameters_file_name = parameters_file_name
        cls.print_output = print_output
        cls.parameters = {}
        cls.parameters["<LINEAR_SOLVER_TYPE>"] = "skyline_lu_factorization"

    def testAfcTkeSteady(self):
        self.parameters["<STABILIZATION_METHOD>"] = "algebraic_flux_corrected"
        self.parameters["<WALL_FRICTION_VELOCITY_CALCULATION_METHOD>"] = "turbulent_kinetic_energy_based"

        self._runTest()

    def testAfcVelocitySteady(self):
        self.parameters["<STABILIZATION_METHOD>"] = "algebraic_flux_corrected"
        self.parameters["<WALL_FRICTION_VELOCITY_CALCULATION_METHOD>"] = "velocity_based"

        self._runTest()

    def testRfcTkeSteady(self):
        self.parameters["<STABILIZATION_METHOD>"] = "residual_based_flux_corrected"
        self.parameters["<WALL_FRICTION_VELOCITY_CALCULATION_METHOD>"] = "turbulent_kinetic_energy_based"

        self._runTest()

    def testRfcVelocitySteady(self):
        self.parameters["<STABILIZATION_METHOD>"] = "residual_based_flux_corrected"
        self.parameters["<WALL_FRICTION_VELOCITY_CALCULATION_METHOD>"] = "velocity_based"

        self._runTest()

    def _runTest(self):
        self.addCleanup(lambda: kratos_utilities.DeleteTimeFiles("."))

        RunParametericTestCase(self.parameters_file_name, self.working_folder,
                               self.parameters, self.print_output)

