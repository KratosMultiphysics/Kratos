# Making KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division
import os

# Import Kratos core and apps
import KratosMultiphysics
from KratosMultiphysics.StructuralMechanicsApplication import *
import KratosMultiphysics.KratosUnittest as KratosUnittest
import structural_response_function_factory

import KratosMultiphysics.kratos_utilities as kratos_utils

try:
    from KratosMultiphysics.HDF5Application import *
    has_hdf5_application = True
except ImportError:
    has_hdf5_application = False

# This utility will control the execution scope in case we need to access files or we depend
# on specific relative locations of the files.

# TODO: Should we move this to KratosUnittest?
class controlledExecutionScope:
    def __init__(self, scope):
        self.currentPath = os.getcwd()
        self.scope = scope

    def __enter__(self):
        os.chdir(self.scope)

    def __exit__(self, type, value, traceback):
        os.chdir(self.currentPath)

def _get_test_working_dir():
    this_file_dir = os.path.dirname(os.path.realpath(__file__))
    return os.path.join(this_file_dir, "response_function_tests")

class StructuralResponseFunctionTestFactory(KratosUnittest.TestCase):

    def setUp(self):
        with controlledExecutionScope(_get_test_working_dir()):
            with open(self.file_name + "_parameters.json",'r') as parameter_file:
                parameters = KratosMultiphysics.Parameters( parameter_file.read())

            # To avoid many prints
            if (parameters["problem_data"]["echo_level"].GetInt() == 0):
                KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)

            self.problem_name = parameters["problem_data"]["problem_name"].GetString()
            model_part = KratosMultiphysics.ModelPart(self.problem_name)

            self.response_function = structural_response_function_factory.CreateResponseFunction("dummy", parameters["kratos_response_settings"], model_part)

            # import model part
            model_part_io = KratosMultiphysics.ModelPartIO(parameters["problem_data"]["problem_name"].GetString())
            model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 3)
            model_part_io.ReadModelPart(model_part)

            # call response function
            self.response_function.Initialize()

    def _calculate_response_and_gradient(self):
        # Within this location context:
        with controlledExecutionScope(_get_test_working_dir()):
            self.response_function.InitializeSolutionStep()
            self.response_function.CalculateValue()
            self.value = self.response_function.GetValue()
            self.response_function.CalculateGradient()
            self.gradient = self.response_function.GetShapeGradient()
            self.response_function.FinalizeSolutionStep()

    def tearDown(self):
        # Within this location context:
        with controlledExecutionScope(_get_test_working_dir()):
            self.response_function.Finalize()

            kratos_utils.DeleteFileIfExisting(self.problem_name + ".post.bin")
            kratos_utils.DeleteFileIfExisting(self.problem_name + ".time")
            kratos_utils.DeleteFileIfExisting(self.problem_name + ".h5")
            kratos_utils.DeleteFileIfExisting(self.problem_name + "-1.0000.h5")
            kratos_utils.DeleteFileIfExisting("response_function_tests.post.lst")

@KratosUnittest.skipUnless(has_hdf5_application,"Missing required application: HDF5Application")
class TestAdjointStrainEnergyResponseFunction(StructuralResponseFunctionTestFactory):
    file_name = "adjoint_strain_energy_response"

    def test_execution(self):
        self._calculate_response_and_gradient()
        self.assertAlmostEqual(self.value, 0.011841310779599155)

        nodeId = 4
        self.assertAlmostEqual(self.gradient[nodeId][0], -0.02211883892695645, 12)
        self.assertAlmostEqual(self.gradient[nodeId][1], 0.0034659667625643337, 12)
        self.assertAlmostEqual(self.gradient[nodeId][2], -3.942528133752394e-07, 12)

@KratosUnittest.skipUnless(has_hdf5_application,"Missing required application: HDF5Application")
class TestAdjointDisplacementResponseFunction(StructuralResponseFunctionTestFactory):
    file_name = "adjoint_displacement_response"

    def test_execution(self):
        self._calculate_response_and_gradient()
        self.assertAlmostEqual(self.value, 0.0023682621559198303)

        nodeId = 4
        self.assertAlmostEqual(self.gradient[nodeId][0], 35162269878424.273, 12)
        self.assertAlmostEqual(self.gradient[nodeId][1], 4507732710819.046, 12)
        self.assertAlmostEqual(self.gradient[nodeId][2], 1125096436.1808963, 12)

@KratosUnittest.skipUnless(has_hdf5_application,"Missing required application: HDF5Application")
class TestAdjointStressResponseFunction(StructuralResponseFunctionTestFactory):
    file_name = "adjoint_stress_response"

    def test_execution(self):
        self._calculate_response_and_gradient()
        self.assertAlmostEqual(self.value, -0.8233392989483465)

        nodeId = 4
        self.assertAlmostEqual(self.gradient[nodeId][0], -8.480870212810133, 12)
        self.assertAlmostEqual(self.gradient[nodeId][1], -2.6870951672193213, 12)
        self.assertAlmostEqual(self.gradient[nodeId][2], -0.0002913758359272386, 12)

class TestMassResponseFunction(StructuralResponseFunctionTestFactory):
    file_name = "mass_response"

    def test_execution(self):
        self._calculate_response_and_gradient()
        self.assertAlmostEqual(self.value, 2355.0)

        self.assertEqual(len(self.gradient.keys()), 9)
        nodeId = 4
        self.assertAlmostEqual(self.gradient[nodeId][0], -1177.5000000397995, 12)
        self.assertAlmostEqual(self.gradient[nodeId][1], 2.2737367544323206e-07, 12)
        self.assertAlmostEqual(self.gradient[nodeId][2], 0.0017662387108430266, 12)

class TestStrainEnergyResponseFunction(StructuralResponseFunctionTestFactory):
    file_name = "strain_energy_response"

    def test_execution(self):
        self._calculate_response_and_gradient()
        self.assertAlmostEqual(self.value, 0.00011841310779599155)

        nodeId = 4
        self.assertAlmostEqual(self.gradient[nodeId][0], -0.00022118866455491083, 12)
        self.assertAlmostEqual(self.gradient[nodeId][1],  0.00003465968099176853, 12)
        self.assertAlmostEqual(self.gradient[nodeId][2], -4.974301305577504e-09, 12)


if __name__ == "__main__":
    suites = KratosUnittest.KratosSuites
    smallSuite = suites['small'] # These tests are executed by the continuous integration tool
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestAdjointStrainEnergyResponseFunction]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestAdjointDisplacementResponseFunction]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestAdjointStressResponseFunction]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestMassResponseFunction]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestStrainEnergyResponseFunction]))
    allSuite = suites['all']
    allSuite.addTests(smallSuite)
    KratosUnittest.runTests(suites)
