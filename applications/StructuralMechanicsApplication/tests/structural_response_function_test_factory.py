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

            model = KratosMultiphysics.Model()
            self.response_function = structural_response_function_factory.CreateResponseFunction("dummy", parameters["kratos_response_settings"], model)

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
        self.assertAlmostEqual(self.value, 0.6062751119154768)

        self.assertAlmostEqual(self.gradient[4][0], -1.1325221565711838)
        self.assertAlmostEqual(self.gradient[4][1], 0.1774595908085034)
        self.assertAlmostEqual(self.gradient[4][2], -0.00012299112980707286)

@KratosUnittest.skipUnless(has_hdf5_application,"Missing required application: HDF5Application")
class TestAdjointDisplacementResponseFunction(StructuralResponseFunctionTestFactory):
    file_name = "adjoint_displacement_response"

    def test_execution(self):
        self._calculate_response_and_gradient()
        model_part = self.response_function.adjoint_analysis.model.GetModelPart("rectangular_plate_structure")
        self.assertAlmostEqual(model_part.Nodes[5].GetSolutionStepValue(ADJOINT_DISPLACEMENT_X), 0.0, 10)
        self.assertAlmostEqual(model_part.Nodes[5].GetSolutionStepValue(ADJOINT_DISPLACEMENT_Y), 0.0, 10)
        self.assertAlmostEqual(model_part.Nodes[5].GetSolutionStepValue(ADJOINT_DISPLACEMENT_Z), 0.012125502238309537)
        self.assertAlmostEqual(model_part.Nodes[4].GetSolutionStepValue(ADJOINT_ROTATION_Y), -0.029186453309188263)

        self.assertAlmostEqual(self.value, 0.12125502238309535)

        self.assertAlmostEqual(self.gradient[4][0], -0.22649691062364785)
        self.assertAlmostEqual(self.gradient[4][1], 0.03549149963347915)
        self.assertAlmostEqual(self.gradient[4][2], -2.4598478882490207e-06, 9)

@KratosUnittest.skipUnless(has_hdf5_application,"Missing required application: HDF5Application")
class TestAdjointStressResponseFunction(StructuralResponseFunctionTestFactory):
    file_name = "adjoint_stress_response"

    def test_execution(self):
        self._calculate_response_and_gradient()

        model_part = self.response_function.adjoint_analysis.model.GetModelPart("rectangular_plate_structure")
        self.assertAlmostEqual(model_part.Nodes[5].GetSolutionStepValue(ADJOINT_DISPLACEMENT_X), 0.0, 10)
        self.assertAlmostEqual(model_part.Nodes[5].GetSolutionStepValue(ADJOINT_DISPLACEMENT_Y), 0.0, 10)
        self.assertAlmostEqual(model_part.Nodes[5].GetSolutionStepValue(ADJOINT_DISPLACEMENT_Z), -0.0823339298948347)
        self.assertAlmostEqual(model_part.Nodes[4].GetSolutionStepValue(ADJOINT_ROTATION_Y), 0.5348048603644553)

        self.assertAlmostEqual(self.value, -0.8233392989483465)

        self.assertAlmostEqual(self.gradient[4][0], 0.3528026402808798)
        self.assertAlmostEqual(self.gradient[4][1], -0.6917210464170941)
        self.assertAlmostEqual(self.gradient[4][2], 0.00011013551613132068)

class TestMassResponseFunction(StructuralResponseFunctionTestFactory):
    file_name = "mass_response"

    def test_execution(self):
        self._calculate_response_and_gradient()
        self.assertAlmostEqual(self.value, 2943.7499999999995)

        self.assertEqual(len(self.gradient.keys()), 9)
        nodeId = 4
        self.assertAlmostEqual(self.gradient[nodeId][0], -1471.874999879219, 5)
        self.assertAlmostEqual(self.gradient[nodeId][1], 0.0)
        self.assertAlmostEqual(self.gradient[nodeId][2], 0.022078165784478184)

class TestStrainEnergyResponseFunction(StructuralResponseFunctionTestFactory):
    file_name = "strain_energy_response"

    def test_execution(self):
        self._calculate_response_and_gradient()
        self.assertAlmostEqual(self.value, 0.00606275111915477)

        nodeId = 4
        self.assertAlmostEqual(self.gradient[nodeId][0], -0.011324859625486555, 9)
        self.assertAlmostEqual(self.gradient[nodeId][1],  0.0017745756664132117, 9)
        self.assertAlmostEqual(self.gradient[nodeId][2], -1.5466215093998671e-07, 9)


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
