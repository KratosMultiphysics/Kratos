# Making KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division
import os

# Import Kratos core and apps
import KratosMultiphysics
from KratosMultiphysics.StructuralMechanicsApplication import *
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.StructuralMechanicsApplication import structural_response_function_factory

import KratosMultiphysics.kratos_utilities as kratos_utils

if kratos_utils.CheckIfApplicationsAvailable("EigenSolversApplication"):
    has_eigensolvers_application = True
else:
    has_eigensolvers_application = False

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
            self.response_function = structural_response_function_factory.CreateResponseFunction("dummy", parameters["response_settings"], model)

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

class TestAdjointStrainEnergyResponseFunction(StructuralResponseFunctionTestFactory):
    file_name = "adjoint_strain_energy_response"

    def test_execution(self):
        self._calculate_response_and_gradient()
        self.assertAlmostEqual(self.value, 0.6062751119154768)

        self.assertAlmostEqual(self.gradient[4][0], -1.1325221565711838)
        self.assertAlmostEqual(self.gradient[4][1], 0.1774595908085034)
        self.assertAlmostEqual(self.gradient[4][2], -0.00012299112980707286)

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

class TestAdjointMaxStressResponseFunction(StructuralResponseFunctionTestFactory):
    file_name = "adjoint_max_stress_response"

    def test_execution(self):
        self._calculate_response_and_gradient()

        model_part = self.response_function.adjoint_analysis.model.GetModelPart("cantilever_beam")
        self.assertAlmostEqual(model_part.Nodes[53].GetSolutionStepValue(ADJOINT_DISPLACEMENT_X), 7.657448651571309, 10)
        self.assertAlmostEqual(model_part.Nodes[53].GetSolutionStepValue(ADJOINT_DISPLACEMENT_Y), -19.9044491754745, 10)
        self.assertAlmostEqual(model_part.Nodes[53].GetSolutionStepValue(ADJOINT_DISPLACEMENT_Z), -8.37326311561973, 10)

        self.assertIsClose(self.value, 1610060.3904999627)

        self.assertIsClose(self.gradient[5][0], 1787255.3702425747)
        self.assertIsClose(self.gradient[5][1], -247.0446103799622, rel_tol=1e-5)
        self.assertIsClose(self.gradient[5][2], -562640.0306970887)

class TestMassResponseFunction(StructuralResponseFunctionTestFactory):
    file_name = "mass_response"

    def test_execution(self):
        self.current_model = KratosMultiphysics.Model()
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
        self.current_model = KratosMultiphysics.Model()
        self._calculate_response_and_gradient()
        self.assertAlmostEqual(self.value, 0.6062751119154768)

        self.assertAlmostEqual(self.gradient[4][0], -1.132485962510845)
        self.assertAlmostEqual(self.gradient[4][1], 0.17745756668175833)
        self.assertAlmostEqual(self.gradient[4][2], -1.5466170818541692e-05)

@KratosUnittest.skipUnless(has_eigensolvers_application,"Missing required application: EigenSolversApplication")
class TestEigenfrequencyResponseFunction(StructuralResponseFunctionTestFactory):
    file_name = "eigenfrequency_response"

    def test_execution(self):
        self._calculate_response_and_gradient()

        self.assertAlmostEqual(self.value, 0.014123803835107267)

        self.assertAlmostEqual(self.gradient[19][0], -2.629125898327006e-09)
        self.assertAlmostEqual(self.gradient[19][1], -0.0070165270383214864)
        self.assertAlmostEqual(self.gradient[19][2], 1.0549911195848093e-08)

if __name__ == "__main__":
    suites = KratosUnittest.KratosSuites
    smallSuite = suites['small'] # These tests are executed by the continuous integration tool
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestAdjointStrainEnergyResponseFunction]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestAdjointDisplacementResponseFunction]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestAdjointStressResponseFunction]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestAdjointMaxStressResponseFunction]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestMassResponseFunction]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestStrainEnergyResponseFunction]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestEigenfrequencyResponseFunction]))
    allSuite = suites['all']
    allSuite.addTests(smallSuite)
    KratosUnittest.runTests(suites)
