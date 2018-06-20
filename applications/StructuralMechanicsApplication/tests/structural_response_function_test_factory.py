# Making KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division
import os

# Import Kratos core and apps
import KratosMultiphysics
from KratosMultiphysics.StructuralMechanicsApplication import *
import KratosMultiphysics.KratosUnittest as KratosUnittest
import structural_response_function_factory

import KratosMultiphysics.kratos_utilities as kratos_utils

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

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

class StructuralResponseFunctionTestFactory(KratosUnittest.TestCase):

    def setUp(self):
        with controlledExecutionScope(self.path):
            print(os.getcwd())
            with open(self.file_name + "_parameters.json",'r') as parameter_file:
                parameters = KratosMultiphysics.Parameters( parameter_file.read())

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
        with controlledExecutionScope(self.path):
            self.response_function.InitializeSolutionStep()
            self.response_function.CalculateValue()
            self.value = self.response_function.GetValue()
            self.response_function.CalculateGradient()
            self.gradient = self.response_function.GetShapeGradient()
            self.response_function.FinalizeSolutionStep()

    def tearDown(self):
        # Within this location context:
        with controlledExecutionScope(self.path):
            self.response_function.Finalize()

            kratos_utils.DeleteFileIfExisting(self.problem_name + ".post.bin")
            kratos_utils.DeleteFileIfExisting(self.problem_name + ".time")
            kratos_utils.DeleteFileIfExisting("response_function_tests.post.lst")

class TestMassResponseFunction(StructuralResponseFunctionTestFactory):
    path = GetFilePath("response_function_tests")
    file_name = "mass_response"

    def test_execution(self):
        self._calculate_response_and_gradient()
        self.assertAlmostEqual(self.value, 1569.9999999999998)

        self.assertEqual(len(self.gradient.keys()), 125)
        self.assertNotEqual(self.gradient[1][0], 0.0)

class TestStrainEnergyResponseFunction(StructuralResponseFunctionTestFactory):
    path = GetFilePath("response_function_tests")
    file_name = "strain_energy_response"

    def test_execution(self):
        self._calculate_response_and_gradient()
        self.assertAlmostEqual(self.value, 8.484005297318718e-05)

        nodeId = 109
        self.assertAlmostEqual(self.gradient[nodeId][0], -1.7336721959838976e-05, 12)
        self.assertAlmostEqual(self.gradient[nodeId][1], 3.6004964666202887e-08, 12)
        self.assertAlmostEqual(self.gradient[nodeId][2], -2.885071090097132e-09, 12)

if __name__ == "__main__":
    suites = KratosUnittest.KratosSuites
    smallSuite = suites['small'] # These tests are executed by the continuous integration tool
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestMassResponseFunction]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestStrainEnergyResponseFunction]))
    allSuite = suites['all']
    allSuite.addTests(smallSuite)
    KratosUnittest.runTests(suites)
