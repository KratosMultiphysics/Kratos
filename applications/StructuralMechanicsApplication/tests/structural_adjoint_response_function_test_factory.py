# Making KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division
import os

# Import Kratos core and apps
import KratosMultiphysics
from KratosMultiphysics.StructuralMechanicsApplication import *
import KratosMultiphysics.KratosUnittest as KratosUnittest
import adjoint_structural_response_function_factory

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
            self.response_function = adjoint_structural_response_function_factory.CreateAdjointResponseFunction("dummy", parameters["kratos_response_settings"], model)

            # call response function
            self.response_function.Initialize()

    def _calculate_response_and_gradient(self):
        # Within this location context:
        with controlledExecutionScope(_get_test_working_dir()):
            self.response_function.InitializeSolutionStep()
            self.response_function.SolveAdjointProblem()
            self.response_function.PerformSensitivityPostprocess()
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
class TestAdjointDisplacementResponseFunction(StructuralResponseFunctionTestFactory):
    file_name = "adjoint_displacement_response"

    def test_execution(self):
        self._calculate_response_and_gradient()
        model_part = self.response_function.adjoint_analysis.model.GetModelPart("rectangular_plate_structure")
        self.assertAlmostEqual(model_part.Nodes[5].GetSolutionStepValue(ADJOINT_DISPLACEMENT_X), 0.0, 10)
        self.assertAlmostEqual(model_part.Nodes[5].GetSolutionStepValue(ADJOINT_DISPLACEMENT_Y), 0.0, 10)
        self.assertAlmostEqual(model_part.Nodes[5].GetSolutionStepValue(ADJOINT_DISPLACEMENT_Z), 0.012125502238309537)
        self.assertAlmostEqual(model_part.Nodes[4].GetSolutionStepValue(ADJOINT_ROTATION_Y), -0.029186453309188263)


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


if __name__ == "__main__":
    suites = KratosUnittest.KratosSuites
    smallSuite = suites['small'] # These tests are executed by the continuous integration tool
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestAdjointDisplacementResponseFunction]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestAdjointStressResponseFunction]))
    allSuite = suites['all']
    allSuite.addTests(smallSuite)
    KratosUnittest.runTests(suites)
