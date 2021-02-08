import KratosMultiphysics
import KratosMultiphysics.ConvectionDiffusionApplication as convdiff
from KratosMultiphysics.ConvectionDiffusionApplication.response_functions import convection_diffusion_response_function_factory
import KratosMultiphysics.KratosUnittest as KratosUnittest

import KratosMultiphysics.kratos_utilities as kratos_utils


class ConvectionDiffusionResponseFunctionTest(KratosUnittest.TestCase):

    def setUp(self):
        self.work_folder = "adjoint_diffusion_test"

        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            with open(self.file_name,'r') as parameter_file:
                parameters = KratosMultiphysics.Parameters( parameter_file.read())

            # To avoid many prints
            KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)

            model = KratosMultiphysics.Model()
            self.response_function = convection_diffusion_response_function_factory.CreateResponseFunction("dummy", parameters["response_settings"], model)

            # call response function
            self.response_function.Initialize()

    def _calculate_response_and_gradient(self):
        # Within this location context:
        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            self.response_function.InitializeSolutionStep()
            self.response_function.CalculateValue()
            self.value = self.response_function.GetValue()
            self.response_function.CalculateGradient()
            self.gradient = self.response_function.GetNodalGradient(KratosMultiphysics.SHAPE_SENSITIVITY)
            self.response_function.FinalizeSolutionStep()

    def tearDown(self):
        # Within this location context:
        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            self.response_function.Finalize()

            kratos_utils.DeleteFileIfExisting("diffusion_test_primal.post.bin")
            kratos_utils.DeleteFileIfExisting("diffusion_test_adjoint.post.bin")
            kratos_utils.DeleteFileIfExisting("adjoint_diffusion_test.time")
            kratos_utils.DeleteFileIfExisting("adjoint_diffusion_test.post.lst")


class TestAdjointPointTemperatureResponseFunction(ConvectionDiffusionResponseFunctionTest):
    file_name = "AdjointResponse.json"

    def test_execution(self):
        self._calculate_response_and_gradient()
        self.assertAlmostEqual(self.value, 290.0)

        self.assertAlmostEqual(self.gradient[6][0], -3.33333333)
        self.assertAlmostEqual(self.gradient[6][1], -7.84313725)
        self.assertAlmostEqual(self.gradient[6][2], 0.0)


if __name__ == "__main__":
    suites = KratosUnittest.KratosSuites
    smallSuite = suites['small'] # These tests are executed by the continuous integration tool
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestAdjointPointTemperatureResponseFunction]))
    allSuite = suites['all']
    allSuite.addTests(smallSuite)
    KratosUnittest.runTests(suites)
