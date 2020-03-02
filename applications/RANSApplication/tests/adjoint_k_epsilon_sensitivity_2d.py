import os
import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils

from KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable
missing_applications_message = [
    "Missing required application(s):",
]
have_required_applications = CheckIfApplicationsAvailable("HDF5Application")
if not have_required_applications:
    missing_applications_message.append("HDF5Application")

if CheckIfApplicationsAvailable("FluidDynamicsApplication"):
    from KratosMultiphysics.FluidDynamicsApplication.adjoint_fluid_analysis import AdjointFluidAnalysis
else:
    msg = "RANSApplication requires FluidDynamicsApplication which is not found."
    msg += " Please install/compile it and try again."
    raise Exception(msg)

from KratosMultiphysics.RANSApplication.finite_difference_utilities import FiniteDifferenceDragSensitivities


def _removeH5Files(model_part_name):
    for name in os.listdir():
        if name.find(model_part_name) == 0:
            kratos_utils.DeleteFileIfExisting(name)


def _readParameters(parameter_file_name):
    with open(parameter_file_name + '_parameters.json', 'r') as parameter_file:
        project_parameters = Kratos.Parameters(parameter_file.read())
        parameter_file.close()
    return project_parameters


@KratosUnittest.skipUnless(have_required_applications,
                           " ".join(missing_applications_message))
class AdjointKEpsilonSensitivity2D(KratosUnittest.TestCase):
    def setUp(self):
        pass

    def testOneElementSteady(self):
        step_size = 0.00000001

        fd_calculation = FiniteDifferenceDragSensitivities(
            "AdjointKEpsilonSensitivity2DTest/one_element_steady_test.mdpa",
            "AdjointKEpsilonSensitivity2DTest/one_element_steady_test_fd.mdpa",
            "AdjointKEpsilonSensitivity2DTest/one_element_steady_test_parameters.json",
            "MainModelPart.Structure_drag.dat",
            [1.0, 0.0, 0.0],
        )

        fd_sensitivity = fd_calculation.ComputeFiniteDifferenceSensitivity(
            [1], step_size)

        # solve adjoint
        test = AdjointFluidAnalysis(
            Kratos.Model(),
            _readParameters(
                'AdjointKEpsilonSensitivity2DTest/one_element_steady_test_adjoint'
            ))
        test.Run()
        Sensitivity = [[]]
        Sensitivity[0].append(
            test._GetSolver().main_model_part.GetNode(1).GetSolutionStepValue(
                Kratos.SHAPE_SENSITIVITY_X))
        Sensitivity[0].append(
            test._GetSolver().main_model_part.GetNode(1).GetSolutionStepValue(
                Kratos.SHAPE_SENSITIVITY_Y))

        self.assertAlmostEqual(Sensitivity[0][0], fd_sensitivity[0][0], 4)
        self.assertAlmostEqual(Sensitivity[0][1], fd_sensitivity[0][1], 4)
        _removeH5Files("MainModelPart")
        kratos_utils.DeleteFileIfExisting(
            "./AdjointKEpsilonSensitivity2DTest/one_element_test.time")
        kratos_utils.DeleteFileIfExisting(
            "./AdjointKEpsilonSensitivity2DTest/one_element_steady_test_fd.mdpa"
        )
        kratos_utils.DeleteFileIfExisting("./Structure_drag.dat")
        kratos_utils.DeleteFileIfExisting(
            "./finite_difference_drag_sensitivities.dat")
        kratos_utils.DeleteFileIfExisting("./one_element.post.bin")
        kratos_utils.DeleteFileIfExisting("./tests.post.lst")
        kratos_utils.DeleteFileIfExisting("./interface_test.post.bin")

    def tearDown(self):
        pass


if __name__ == '__main__':
    KratosUnittest.main()
