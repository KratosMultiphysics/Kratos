import os
from KratosMultiphysics import *
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.FluidDynamicsApplication
import KratosMultiphysics.kratos_utilities as kratos_utils

from  KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable
missing_applications_message = ["Missing required application(s):",]
have_required_applications = CheckIfApplicationsAvailable("HDF5Application")
if have_required_applications:
    import KratosMultiphysics.HDF5Application as kh5
else:
    missing_applications_message.append("HDF5Application")

from KratosMultiphysics.FluidDynamicsApplication.fluid_dynamics_analysis import FluidDynamicsAnalysis
from KratosMultiphysics.FluidDynamicsApplication.adjoint_fluid_analysis import AdjointFluidAnalysis
from KratosMultiphysics.FluidDynamicsApplication.finite_difference_drag_shape_sensitivity_analysis import FiniteDifferenceDragShapeSensitivityAnalysis

class ControlledExecutionScope:
    def __init__(self, scope):
        self.currentPath = os.getcwd()
        self.scope = scope

    def __enter__(self):
        os.chdir(self.scope)

    def __exit__(self, type, value, traceback):
        os.chdir(self.currentPath)

@KratosUnittest.skipUnless(have_required_applications," ".join(missing_applications_message))
class AdjointVMSSensitivity2D(KratosUnittest.TestCase):
    def setUp(self):
        pass

    def _removeH5Files(self, model_part_name):
        for name in os.listdir():
            if name.find(model_part_name) == 0:
                kratos_utils.DeleteFileIfExisting(name)

    def _readParameters(self, parameter_file_name):
        with open(parameter_file_name + '_parameters.json', 'r') as parameter_file:
            project_parameters = Parameters(parameter_file.read())
            parameter_file.close()
        return project_parameters

    def _createFluidTest(self, parameter_file_name):
        test = FluidDynamicsAnalysis(Model(), self._readParameters(parameter_file_name))
        return test

    def _createAdjointTest(self, parameter_file_name):
        test = AdjointFluidAnalysis(Model(), self._readParameters(parameter_file_name))
        return test

    def testOneElement(self):
        with ControlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            # calculate sensitivity by finite difference
            step_size = 0.00000001
            fd_analysis = FiniteDifferenceDragShapeSensitivityAnalysis(2, [1], './AdjointVMSSensitivity2DTest/one_element_test', lambda : self._createFluidTest("AdjointVMSSensitivity2DTest/one_element_test"))
            FDSensitivity = fd_analysis.ComputeSensitivity([1.0, 0.0, 0.0], './MainModelPart.Structure_drag.dat', step_size)

            # calculate adjoint sensitivity
            test = AdjointFluidAnalysis(Model(), self._readParameters('AdjointVMSSensitivity2DTest/one_element_test_adjoint'))
            test.Run()
            Sensitivity = [[]]
            Sensitivity[0].append(test._GetSolver().main_model_part.GetNode(1).GetSolutionStepValue(SHAPE_SENSITIVITY_X))
            Sensitivity[0].append(test._GetSolver().main_model_part.GetNode(1).GetSolutionStepValue(SHAPE_SENSITIVITY_Y))

            self.assertAlmostEqual(Sensitivity[0][0], FDSensitivity[0][0], 4)
            self.assertAlmostEqual(Sensitivity[0][1], FDSensitivity[0][1], 4)

            # cleanup
            kratos_utils.DeleteFileIfExisting("./AdjointVMSSensitivity2DTest/one_element_test.time")
            kratos_utils.DeleteFileIfExisting("./Structure_drag.dat")
            kratos_utils.DeleteFileIfExisting("./one_element.post.bin")

    def testCylinder(self):
        with ControlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            # calculate sensitivity by finite difference
            step_size = 0.00000001
            fd_analysis = FiniteDifferenceDragShapeSensitivityAnalysis(2, [1968], './AdjointVMSSensitivity2DTest/cylinder_test', lambda : self._createFluidTest("AdjointVMSSensitivity2DTest/cylinder_test"))
            FDSensitivity = fd_analysis.ComputeSensitivity([1.0, 0.0, 0.0], 'MainModelPart.NoSlip2D_Cylinder_drag.dat', step_size)

            # calculate adjoint sensitivity
            test = self._createAdjointTest('AdjointVMSSensitivity2DTest/cylinder_test_adjoint')
            test.Run()
            Sensitivity = [[]]
            Sensitivity[0].append(test._GetSolver().main_model_part.GetNode(1968).GetSolutionStepValue(SHAPE_SENSITIVITY_X))
            Sensitivity[0].append(test._GetSolver().main_model_part.GetNode(1968).GetSolutionStepValue(SHAPE_SENSITIVITY_Y))

            self.assertAlmostEqual(Sensitivity[0][0], FDSensitivity[0][0], 5)
            self.assertAlmostEqual(Sensitivity[0][1], FDSensitivity[0][1], 5)

            # cleanup
            kratos_utils.DeleteFileIfExisting("./AdjointVMSSensitivity2DTest/cylinder_test.time")
            kratos_utils.DeleteFileIfExisting("./NoSlip2D_Cylinder_drag.dat")
            kratos_utils.DeleteFileIfExisting("./cylinder_test.post.bin")

    def testSteadyCylinder(self):
        with ControlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            # calculate sensitivity by finite difference
            step_size = 0.00000001
            fd_analysis = FiniteDifferenceDragShapeSensitivityAnalysis(2, [1968], './AdjointVMSSensitivity2DTest/steady_cylinder_test', lambda : self._createFluidTest("AdjointVMSSensitivity2DTest/steady_cylinder_test"))
            FDSensitivity = fd_analysis.ComputeSensitivity([1.0, 0.0, 0.0], 'MainModelPart.NoSlip2D_Cylinder_drag.dat', step_size)

            # calculate adjoint sensitivity
            test = self._createAdjointTest('AdjointVMSSensitivity2DTest/steady_cylinder_test_adjoint')
            test.Run()
            Sensitivity = [[]]
            Sensitivity[0].append(test._GetSolver().main_model_part.GetNode(1968).GetSolutionStepValue(SHAPE_SENSITIVITY_X))
            Sensitivity[0].append(test._GetSolver().main_model_part.GetNode(1968).GetSolutionStepValue(SHAPE_SENSITIVITY_Y))

            self.assertAlmostEqual(Sensitivity[0][0], FDSensitivity[0][0], 4)
            self.assertAlmostEqual(Sensitivity[0][1], FDSensitivity[0][1], 2)

            # cleanup
            kratos_utils.DeleteFileIfExisting("./AdjointVMSSensitivity2DTest/steady_cylinder_test.time")
            kratos_utils.DeleteFileIfExisting("./NoSlip2D_Cylinder_drag.dat")
            kratos_utils.DeleteFileIfExisting("./steady_cylinder_test.post.bin")

    def tearDown(self):
        self._removeH5Files("MainModelPart")
        kratos_utils.DeleteFileIfExisting("./tests.post.lst")

if __name__ == '__main__':
    KratosUnittest.main()
