import os
from KratosMultiphysics import *
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.FluidDynamicsApplication

missing_applications_message = ["Missing required application(s):",]
have_required_applications = True

try:
    import KratosMultiphysics.AdjointFluidApplication
except ImportError:
    have_required_applications = False
    missing_applications_message.append("AdjointFluidApplication")

try:
    import KratosMultiphysics.HDF5Application as kh5
except ImportError:
    have_required_applications = False
    missing_applications_message.append("HDF5Application")

from fluid_dynamics_analysis import FluidDynamicsAnalysis

if have_required_applications:
    from adjoint_fluid_analysis import AdjointFluidAnalysis

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

    def _removeFile(self, file_path):
        if os.path.isfile(file_path):
            os.remove(file_path)

    def _removeH5Files(self, model_part_name):
        for name in os.listdir():
            if name.find(model_part_name) == 0:
                self._removeFile(name)

    def _readNodalCoordinates(self,node_id,model_part_file_name):
        with open(model_part_file_name + '.mdpa', 'r') as model_part_file:
            lines = model_part_file.readlines()
        lines = lines[lines.index('Begin Nodes\n'):lines.index('End Nodes\n')]
        line = lines[node_id] # assumes consecutive node numbering starting with 1
        components = line.split()
        if int(components[0]) != node_id:
            raise RuntimeError('Error parsing file ' + model_part_file_name)
        return [float(components[i]) for i in range(1,4)]

    def _writeNodalCoordinates(self,node_id,coords,model_part_file_name):
        with open(model_part_file_name + '.mdpa', 'r') as model_part_file:
            lines = model_part_file.readlines()
        node_lines = lines[lines.index('Begin Nodes\n'):lines.index('End Nodes\n')]
        old_line = node_lines[node_id] # assumes consecutive node numbering starting with 1
        components = old_line.split()
        if int(components[0]) != node_id:
            raise RuntimeError('Error parsing file ' + model_part_file_name)
        new_line = '{:5d}'.format(node_id) + ' ' \
             + '{:19.10f}'.format(coords[0]) + ' ' \
             + '{:19.10f}'.format(coords[1]) + ' ' \
             + '{:19.10f}'.format(coords[2]) + '\n'
        lines[lines.index(old_line)] = new_line
        with open(model_part_file_name + '.mdpa', 'w') as model_part_file:
            model_part_file.writelines(lines)

    def _get_time_averaged_drag(self,direction,drag_file_name):
        with open(drag_file_name, 'r') as drag_file:
            lines = drag_file.readlines()
        lines = lines[1:] # remove header
        if len(lines) > 1:
            Dt = float(lines[1].split()[0]) - float(lines[0].split()[0])
        else:
            Dt = 1.0
        dx = 0.0
        dy = 0.0
        dz = 0.0
        for line in lines:
            components = line.split()
            dx = dx + float(components[1])
            dy = dy + float(components[2])
            dz = dz + float(components[3])
        dx = dx * Dt
        dy = dy * Dt
        dz = dz * Dt
        drag = direction[0] * dx + direction[1] * dy + direction[2] * dz
        return drag

    def _computeFiniteDifferenceDragSensitivity(self,node_ids,step_size,model_part_file_name,drag_direction,drag_file_name):
        sensitivity = []
        # unperturbed drag
        self.solve(model_part_file_name)
        drag0 = self._get_time_averaged_drag(drag_direction,drag_file_name)
        for node_id in node_ids:
            node_sensitivity = []
            coord = self._readNodalCoordinates(node_id,model_part_file_name)
            # X + h
            perturbed_coord = [coord[0] + step_size, coord[1], coord[2]]
            self._writeNodalCoordinates(node_id,perturbed_coord,model_part_file_name)
            self.solve(model_part_file_name)
            drag = self._get_time_averaged_drag(drag_direction,drag_file_name)
            node_sensitivity.append((drag - drag0) / step_size)
            # Y + h
            perturbed_coord = [coord[0], coord[1] + step_size, coord[2]]
            self._writeNodalCoordinates(node_id,perturbed_coord,model_part_file_name)
            self.solve(model_part_file_name)
            drag = self._get_time_averaged_drag(drag_direction,drag_file_name)
            node_sensitivity.append((drag - drag0) / step_size)
            sensitivity.append(node_sensitivity)
            # return mdpa file to unperturbed state
            self._writeNodalCoordinates(node_id,coord,model_part_file_name)
        return sensitivity

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

    def solve(self, parameter_file_name):
        test = self._createFluidTest(parameter_file_name)
        test.Run()

    def testOneElement(self):
        with ControlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            # solve fluid
            self.solve('AdjointVMSSensitivity2DTest/one_element_test')
            # solve adjoint
            test = AdjointFluidAnalysis(Model(), self._readParameters('AdjointVMSSensitivity2DTest/one_element_test_adjoint'))
            test.Run()
            Sensitivity = [[]]
            Sensitivity[0].append(test._GetSolver().main_model_part.GetNode(1).GetSolutionStepValue(SHAPE_SENSITIVITY_X))
            Sensitivity[0].append(test._GetSolver().main_model_part.GetNode(1).GetSolutionStepValue(SHAPE_SENSITIVITY_Y))

            # calculate sensitivity by finite difference
            step_size = 0.00000001
            FDSensitivity = self._computeFiniteDifferenceDragSensitivity([1],step_size,'./AdjointVMSSensitivity2DTest/one_element_test',[1.0,0.0,0.0],'./AdjointVMSSensitivity2DTest/one_element_test.dat')
            self.assertAlmostEqual(Sensitivity[0][0], FDSensitivity[0][0], 4)
            self.assertAlmostEqual(Sensitivity[0][1], FDSensitivity[0][1], 4)
            self._removeH5Files("MainModelPart")
            self._removeFile("./AdjointVMSSensitivity2DTest/one_element_test.dat")
            self._removeFile("./AdjointVMSSensitivity2DTest/one_element_test.time")
            self._removeFile("./one_element.post.bin")
            self._removeFile("./tests.post.lst")

    def testCylinder(self):
        with ControlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            # solve fluid
            self.solve('AdjointVMSSensitivity2DTest/cylinder_test')
            # solve adjoint
            test = self._createAdjointTest('AdjointVMSSensitivity2DTest/cylinder_test_adjoint')
            test.Run()
            Sensitivity = [[]]
            Sensitivity[0].append(test._GetSolver().main_model_part.GetNode(1968).GetSolutionStepValue(SHAPE_SENSITIVITY_X))
            Sensitivity[0].append(test._GetSolver().main_model_part.GetNode(1968).GetSolutionStepValue(SHAPE_SENSITIVITY_Y))

            # calculate sensitivity by finite difference
            step_size = 0.00000001
            FDSensitivity = self._computeFiniteDifferenceDragSensitivity([1968],step_size,'./AdjointVMSSensitivity2DTest/cylinder_test',[1.0,0.0,0.0],'./AdjointVMSSensitivity2DTest/cylinder_test.dat')
            self.assertAlmostEqual(Sensitivity[0][0], FDSensitivity[0][0], 5)
            self.assertAlmostEqual(Sensitivity[0][1], FDSensitivity[0][1], 5)
            self._removeH5Files("MainModelPart")
            self._removeFile("./AdjointVMSSensitivity2DTest/cylinder_test.dat")
            self._removeFile("./AdjointVMSSensitivity2DTest/cylinder_test.time")
            self._removeFile("./AdjointVMSSensitivity2DTest/cylinder_test_probe1.dat")
            self._removeFile("./AdjointVMSSensitivity2DTest/cylinder_test_probe2.dat")
            self._removeFile("./AdjointVMSSensitivity2DTest/cylinder_test_adjoint_probe1.dat")
            self._removeFile("./AdjointVMSSensitivity2DTest/cylinder_test_adjoint_probe2.dat")
            self._removeFile("./AdjointVMSSensitivity2DTest/cylinder_test_adjoint_probe3.dat")
            self._removeFile("./cylinder_test.post.bin")

    def tearDown(self):
        pass

if __name__ == '__main__':
    KratosUnittest.main()
