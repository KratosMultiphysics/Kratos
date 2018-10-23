import os
from KratosMultiphysics import *
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.FluidDynamicsApplication
import KratosMultiphysics.kratos_utilities as kratos_utils

missing_applications_message = ["Missing required application(s):",]
have_required_applications = True

try:
    import KratosMultiphysics.HDF5Application as kh5
except ImportError:
    have_required_applications = False
    missing_applications_message.append("HDF5Application")

from fluid_dynamics_analysis import FluidDynamicsAnalysis
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

    def _removeH5Files(self, model_part_name):
        for name in os.listdir():
            if name.find(model_part_name) == 0:
                kratos_utils.DeleteFileIfExisting(name)

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

    def readDrag(self, filename):
        with open(filename, "r") as file_input:
            lines = file_input.readlines()
        file_input.close()

        found_headers = False

        time_steps = []
        reaction = []

        for line in lines:
            _data = line.strip().split()
            if not found_headers:
                if len(line)>4:
                    if line[:6]=="# Time":
                        found_headers = True
                        continue

            if line.strip()=="":
                continue

            if found_headers:
                _time = float(_data[0])
                _fx = float(_data[1])
                _fy = float(_data[2])
                _fz = float(_data[3])
                time_steps.append(_time)
                reaction.append([_fx, _fy, _fz])


        return time_steps, reaction

    def _getTimeAveragedDrag(self,direction,drag_file_name):
        _time_steps, reaction = self.readDrag(drag_file_name)

        total_drag = 0.0

        for i in range(0, len(_time_steps)):
            _time = _time_steps[i]
            if i==0:
                start_time = _time
                previous_time = start_time
                current_time = previous_time

            total_drag += reaction[i][0]*direction[0]+reaction[i][1]*direction[1]+reaction[i][2]*direction[2]
            previous_time = current_time
            current_time = _time

        if len(_time_steps) > 1:
            delta_t = current_time - previous_time
            total_drag *= (delta_t)

        return total_drag

    def _computeFiniteDifferenceDragSensitivity(self,node_ids,step_size,model_part_file_name,drag_direction,drag_file_name):
        sensitivity = []
        # unperturbed drag
        self.solve(model_part_file_name)
        drag0 = self._getTimeAveragedDrag(drag_direction,drag_file_name)
        for node_id in node_ids:
            node_sensitivity = []
            coord = self._readNodalCoordinates(node_id,model_part_file_name)
            # X + h
            perturbed_coord = [coord[0] + step_size, coord[1], coord[2]]
            self._writeNodalCoordinates(node_id,perturbed_coord,model_part_file_name)
            self.solve(model_part_file_name)
            drag = self._getTimeAveragedDrag(drag_direction,drag_file_name)
            node_sensitivity.append((drag - drag0) / step_size)
            # Y + h
            perturbed_coord = [coord[0], coord[1] + step_size, coord[2]]
            self._writeNodalCoordinates(node_id,perturbed_coord,model_part_file_name)
            self.solve(model_part_file_name)
            drag = self._getTimeAveragedDrag(drag_direction,drag_file_name)
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
            FDSensitivity = self._computeFiniteDifferenceDragSensitivity([1],step_size,'./AdjointVMSSensitivity2DTest/one_element_test',[1.0,0.0,0.0],'./Structure_drag.dat')
            self.assertAlmostEqual(Sensitivity[0][0], FDSensitivity[0][0], 4)
            self.assertAlmostEqual(Sensitivity[0][1], FDSensitivity[0][1], 4)
            self._removeH5Files("MainModelPart")
            kratos_utils.DeleteFileIfExisting("./AdjointVMSSensitivity2DTest/one_element_test.dat")
            kratos_utils.DeleteFileIfExisting("./AdjointVMSSensitivity2DTest/one_element_test.time")
            kratos_utils.DeleteFileIfExisting("./Structure_drag.dat")
            kratos_utils.DeleteFileIfExisting("./one_element.post.bin")
            kratos_utils.DeleteFileIfExisting("./tests.post.lst")

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
            FDSensitivity = self._computeFiniteDifferenceDragSensitivity([1968],step_size,'./AdjointVMSSensitivity2DTest/cylinder_test',[1.0,0.0,0.0],'NoSlip2D_Cylinder_drag.dat')
            self.assertAlmostEqual(Sensitivity[0][0], FDSensitivity[0][0], 5)
            self.assertAlmostEqual(Sensitivity[0][1], FDSensitivity[0][1], 5)
            self._removeH5Files("MainModelPart")
            kratos_utils.DeleteFileIfExisting("./AdjointVMSSensitivity2DTest/cylinder_test.dat")
            kratos_utils.DeleteFileIfExisting("./AdjointVMSSensitivity2DTest/cylinder_test.time")
            kratos_utils.DeleteFileIfExisting("./AdjointVMSSensitivity2DTest/cylinder_test_probe1.dat")
            kratos_utils.DeleteFileIfExisting("./AdjointVMSSensitivity2DTest/cylinder_test_probe2.dat")
            kratos_utils.DeleteFileIfExisting("./AdjointVMSSensitivity2DTest/cylinder_test_adjoint_probe1.dat")
            kratos_utils.DeleteFileIfExisting("./AdjointVMSSensitivity2DTest/cylinder_test_adjoint_probe2.dat")
            kratos_utils.DeleteFileIfExisting("./AdjointVMSSensitivity2DTest/cylinder_test_adjoint_probe3.dat")
            kratos_utils.DeleteFileIfExisting("./NoSlip2D_Cylinder_drag.dat")
            kratos_utils.DeleteFileIfExisting("./cylinder_test.post.bin")

    def testSteadyCylinder(self):
        with ControlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            # solve fluid
            self.solve('AdjointVMSSensitivity2DTest/steady_cylinder_test')
            # solve adjoint
            test = self._createAdjointTest('AdjointVMSSensitivity2DTest/steady_cylinder_test_adjoint')
            test.Run()
            Sensitivity = [[]]
            Sensitivity[0].append(test._GetSolver().main_model_part.GetNode(1968).GetSolutionStepValue(SHAPE_SENSITIVITY_X))
            Sensitivity[0].append(test._GetSolver().main_model_part.GetNode(1968).GetSolutionStepValue(SHAPE_SENSITIVITY_Y))

            # calculate sensitivity by finite difference
            step_size = 0.00000001
            FDSensitivity = self._computeFiniteDifferenceDragSensitivity([1968],step_size,'./AdjointVMSSensitivity2DTest/steady_cylinder_test',[1.0,0.0,0.0],'NoSlip2D_Cylinder_drag.dat')
            self.assertAlmostEqual(Sensitivity[0][0], FDSensitivity[0][0], 4)
            self.assertAlmostEqual(Sensitivity[0][1], FDSensitivity[0][1], 2)
            self._removeH5Files("MainModelPart")
            kratos_utils.DeleteFileIfExisting("./AdjointVMSSensitivity2DTest/steady_cylinder_test.dat")
            kratos_utils.DeleteFileIfExisting("./NoSlip2D_Cylinder_drag.dat")
            kratos_utils.DeleteFileIfExisting("./steady_cylinder_test.time")

    def tearDown(self):
        pass

if __name__ == '__main__':
    KratosUnittest.main()
