import os
from KratosMultiphysics import *
import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_MainKratos

class ControlledExecutionScope:
    def __init__(self, scope):
        self.currentPath = os.getcwd()
        self.scope = scope

    def __enter__(self):
        os.chdir(self.scope)

    def __exit__(self, type, value, traceback):
        os.chdir(self.currentPath)

class TestCase(KratosUnittest.TestCase):

    def setUp(self):
        pass

    def readNodalCoordinates(self,node_id,model_part_file_name):
        with open(model_part_file_name + '.mdpa', 'r') as model_part_file:
            lines = model_part_file.readlines()
        lines = lines[lines.index('Begin Nodes\n'):lines.index('End Nodes\n')]
        line = lines[node_id] # assumes consecutive node numbering starting with 1
        components = line.split()
        if int(components[0]) != node_id:
            raise RuntimeError('Error parsing file ' + model_part_file_name)
        return [float(components[i]) for i in range(1,4)]

    def writeNodalCoordinates(self,node_id,coords,model_part_file_name):
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

    def getTimeAveragedDrag(self,direction,drag_file_name):
        with open(drag_file_name, 'r') as drag_file:
            lines = drag_file.readlines()
        lines = lines[1:] # remove header
        if len(lines) > 1:
            Dt = float(lines[1].split()[0]) - float(lines[0].split()[0])
            T  = float(lines[-1].split()[0]) - float(lines[0].split()[0]) + Dt
        else:
            Dt = 1.0
            T = 1.0
        dx = 0.0
        dy = 0.0
        dz = 0.0
        for line in lines:
            components = line.split()
            dx = dx + float(components[1])
            dy = dy + float(components[2])
            dz = dz + float(components[3])
        dx = dx * Dt / T
        dy = dy * Dt / T
        dz = dz * Dt / T
        drag = direction[0] * dx + direction[1] * dy + direction[2] * dz
        return drag

    def ComputeFiniteDifferenceDragSensitivity(self,node_ids,step_size,model_part_file_name,drag_direction,drag_file_name):
        sensitivity = []
        # unperturbed drag
        self.solve(model_part_file_name)
        drag0 = self.getTimeAveragedDrag(drag_direction,drag_file_name)
        for node_id in node_ids:
            node_sensitivity = []
            coord = self.readNodalCoordinates(node_id,model_part_file_name)
            # X + h
            perturbed_coord = [coord[0] + step_size, coord[1], coord[2]]
            self.writeNodalCoordinates(node_id,perturbed_coord,model_part_file_name)
            self.solve(model_part_file_name)
            drag = self.getTimeAveragedDrag(drag_direction,drag_file_name)
            node_sensitivity.append((drag - drag0) / step_size)
            # Y + h
            perturbed_coord = [coord[0], coord[1] + step_size, coord[2]]
            self.writeNodalCoordinates(node_id,perturbed_coord,model_part_file_name)
            self.solve(model_part_file_name)
            drag = self.getTimeAveragedDrag(drag_direction,drag_file_name)
            node_sensitivity.append((drag - drag0) / step_size)
            sensitivity.append(node_sensitivity)
            # return mdpa file to unperturbed state
            self.writeNodalCoordinates(node_id,coord,model_part_file_name)
        return sensitivity

    def createTest(self, parameter_file_name):
        with open(parameter_file_name + '_parameters.json', 'r') as parameter_file:
            project_parameters = Parameters(parameter_file.read())
            parameter_file.close()
        test = test_MainKratos.MainKratos(project_parameters)
        return test

    def solve(self, parameter_file_name):
        test = self.createTest(parameter_file_name)
        test.Solve()

    def test_OneElement(self):
        with ControlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            # solve fluid
            self.solve('test_vms_sensitivity_2d/one_element_test')
            # solve adjoint
            test = self.createTest('test_vms_sensitivity_2d/one_element_test_adjoint')
            test.Solve()
            Sensitivity = [[]]
            Sensitivity[0].append(test.main_model_part.GetNode(1).GetSolutionStepValue(SHAPE_SENSITIVITY_X))
            Sensitivity[0].append(test.main_model_part.GetNode(1).GetSolutionStepValue(SHAPE_SENSITIVITY_Y))

            # calculate sensitivity by finite difference
            step_size = 0.00000001
            FDSensitivity = self.ComputeFiniteDifferenceDragSensitivity([1],step_size,'./test_vms_sensitivity_2d/one_element_test',[1.0,0.0,0.0],'./test_vms_sensitivity_2d/one_element_test.dat')
            self.assertAlmostEqual(Sensitivity[0][0], FDSensitivity[0][0], 5)
            self.assertAlmostEqual(Sensitivity[0][1], FDSensitivity[0][1], 5)
            # remove hdf5 file
            if "one_element_test_0.h5" in os.listdir("./test_vms_sensitivity_2d"):
                os.remove("./test_vms_sensitivity_2d/one_element_test_0.h5")
            # remove drag file
            if "one_element_test.dat" in os.listdir("./test_vms_sensitivity_2d"):
                os.remove("./test_vms_sensitivity_2d/one_element_test.dat")

    def test_Cylinder(self):
        with ControlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            # solve fluid
            self.solve('test_vms_sensitivity_2d/cylinder_test')
            # solve adjoint
            test = self.createTest('test_vms_sensitivity_2d/cylinder_test_adjoint')
            test.Solve()
            Sensitivity = [[]]
            Sensitivity[0].append(test.main_model_part.GetNode(1968).GetSolutionStepValue(SHAPE_SENSITIVITY_X))
            Sensitivity[0].append(test.main_model_part.GetNode(1968).GetSolutionStepValue(SHAPE_SENSITIVITY_Y))

            # calculate sensitivity by finite difference
            step_size = 0.00000001
            FDSensitivity = self.ComputeFiniteDifferenceDragSensitivity([1968],step_size,'./test_vms_sensitivity_2d/cylinder_test',[1.0,0.0,0.0],'./test_vms_sensitivity_2d/cylinder_test.dat')
            self.assertAlmostEqual(Sensitivity[0][0], FDSensitivity[0][0], 5)
            self.assertAlmostEqual(Sensitivity[0][1], FDSensitivity[0][1], 5)
            # remove hdf5 file
            if "cylinder_test_0.h5" in os.listdir("./test_vms_sensitivity_2d"):
                os.remove("./test_vms_sensitivity_2d/cylinder_test_0.h5")
            # remove drag file
            if "cylinder_test.dat" in os.listdir("./test_vms_sensitivity_2d"):
                os.remove("./test_vms_sensitivity_2d/cylinder_test.dat")

    def test_SteadyCylinder(self):
        with ControlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            # solve fluid
            self.solve('test_vms_sensitivity_2d/steady_cylinder_test')
            # solve adjoint
            test = self.createTest('test_vms_sensitivity_2d/steady_cylinder_test_adjoint')
            test.Solve()
            Sensitivity = [[]]
            Sensitivity[0].append(test.main_model_part.GetNode(1968).GetSolutionStepValue(SHAPE_SENSITIVITY_X))
            Sensitivity[0].append(test.main_model_part.GetNode(1968).GetSolutionStepValue(SHAPE_SENSITIVITY_Y))

            # calculate sensitivity by finite difference
            step_size = 0.00000001
            FDSensitivity = self.ComputeFiniteDifferenceDragSensitivity([1968],step_size,'./test_vms_sensitivity_2d/steady_cylinder_test',[1.0,0.0,0.0],'./test_vms_sensitivity_2d/steady_cylinder_test.dat')
            self.assertAlmostEqual(Sensitivity[0][0], FDSensitivity[0][0], 4)
            self.assertAlmostEqual(Sensitivity[0][1], FDSensitivity[0][1], 2)
            # remove hdf5 file
            if "steady_cylinder_test_0.h5" in os.listdir("./test_vms_sensitivity_2d"):
                os.remove("./test_vms_sensitivity_2d/steady_cylinder_test_0.h5")
            # remove drag file
            if "steady_cylinder_test.dat" in os.listdir("./test_vms_sensitivity_2d"):
                os.remove("./test_vms_sensitivity_2d/steady_cylinder_test.dat")

    def tearDown(self):
        pass

if __name__ == '__main__':
    KratosUnittest.main()
