import os
from KratosMultiphysics import *
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils
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

    def createTest(self, parameter_file_name):
        with open(parameter_file_name + '_parameters.json', 'r') as parameter_file:
            project_parameters = Parameters(parameter_file.read())
            parameter_file.close()
        test = test_MainKratos.MainKratos(project_parameters)
        return test

    def test_Rectangle_2D3N(self):
        with ControlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            test = self.createTest('test_structural_mesh_motion_2d/mpi_rectangle_2D3N_test')
            test.Solve()
            kratos_utils.DeleteFileIfExisting("./test_mdpa_files/rectangle_2D3N_test.time")

    def test_Rectangle_2D4N(self):
        with ControlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            test = self.createTest('test_structural_mesh_motion_2d/mpi_rectangle_2D4N_test')
            test.Solve()
            kratos_utils.DeleteFileIfExisting("./test_mdpa_files/rectangle_2D4N_test.time")


if __name__ == '__main__':
    KratosUnittest.main()
