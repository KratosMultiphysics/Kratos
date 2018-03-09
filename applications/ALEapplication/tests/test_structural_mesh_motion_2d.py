import os
import KratosMultiphysics
import KratosMultiphysics.ALEApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils
import ale_analysis

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
            self.project_parameters = KratosMultiphysics.Parameters(parameter_file.read())

    def test_Rectangle_2D3N(self):
        with ControlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            self.createTest('test_structural_mesh_motion_2d/rectangle_2D3N_test')
            ale_analysis.ALEAnalysis(self.project_parameters).Run()
            kratos_utils.DeleteFileIfExisting("./test_mdpa_files/rectangle_2D3N_test.time")

    def test_Rectangle_2D4N(self):
        with ControlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            self.createTest('test_structural_mesh_motion_2d/rectangle_2D4N_test')
            ale_analysis.ALEAnalysis(self.project_parameters).Run()
            kratos_utils.DeleteFileIfExisting("./test_mdpa_files/rectangle_2D4N_test.time")


if __name__ == '__main__':
    KratosUnittest.main()
