import os
import KratosMultiphysics
import KratosMultiphysics.MeshMovingApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils
import mesh_moving_test_case

class TestCase(mesh_moving_test_case.MeshMovingTestCase):
    def test_Rectangle_2D3N(self):
        with mesh_moving_test_case.ControlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            self.createTest('test_laplacian_mesh_motion_2d/rectangle_2D3N_test')
            self.runTest()
            kratos_utils.DeleteFileIfExisting("./test_mdpa_files/rectangle_2D3N_test.time")

    def test_Rectangle_2D4N(self):
        with mesh_moving_test_case.ControlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            self.createTest('test_laplacian_mesh_motion_2d/rectangle_2D4N_test')
            self.runTest()
            kratos_utils.DeleteFileIfExisting("./test_mdpa_files/rectangle_2D4N_test.time")


if __name__ == '__main__':
    KratosUnittest.main()
