import KratosMultiphysics as KM
from  mesh_moving_test_case import MeshMovingTestCase

class TestCase(MeshMovingTestCase):

    def test_laplacian_3D4N(self):
        # General Settings for the test
        self.domain_size = 3
        self.number_of_nodes_per_elements = 4
        self.solver_type = "laplacian"

        mesh_vel_calc_settings = KM.Parameters("""{
            "calculate_mesh_velocity"   : true,
            "mesh_velocity_calculation" : {
                "time_scheme" : "newmark"
            }
        }""")

        self.print_reference_results = False

        # to suppress many prints
        self.print_logger_info = False

        # Set to true to get post-process files for the test
        self.print_vtk_output = False
        self.print_gid_output = False

        self.executeTest(mesh_vel_calc_settings)

    def test_laplacian_3D8N(self):
        # General Settings for the test
        self.domain_size = 3
        self.number_of_nodes_per_elements = 8
        self.solver_type = "laplacian"

        mesh_vel_calc_settings = KM.Parameters("""{
            "calculate_mesh_velocity"   : true,
            "mesh_velocity_calculation" : {
                "time_scheme" : "bossak",
                "alpha_m" : -0.3
            }
        }""")

        self.print_reference_results = False

        # to suppress many prints
        self.print_logger_info = False

        # Set to true to get post-process files for the test
        self.print_vtk_output = False
        self.print_gid_output = False

        self.executeTest(mesh_vel_calc_settings)

if __name__ == '__main__':
    import KratosMultiphysics.KratosUnittest as KratosUnittest
    KratosUnittest.main()
