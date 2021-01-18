import KratosMultiphysics as KM
from  mesh_moving_test_case import MeshMovingTestCase

class TestCase(MeshMovingTestCase):

    def test_laplacian_2D3N(self):
        # General Settings for the test
        self.domain_size = 2
        self.number_of_nodes_per_elements = 3
        self.solver_type = "laplacian"

        mesh_vel_calc_settings = KM.Parameters("""{
            "calculate_mesh_velocity"   : true,
            "mesh_velocity_calculation" : {
                "time_scheme" : "bdf1"
            }
        }""")

        self.print_reference_results = False

        # to suppress many prints
        self.print_logger_info = False

        # Set to true to get post-process files for the test
        self.print_vtk_output = False
        self.print_gid_output = False

        self.executeTest(mesh_vel_calc_settings)

    def test_laplacian_2D4N(self):
        # General Settings for the test
        self.domain_size = 2
        self.number_of_nodes_per_elements = 4
        self.solver_type = "laplacian"

        mesh_vel_calc_settings = KM.Parameters("""{
            "calculate_mesh_velocity"   : true,
            "mesh_velocity_calculation" : {
                "time_scheme" : "bdf2"
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
