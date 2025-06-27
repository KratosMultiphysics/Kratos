import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_helper

import os


class KratosGeoMechanicsHydraulicHeads(KratosUnittest.TestCase):
    def calculate_head(self, head_bottom, head_top, y):
        # Assuming a layer thickness of 1.0
        return head_bottom + (head_top - head_bottom) * y

    def run_and_check_nodal_hydraulic_heads(self, test_no, head_bottom, head_top):
        import KratosMultiphysics.GeoMechanicsApplication.run_geo_flow as run_geo_flow

        test_path = test_helper.get_file_path("test_nodal_hydraulic_head")
        status = run_geo_flow.run_flow_analysis(
            test_path, f"ProjectParameters_{test_no}.json", ""
        )

        self.assertEqual(status, 0)

        reader = test_helper.GiDOutputFileReader()
        output_data = reader.read_output_from(
            os.path.join(test_path, f"test_nodal_hydraulic_head_{test_no}.post.res")
        )
        time = 1.0
        numerical_hydraulic_heads = reader.nodal_values_at_time(
            "HYDRAULIC_HEAD", time, output_data
        )

        post_msh_file_path = os.path.join(
            test_path, f"test_nodal_hydraulic_head_{test_no}.post.msh"
        )
        nodal_coordinates = test_helper.read_coordinates_from_post_msh_file(
            post_msh_file_path
        )

        analytical_hydraulic_heads = [
            self.calculate_head(head_bottom, head_top, coord[1])
            for coord in nodal_coordinates
        ]

        self.assertVectorAlmostEqual(
            numerical_hydraulic_heads, analytical_hydraulic_heads, places=5
        )

    def test_hydraulic_heads_1(self):
        test_no = 1
        head_bottom = 1.0
        head_top = 1.0
        self.run_and_check_nodal_hydraulic_heads(test_no, head_bottom, head_top)

    def test_hydraulic_heads_2(self):
        test_no = 2
        head_bottom = 0.0
        head_top = 1.0
        self.run_and_check_nodal_hydraulic_heads(test_no, head_bottom, head_top)

    def test_hydraulic_heads_3(self):
        test_no = 3
        head_bottom = -1.0
        head_top = 1.0
        self.run_and_check_nodal_hydraulic_heads(test_no, head_bottom, head_top)

    def test_hydraulic_heads_4(self):
        test_no = 4
        head_bottom = 1.0
        head_top = 0.0
        self.run_and_check_nodal_hydraulic_heads(test_no, head_bottom, head_top)


if __name__ == "__main__":
    KratosUnittest.main()
