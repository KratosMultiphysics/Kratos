import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_helper

import os


class KratosGeoMechanicsHydraulicHeads(KratosUnittest.TestCase):
    def calculate_head(self, head_bottom, head_top, y):
        # Assuming a layer thickness of 1.0
        return head_bottom + (head_top - head_bottom) * y

    def test_hydraulic_heads(self):
        import KratosMultiphysics.GeoMechanicsApplication.run_geo_flow as run_geo_flow

        test_path = test_helper.get_file_path(
            "test_head_extrapolation_custom_workflow_flow"
        )
        status = run_geo_flow.run_flow_analysis(
            test_path, "ProjectParameters_1.json", ""
        )

        self.assertEqual(status, 0)

        reader = test_helper.GiDOutputFileReader()
        output_data = reader.read_output_from(
            os.path.join(test_path, "test_head_extrapolate_1.post.res")
        )
        time = 1.0
        numerical_hydraulic_heads = reader.nodal_values_at_time(
            "HYDRAULIC_HEAD", time, output_data
        )

        post_msh_file_path = os.path.join(test_path, "test_head_extrapolate_1.post.msh")
        nodal_coordinates = test_helper.read_coordinates_from_post_msh_file(
            post_msh_file_path
        )

        head_bottom = 1.0
        head_top = 1.0
        analytical_hydraulic_heads = [
            self.calculate_head(head_bottom, head_top, coord[1])
            for coord in nodal_coordinates
        ]

        self.assertVectorAlmostEqual(
            numerical_hydraulic_heads, analytical_hydraulic_heads
        )


if __name__ == "__main__":
    KratosUnittest.main()
