import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.GeoMechanicsApplication.gid_output_file_reader import GiDOutputFileReader
import test_helper

import os
import parameterized


class KratosGeoMechanicsHydraulicHeads(KratosUnittest.TestCase):
    def setUp(self):
        super().setUp()

        self.test_path = test_helper.get_file_path("test_nodal_hydraulic_head")

    @staticmethod
    def calculate_head(head_bottom, head_top, y):
        # Assuming a layer thickness of 1.0
        return head_bottom + (head_top - head_bottom) * y

    @parameterized.parameterized.expand(
        [
            (1, 1.0, 1.0),
            (2, 0.0, 1.0),
            (3, -1.0, 1.0),
            (4, 1.0, 0.0),
        ]
    )
    def test_nodal_hydraulic_heads(self, test_no, head_bottom, head_top):
        import KratosMultiphysics.GeoMechanicsApplication.run_geo_flow as run_geo_flow

        status = run_geo_flow.run_flow_analysis(
            self.test_path, f"ProjectParameters_{test_no}.json", ""
        )

        self.assertEqual(status, 0)

        reader = GiDOutputFileReader()
        output_data = reader.read_output_from(
            os.path.join(
                self.test_path, f"test_nodal_hydraulic_head_{test_no}.post.res"
            )
        )
        time = 1.0
        numerical_hydraulic_heads = reader.nodal_values_at_time(
            "HYDRAULIC_HEAD", time, output_data
        )

        post_msh_file_path = os.path.join(
            self.test_path, f"test_nodal_hydraulic_head_{test_no}.post.msh"
        )
        nodal_coordinates = test_helper.read_coordinates_from_post_msh_file(
            post_msh_file_path
        )
        # Make sure that we check all nodes
        self.assertEqual(len(nodal_coordinates), 251)

        analytical_hydraulic_heads = [
            self.calculate_head(head_bottom, head_top, coord[1])
            for coord in nodal_coordinates
        ]

        self.assertVectorAlmostEqual(
            numerical_hydraulic_heads, analytical_hydraulic_heads, places=5
        )


if __name__ == "__main__":
    KratosUnittest.main()
