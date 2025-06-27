import os
import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_helper


class KratosGeoMechanicsMovingColumnWithFixedPressureAbovePhreaticLine(KratosUnittest.TestCase):
    """
    This test class checks the settlement workflow using the C++ route.
    """


    def test_fixed_pressure_above_phreatic_line_via_cpp_workflow(self):
        import KratosMultiphysics.GeoMechanicsApplication.run_geo_settlement as run_geo_settlement
        test_root = test_helper.get_file_path(
            "moving_column_with_fixed_pressure_above_phreatic_line"
        )
        status = run_geo_settlement.run_stages(
            test_root, ["ProjectParameters.json"]
        )
        reader = test_helper.GiDOutputFileReader()
        output_data = reader.read_output_from(
            os.path.join(test_root, "output.post.res")
        )

        self.assertEqual(status, 0)
        self.assertAlmostEqual(
            reader.nodal_values_at_time("WATER_PRESSURE", 21600.0, output_data, [27])[
                0
            ],
            0.0,
            5,
        )
        self.assertAlmostEqual(
            reader.nodal_values_at_time("WATER_PRESSURE", 22200.0, output_data, [27])[
                0
            ],
            -3678.75,
            5,
        )


if __name__ == "__main__":
    KratosUnittest.main()
