import os
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.GeoMechanicsApplication.gid_output_file_reader import GiDOutputFileReader
import test_helper


class KratosGeoMechanicsMovingColumnWithFixedPressureAbovePhreaticLine(
    KratosUnittest.TestCase
):
    """
    This test class checks a moving column with fixed pressure above the phreatic line
    using the C++ route. The column is moving downwards and upwards, crossing the phreatic line.
    The assertion is done that the water pressure at a specific node switches from a free value
    (below the phreatic line) to a fixed value (above the phreatic line) and vice versa when the "move_mesh" flag is set
    to `true`.
    The validation is done for both cases: with and without moving mesh (other than this option, the two
    cases are identical).
    """

    def test_fixed_pressure_above_phreatic_line_via_cpp_workflow_with_move_mesh(self):
        import KratosMultiphysics.GeoMechanicsApplication.run_geo_settlement as run_geo_settlement

        test_folder = os.path.join(
            test_helper.get_file_path(
                "moving_column_with_fixed_pressure_above_phreatic_line"
            ),
            "with_move_mesh",
        )
        status = run_geo_settlement.run_stages(test_folder, ["ProjectParameters.json"])
        self.assertEqual(status, 0)

        self.assert_results(test_folder, expected_value1 = -4479.17, expected_value2 = -4479.17)

    def test_fixed_pressure_above_phreatic_line_via_cpp_workflow_without_move_mesh(
        self,
    ):
        import KratosMultiphysics.GeoMechanicsApplication.run_geo_settlement as run_geo_settlement

        test_folder = os.path.join(
            test_helper.get_file_path(
                "moving_column_with_fixed_pressure_above_phreatic_line"
            ),
            "without_move_mesh",
        )
        status = run_geo_settlement.run_stages(test_folder, ["ProjectParameters.json"])
        self.assertEqual(status, 0)

        # Since "move_mesh" has been set to `false`, the nodal positions remain fixed. Consequently, the nodes with
        # fixed water pressures will not change, since we consider the undeformed mesh. Prescribing a displacement
        # field won't have any effect on that.
        self.assert_results(test_folder, expected_value1 = 0.0, expected_value2 = 0.0)

    def assert_results(self, test_folder, expected_value1, expected_value2):
        reader = GiDOutputFileReader()
        output_data = reader.read_output_from(
            os.path.join(test_folder, "output.post.res")
        )
        # The Cross-over time is at 21600.0 seconds (the time it takes for the column to move 0.5 m downwards),
        self.assertAlmostEqual(
            reader.nodal_values_at_time("WATER_PRESSURE", 25200.0, output_data, [27])[
                0
            ],
            0.0,
            5,
        )
        self.assertAlmostEqual(
            reader.nodal_values_at_time("WATER_PRESSURE", 28800.0, output_data, [27])[
                0
            ],
            expected_value1,
            5,
        )
        # When it's moving upwards again, we have another theoretical cross-over time at 151200 seconds. However,
        # since the DoF is fixed based on the actual nodal position (which is potentially updated at the end of the time
        # step depending on the "move_mesh" flag), the real cross-over happens one step later.
        self.assertAlmostEqual(
            reader.nodal_values_at_time("WATER_PRESSURE", 151200, output_data, [27])[0],
            expected_value2,
            5,
        )
        self.assertAlmostEqual(
            reader.nodal_values_at_time("WATER_PRESSURE", 154800, output_data, [27])[0],
            0.0,
            5,
        )


if __name__ == "__main__":
    KratosUnittest.main()
