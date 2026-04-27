import os

import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_helper
from KratosMultiphysics.GeoMechanicsApplication.gid_output_file_reader import (
    GiDOutputFileReader,
)
import test_upw_interface


class KratosGeoMechanicsExcessPorePressureTests(test_upw_interface.KratosGeoMechanicsUPwInterfaceTests):
    """
    Integration tests for excess pore pressure scenarios based on UPw interface setup.
    """

    # Suppress inherited test methods from the parent class — only run the tests defined here.
    test_vertical_interface_matches_column_on_shared_unique_nodes = None
    test_horizontal_interface_matches_column_on_shared_nodes = None
    test_horizontal_interface_diff_order_matches_column_diff_order = None

    _ROOT = "excess_pore_pressure"

    def _run_case(self, case_name, use_interface_parameters=False):
        case_path = test_helper.get_file_path(os.path.join(self._ROOT, case_name))
        if use_interface_parameters:
            output_data = self._run_two_stage_interface_case(case_path)
        else:
            output_data = self._run_two_stage_soil_case(case_path)

        time_steps = self._all_times(output_data)
        self.assertTrue(time_steps)
        return output_data, time_steps

    def _assert_has_excess_pore_pressure(self, output_data, time_steps, threshold=1.0):
        max_abs_pressure = 0.0
        for time_step in time_steps:
            nodal_pressures = GiDOutputFileReader.nodal_values_at_time(
                "WATER_PRESSURE", time_step, output_data
            )
            if nodal_pressures:
                max_abs_pressure = max(
                    max_abs_pressure, max(abs(value) for value in nodal_pressures)
                )

        self.assertGreater(
            max_abs_pressure,
            threshold,
            msg=(
                "No meaningful excess pore pressure detected. "
                f"Maximum absolute WATER_PRESSURE was {max_abs_pressure:.12e} Pa"
            ),
        )

    def _sum_vertical_reaction_on_prescribed_top_displacement(
        self, output_data, time_step, prescribed_top_displacement=-0.01, tolerance=1.0e-10
    ):
        displacements = GiDOutputFileReader.nodal_values_at_time(
            "TOTAL_DISPLACEMENT", time_step, output_data
        )
        reactions = GiDOutputFileReader.nodal_values_at_time(
            "REACTION", time_step, output_data
        )

        top_node_indices = [
            index
            for index, displacement in enumerate(displacements)
            if abs(displacement[1] - prescribed_top_displacement) < tolerance
        ]
        self.assertTrue(
            top_node_indices,
            msg="No nodes found with prescribed top displacement for reaction check",
        )

        return sum(reactions[index][1] for index in top_node_indices)

    def test_same_order_column_generates_excess_pore_pressure(self):
        output_data, time_steps = self._run_case("column")
        self._assert_has_excess_pore_pressure(output_data, time_steps)

    def test_diff_order_column_generates_excess_pore_pressure(self):
        output_data, time_steps = self._run_case("column_diff_order_elements")
        self._assert_has_excess_pore_pressure(output_data, time_steps)

    def test_horizontal_interface_generates_excess_pore_pressure(self):
        output_data, time_steps = self._run_case(
            "column_horizontal_interface", use_interface_parameters=True
        )
        self._assert_has_excess_pore_pressure(output_data, time_steps)

    def test_diff_order_column_reaction_matches_same_order_under_prescribed_displacement(self):
        same_order_output, same_order_times = self._run_case("column")
        diff_order_output, diff_order_times = self._run_case("column_diff_order_elements")

        same_order_reaction_y = self._sum_vertical_reaction_on_prescribed_top_displacement(
            same_order_output, same_order_times[-1]
        )
        diff_order_reaction_y = self._sum_vertical_reaction_on_prescribed_top_displacement(
            diff_order_output, diff_order_times[-1]
        )

        self.assertAlmostEqual(
            diff_order_reaction_y,
            same_order_reaction_y,
            delta=1.0e-10 * max(1.0, abs(diff_order_reaction_y), abs(same_order_reaction_y)),
            msg=(
                "Under prescribed top displacement, diff-order and same-order should produce matching global reaction. "
                f"same-order={same_order_reaction_y:.12e}, "
                f"diff-order={diff_order_reaction_y:.12e}, "
                f"|delta|={abs(diff_order_reaction_y - same_order_reaction_y):.12e}"
            ),
        )


if __name__ == "__main__":
    KratosUnittest.main()
