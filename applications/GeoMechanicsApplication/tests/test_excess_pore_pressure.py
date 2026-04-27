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
            stage2_output_data = self._run_two_stage_interface_case(case_path)
        else:
            stage2_output_data = self._run_two_stage_soil_case(case_path)

        stage1_output_data = self._read_output(case_path, "stage1")
        stage1_time_steps = self._all_times(stage1_output_data)
        stage2_time_steps = self._all_times(stage2_output_data)
        self.assertTrue(stage1_time_steps)
        self.assertTrue(stage2_time_steps)
        return stage1_output_data, stage2_output_data, stage1_time_steps, stage2_time_steps

    def _assert_has_excess_pore_pressure(
        self,
        stage1_output_data,
        stage2_output_data,
        stage1_time_steps,
        stage2_time_steps,
        threshold=1.0,
    ):
        stage1_reference_pressures = GiDOutputFileReader.nodal_values_at_time(
            "WATER_PRESSURE", stage1_time_steps[-1], stage1_output_data
        )

        max_abs_pressure_increment = 0.0
        for time_step in stage2_time_steps:
            nodal_pressures = GiDOutputFileReader.nodal_values_at_time(
                "WATER_PRESSURE", time_step, stage2_output_data
            )
            nodal_pressure_increments = [
                pressure - reference_pressure
                for pressure, reference_pressure in zip(
                    nodal_pressures, stage1_reference_pressures
                )
            ]
            if nodal_pressure_increments:
                max_abs_pressure_increment = max(
                    max_abs_pressure_increment,
                    max(abs(value) for value in nodal_pressure_increments),
                )

        self.assertGreater(
            max_abs_pressure_increment,
            threshold,
            msg=(
                "No meaningful excess pore pressure increment detected between stage 1 and stage 2. "
                f"Maximum absolute WATER_PRESSURE increment was {max_abs_pressure_increment:.12e} Pa"
            ),
        )

    @staticmethod
    def _max_abs_matrix_result_increment(
        result_name, stage1_output_data, stage2_output_data, stage1_time, stage2_time
    ):
        stage1_items = [
            item
            for item in stage1_output_data.get("results", {}).get(result_name, [])
            if item.get("time") == stage1_time
        ]
        stage2_items = [
            item
            for item in stage2_output_data.get("results", {}).get(result_name, [])
            if item.get("time") == stage2_time
        ]

        max_abs_increment = 0.0
        for stage1_item, stage2_item in zip(stage1_items, stage2_items):
            for stage1_element_values, stage2_element_values in zip(
                stage1_item.get("values", []), stage2_item.get("values", [])
            ):
                for stage1_gp_values, stage2_gp_values in zip(
                    stage1_element_values.get("value", []),
                    stage2_element_values.get("value", []),
                ):
                    for stage1_component, stage2_component in zip(
                        stage1_gp_values, stage2_gp_values
                    ):
                        max_abs_increment = max(
                            max_abs_increment, abs(stage2_component - stage1_component)
                        )

        return max_abs_increment

    def _assert_has_total_stress_increment(
        self,
        stage1_output_data,
        stage2_output_data,
        stage1_time_steps,
        stage2_time_steps,
        threshold=1.0,
    ):
        max_abs_total_stress_increment = self._max_abs_matrix_result_increment(
            "TOTAL_STRESS_TENSOR",
            stage1_output_data,
            stage2_output_data,
            stage1_time_steps[-1],
            stage2_time_steps[-1],
        )

        self.assertGreater(
            max_abs_total_stress_increment,
            threshold,
            msg=(
                "No meaningful TOTAL_STRESS_TENSOR increment detected between stage 1 and stage 2. "
                f"Maximum absolute TOTAL_STRESS_TENSOR increment was {max_abs_total_stress_increment:.12e}"
            ),
        )

    def test_same_order_column_generates_excess_pore_pressure(self):
        stage1_output_data, stage2_output_data, stage1_time_steps, stage2_time_steps = self._run_case("column")
        self._assert_has_excess_pore_pressure(
            stage1_output_data, stage2_output_data, stage1_time_steps, stage2_time_steps
        )
        self._assert_has_total_stress_increment(
            stage1_output_data, stage2_output_data, stage1_time_steps, stage2_time_steps
        )

    def test_diff_order_column_generates_excess_pore_pressure(self):
        stage1_output_data, stage2_output_data, stage1_time_steps, stage2_time_steps = self._run_case(
            "column_diff_order_elements"
        )
        self._assert_has_excess_pore_pressure(
            stage1_output_data, stage2_output_data, stage1_time_steps, stage2_time_steps
        )
        self._assert_has_total_stress_increment(
            stage1_output_data, stage2_output_data, stage1_time_steps, stage2_time_steps
        )

    def test_horizontal_interface_generates_excess_pore_pressure(self):
        stage1_output_data, stage2_output_data, stage1_time_steps, stage2_time_steps = self._run_case(
            "column_horizontal_interface", use_interface_parameters=True
        )
        self._assert_has_excess_pore_pressure(
            stage1_output_data, stage2_output_data, stage1_time_steps, stage2_time_steps
        )
        self._assert_has_total_stress_increment(
            stage1_output_data, stage2_output_data, stage1_time_steps, stage2_time_steps
        )

if __name__ == "__main__":
    KratosUnittest.main()
