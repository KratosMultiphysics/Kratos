import os

import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_helper
import KratosMultiphysics.GeoMechanicsApplication.run_multiple_stages as run_multiple_stages
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

    def _run_soil_case(self, case_name):
        case_path = test_helper.get_file_path(os.path.join(self._ROOT, case_name))
        output_data = self._run_two_stage_soil_case(case_path)
        return case_path, output_data, self._all_times(output_data)

    def _run_soil_case_with_input_folder(self, case_name, input_folder_name):
        case_path = test_helper.get_file_path(os.path.join(self._ROOT, case_name))
        input_path = test_helper.get_file_path(os.path.join(self._ROOT, input_folder_name))
        run_multiple_stages.run_stages(
            case_path,
            2,
            "ProjectParameters_stage{}.json",
            input_path=input_path,
        )
        output_data = self._read_output(case_path, "stage2")
        return case_path, output_data, self._all_times(output_data)

    def _assert_matrix_result_match_at_time(
        self,
        result_name,
        first_output,
        second_output,
        first_time,
        second_time,
        time_index,
        tolerance,
    ):
        first_items = self._result_items_at_time(result_name, first_time, first_output)
        second_items = self._result_items_at_time(result_name, second_time, second_output)

        self.assertEqual(
            len(first_items),
            len(second_items),
            msg=(
                f"Different number of {result_name} result blocks at time index {time_index}: "
                f"first={len(first_items)}, second={len(second_items)}"
            ),
        )

        max_abs_difference = 0.0
        for first_item, second_item in zip(first_items, second_items):
            first_values_by_element = self._element_id_to_integration_values_map(first_item)
            second_values_by_element = self._element_id_to_integration_values_map(second_item)

            self.assertEqual(
                sorted(first_values_by_element.keys()),
                sorted(second_values_by_element.keys()),
                msg=(
                    f"Different {result_name} element IDs at time index {time_index} "
                    f"(time={first_time:.12g})"
                ),
            )

            for element_id in sorted(first_values_by_element.keys()):
                first_integration_values = first_values_by_element[element_id]
                second_integration_values = second_values_by_element[element_id]

                self.assertEqual(
                    len(first_integration_values),
                    len(second_integration_values),
                    msg=(
                        f"Different number of {result_name} integration points for element {element_id} "
                        f"at time index {time_index}"
                    ),
                )

                for integration_point_index, (
                    first_vector,
                    second_vector,
                ) in enumerate(zip(first_integration_values, second_integration_values)):
                    self.assertEqual(len(first_vector), len(second_vector))
                    for component_index, (first_component, second_component) in enumerate(
                        zip(first_vector, second_vector)
                    ):
                        max_abs_difference = max(
                            max_abs_difference, abs(first_component - second_component)
                        )

        self.assertLess(
            max_abs_difference,
            tolerance,
            msg=(
                f"{result_name} mismatch at time index {time_index} "
                f"(time={first_time:.12g}). Maximum absolute component difference was "
                f"{max_abs_difference:.12e}"
            ),
        )

    def _assert_case_matches_reference_solution(
        self,
        reference_case_name,
        comparison_case_name,
        comparison_input_folder_name,
        min_shared_unique_pairs,
        displacement_tolerance,
        pressure_tolerance,
        flux_tolerance,
        total_stress_tolerance,
    ):
        reference_case_path, reference_output, reference_times = self._run_soil_case(reference_case_name)
        comparison_case_path, comparison_output, comparison_times = self._run_soil_case_with_input_folder(
            comparison_case_name, comparison_input_folder_name
        )

        reference_node_ids, comparison_node_ids = self._shared_unique_node_ids(
            reference_case_path,
            comparison_case_path,
            min_pairs=min_shared_unique_pairs,
        )

        self._assert_nodal_match(
            reference_output,
            comparison_output,
            reference_times,
            comparison_times,
            reference_node_ids,
            comparison_node_ids,
            displacement_tolerance=displacement_tolerance,
            pressure_tolerance=pressure_tolerance,
            flux_tolerance=flux_tolerance,
        )

        for time_index, (reference_time, comparison_time) in enumerate(
            zip(reference_times, comparison_times)
        ):
            self._assert_matrix_result_match_at_time(
                "TOTAL_STRESS_TENSOR",
                reference_output,
                comparison_output,
                reference_time,
                comparison_time,
                time_index,
                total_stress_tolerance,
            )

    def _assert_soil_case_match(
        self,
        first_case_name,
        second_case_name,
        min_shared_unique_pairs,
        displacement_tolerance,
        pressure_tolerance,
    ):
        first_case_path, first_output, first_times = self._run_soil_case(first_case_name)
        second_case_path, second_output, second_times = self._run_soil_case(second_case_name)

        first_node_ids, second_node_ids = self._shared_unique_node_ids(
            first_case_path,
            second_case_path,
            min_pairs=min_shared_unique_pairs,
        )

        self.assertEqual(
            len(first_times),
            len(second_times),
            msg=(
                "Different number of output time steps found: "
                f"first={len(first_times)}, second={len(second_times)}"
            ),
        )

        for time_index, (first_time, second_time) in enumerate(
            zip(first_times, second_times)
        ):
            self.assertTrue(
                test_upw_interface.math.isclose(first_time, second_time),
                msg=(
                    f"Time step mismatch at index {time_index}: "
                    f"first time={first_time:.12g}, second time={second_time:.12g}"
                ),
            )

            self._assert_displacement_match_at_time(
                first_output,
                second_output,
                first_time,
                second_time,
                first_node_ids,
                second_node_ids,
                time_index,
                displacement_tolerance,
            )
            self._assert_pressure_match_at_time(
                first_output,
                second_output,
                first_time,
                second_time,
                first_node_ids,
                second_node_ids,
                time_index,
                pressure_tolerance,
            )

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

    def test_same_order_and_diff_order_columns_match_on_shared_unique_nodes(self):
        self._assert_soil_case_match(
            first_case_name="column",
            second_case_name="column_diff_order_elements",
            min_shared_unique_pairs=3,
            displacement_tolerance=1e-6,
            pressure_tolerance=0.11,
        )

    def test_newton_raphson_column_matches_original_solution(self):
        self._assert_case_matches_reference_solution(
            reference_case_name="column",
            comparison_case_name="column_newton_raphson",
            comparison_input_folder_name="common_newton_raphson",
            min_shared_unique_pairs=6,
            displacement_tolerance=1e-6,
            pressure_tolerance=0.11,
            flux_tolerance=5e-6,
            total_stress_tolerance=1e-5,
        )

    def test_newton_raphson_diff_order_column_matches_original_solution(self):
        self._assert_case_matches_reference_solution(
            reference_case_name="column_diff_order_elements",
            comparison_case_name="column_diff_order_elements_newton_raphson",
            comparison_input_folder_name="common_newton_raphson",
            min_shared_unique_pairs=15,
            displacement_tolerance=1e-6,
            pressure_tolerance=0.11,
            flux_tolerance=5e-6,
            total_stress_tolerance=1e-5,
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
