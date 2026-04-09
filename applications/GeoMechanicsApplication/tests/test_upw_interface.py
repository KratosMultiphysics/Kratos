import math
import os

import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_helper
import KratosMultiphysics.GeoMechanicsApplication.run_multiple_stages as run_multiple_stages
from KratosMultiphysics.GeoMechanicsApplication.gid_output_file_reader import (
    GiDOutputFileReader,
)


class KratosGeoMechanicsUPwInterfaceTests(KratosUnittest.TestCase):
    """
    This class contains benchmark tests which check if interfaces between soil and structural elements are
    calculated correctly
    """

    def _run_two_stage_case(self, case_path, project_parameters_filename_pattern):
        run_multiple_stages.run_stages(
            case_path, 2, project_parameters_filename_pattern, input_path="../common"
        )
        return self._read_output(case_path, "stage2")

    def _run_two_stage_interface_case(self, case_path):
        return self._run_two_stage_case(
            case_path, "ProjectParameters_interface_stage{}.json"
        )

    def _run_two_stage_soil_case(self, case_path):
        return self._run_two_stage_case(case_path, "ProjectParameters_stage{}.json")

    def _read_output(self, case_path, stage_name):
        output_file_name = os.path.join(case_path, f"{stage_name}.post.res")
        return GiDOutputFileReader().read_output_from(output_file_name)

    @staticmethod
    def _all_times(output_data):
        times = GiDOutputFileReader.get_time_steps_from_first_valid_result(output_data)
        if not times:
            raise RuntimeError("No valid output time steps found")
        return times

    @staticmethod
    def _top_vector_component_differences(
        vectors_a,
        vectors_b,
        ids_a,
        ids_b,
        component_index,
        top_n=5,
    ):
        differences = [
            (
                abs(vector_a[component_index] - vector_b[component_index]),
                id_a,
                id_b,
                vector_a[component_index],
                vector_b[component_index],
            )
            for id_a, id_b, vector_a, vector_b in zip(
                ids_a, ids_b, vectors_a, vectors_b
            )
        ]
        return sorted(differences, key=lambda item: item[0], reverse=True)[:top_n]

    @staticmethod
    def _top_scalar_differences(values_a, values_b, ids_a, ids_b, top_n=5):
        differences = [
            (
                abs(value_a - value_b),
                id_a,
                id_b,
                value_a,
                value_b,
            )
            for id_a, id_b, value_a, value_b in zip(ids_a, ids_b, values_a, values_b)
        ]
        return sorted(differences, key=lambda item: item[0], reverse=True)[:top_n]

    @staticmethod
    def _format_nodal_difference_summary(differences):
        if not differences:
            return "No nodal pairs were compared"

        return " | ".join(
            (
                f"interface node {interface_node_id}, soil node {soil_node_id}: "
                f"interface={interface_value:.12e}, soil={soil_value:.12e}, "
                f"|delta|={difference:.12e}"
            )
            for (
                difference,
                interface_node_id,
                soil_node_id,
                interface_value,
                soil_value,
            ) in differences
        )

    @staticmethod
    def _format_flux_difference_summary(differences):
        if not differences:
            return "No overlapping flux components were compared"

        return " | ".join(
            (
                f"element {element_id}, gp {integration_point_index}, "
                f"component {component_index}: "
                f"interface={interface_value:.12e}, soil={soil_value:.12e}, "
                f"|delta|={difference:.12e}"
            )
            for (
                difference,
                element_id,
                integration_point_index,
                component_index,
                interface_value,
                soil_value,
            ) in differences
        )

    @staticmethod
    def _result_items_at_time(result_item_name, time, output_data):
        return [
            item
            for item in output_data.get("results", {}).get(result_item_name, [])
            if math.isclose(item["time"], time)
        ]

    @staticmethod
    def _element_id_to_integration_values_map(result_item):
        return {
            value["element"]: value["value"] for value in result_item.get("values", [])
        }

    @staticmethod
    def _read_mdpa_nodes(mdpa_file_path):
        node_id_to_coordinates = {}
        reading_nodes = False

        with open(mdpa_file_path, "r", encoding="utf-8") as mdpa_file:
            for line in mdpa_file:
                stripped_line = line.strip()
                if stripped_line == "Begin Nodes":
                    reading_nodes = True
                    continue
                if stripped_line == "End Nodes":
                    break

                if reading_nodes and stripped_line:
                    columns = stripped_line.split()
                    if len(columns) >= 4:
                        node_id_to_coordinates[int(columns[0])] = (
                            float(columns[1]),
                            float(columns[2]),
                            float(columns[3]),
                        )

        return node_id_to_coordinates

    @staticmethod
    def _coordinate_to_node_ids_map(node_id_to_coordinates):
        coordinate_to_node_ids_map = {}
        for node_id, coordinates in node_id_to_coordinates.items():
            coordinate_key = tuple(round(value, 12) for value in coordinates)
            coordinate_to_node_ids_map.setdefault(coordinate_key, []).append(node_id)
        return coordinate_to_node_ids_map

    @staticmethod
    def _shared_unique_node_id_pairs(first_mdpa_file_path, second_mdpa_file_path):
        first_nodes = KratosGeoMechanicsUPwInterfaceTests._read_mdpa_nodes(
            first_mdpa_file_path
        )
        second_nodes = KratosGeoMechanicsUPwInterfaceTests._read_mdpa_nodes(
            second_mdpa_file_path
        )

        first_coordinate_to_node_ids_map = (
            KratosGeoMechanicsUPwInterfaceTests._coordinate_to_node_ids_map(first_nodes)
        )
        second_coordinate_to_node_ids_map = (
            KratosGeoMechanicsUPwInterfaceTests._coordinate_to_node_ids_map(
                second_nodes
            )
        )

        shared_coordinates = sorted(
            set(first_coordinate_to_node_ids_map.keys()).intersection(
                second_coordinate_to_node_ids_map.keys()
            )
        )

        node_id_pairs = []
        for coordinates in shared_coordinates:
            first_node_ids = first_coordinate_to_node_ids_map[coordinates]
            second_node_ids = second_coordinate_to_node_ids_map[coordinates]

            if len(first_node_ids) == 1 and len(second_node_ids) == 1:
                node_id_pairs.append((first_node_ids[0], second_node_ids[0]))

        return node_id_pairs

    @staticmethod
    def _interface_node_pairs_from_element_line(element_line):
        columns = element_line.split()
        if len(columns) < 4:
            return ()

        element_node_ids = [int(value) for value in columns[2:]]
        if len(element_node_ids) < 2 or len(element_node_ids) % 2 != 0:
            return ()

        half = int(len(element_node_ids) / 2)
        return zip(element_node_ids[:half], element_node_ids[half:])

    @staticmethod
    def _interface_unique_node_id_pairs(mdpa_file_path):
        node_id_pairs = []
        reading_interface_elements = False

        with open(mdpa_file_path, "r", encoding="utf-8") as mdpa_file:
            for stripped_line in map(str.strip, mdpa_file):

                if stripped_line.startswith("Begin Elements"):
                    element_name = stripped_line[len("Begin Elements") :].strip()
                    reading_interface_elements = "Interface" in element_name
                    continue

                if stripped_line == "End Elements":
                    reading_interface_elements = False
                    continue

                if not reading_interface_elements:
                    continue

                for (
                    pair
                ) in KratosGeoMechanicsUPwInterfaceTests._interface_node_pairs_from_element_line(
                    stripped_line
                ):
                    if pair not in node_id_pairs:
                        node_id_pairs.append(pair)

        return node_id_pairs

    def _shared_unique_node_ids(self, interface_path, soil_path, min_pairs):
        node_id_pairs = self._shared_unique_node_id_pairs(
            os.path.join(interface_path, "column.mdpa"),
            os.path.join(soil_path, "column.mdpa"),
        )
        self.assertGreaterEqual(len(node_id_pairs), min_pairs)

        interface_node_ids, soil_node_ids = zip(*node_id_pairs)

        return interface_node_ids, soil_node_ids

    def _all_interface_to_soil_node_ids(self, interface_path, soil_path, min_pairs):
        interface_mdpa_file_path = os.path.join(interface_path, "column.mdpa")
        soil_mdpa_file_path = os.path.join(soil_path, "column.mdpa")

        interface_node_id_pairs = self._interface_unique_node_id_pairs(
            interface_mdpa_file_path
        )
        self.assertGreaterEqual(len(interface_node_id_pairs), min_pairs)

        interface_node_ids = list(
            dict.fromkeys(
                node_id
                for interface_node_id_pair in interface_node_id_pairs
                for node_id in interface_node_id_pair
            )
        )

        nodes_of_interface_model = self._read_mdpa_nodes(interface_mdpa_file_path)
        nodes_of_soil_model = self._read_mdpa_nodes(soil_mdpa_file_path)
        soil_coordinate_to_node_ids_map = self._coordinate_to_node_ids_map(
            nodes_of_soil_model
        )

        soil_node_ids = []
        for interface_node_id in interface_node_ids:
            self.assertIn(interface_node_id, nodes_of_interface_model)

            coordinate_key = tuple(
                round(value, 12)
                for value in nodes_of_interface_model[interface_node_id]
            )
            matching_soil_node_ids = soil_coordinate_to_node_ids_map.get(
                coordinate_key, []
            )
            self.assertEqual(
                len(matching_soil_node_ids),
                1,
                msg=(
                    "Expected exactly one matching soil node for interface node "
                    f"{interface_node_id} at coordinates {coordinate_key}"
                ),
            )
            soil_node_ids.append(matching_soil_node_ids[0])

        return interface_node_ids, soil_node_ids

    def _interface_node_ids_by_side(self, interface_path, min_pairs):
        node_id_pairs = self._interface_unique_node_id_pairs(
            os.path.join(interface_path, "column.mdpa")
        )
        self.assertGreaterEqual(len(node_id_pairs), min_pairs)

        side_a_node_ids = [pair[0] for pair in node_id_pairs]
        side_b_node_ids = [pair[1] for pair in node_id_pairs]

        return side_a_node_ids, side_b_node_ids

    def _assert_displacement_match_at_time(
        self,
        interface_output,
        soil_output,
        interface_time,
        soil_time,
        interface_node_ids,
        soil_node_ids,
        time_index,
        displacement_tolerance,
    ):
        interface_displacements = GiDOutputFileReader.nodal_values_at_time(
            "TOTAL_DISPLACEMENT",
            interface_time,
            interface_output,
            interface_node_ids,
        )
        soil_displacements = GiDOutputFileReader.nodal_values_at_time(
            "TOTAL_DISPLACEMENT", soil_time, soil_output, soil_node_ids
        )
        top_displacement_differences = self._top_vector_component_differences(
            interface_displacements,
            soil_displacements,
            interface_node_ids,
            soil_node_ids,
            component_index=1,
        )
        max_displacement_difference_y = (
            top_displacement_differences[0][0] if top_displacement_differences else 0.0
        )
        self.assertLess(
            max_displacement_difference_y,
            displacement_tolerance,
            msg=(
                f"TOTAL_DISPLACEMENT_Y mismatch at time index {time_index} "
                f"(time={interface_time:.12g}). Top nodal differences: "
                f"{self._format_nodal_difference_summary(top_displacement_differences)}"
            ),
        )

    def _assert_pressure_match_at_time(
        self,
        interface_output,
        soil_output,
        interface_time,
        soil_time,
        interface_node_ids,
        soil_node_ids,
        time_index,
        pressure_tolerance,
    ):
        interface_pressures = GiDOutputFileReader.nodal_values_at_time(
            "WATER_PRESSURE", interface_time, interface_output, interface_node_ids
        )
        soil_pressures = GiDOutputFileReader.nodal_values_at_time(
            "WATER_PRESSURE", soil_time, soil_output, soil_node_ids
        )
        top_pressure_differences = self._top_scalar_differences(
            interface_pressures,
            soil_pressures,
            interface_node_ids,
            soil_node_ids,
        )
        max_pressure_difference = (
            top_pressure_differences[0][0] if top_pressure_differences else 0.0
        )
        self.assertLess(
            max_pressure_difference,
            pressure_tolerance,
            msg=(
                f"WATER_PRESSURE mismatch at time index {time_index} "
                f"(time={interface_time:.12g}). Top nodal differences: "
                f"{self._format_nodal_difference_summary(top_pressure_differences)}"
            ),
        )

    def _assert_flux_match_at_time(
        self,
        interface_output,
        soil_output,
        interface_time,
        soil_time,
        time_index,
        flux_tolerance,
    ):
        interface_flux_items = self._result_items_at_time(
            "FLUID_FLUX_VECTOR", interface_time, interface_output
        )
        soil_flux_items = self._result_items_at_time(
            "FLUID_FLUX_VECTOR", soil_time, soil_output
        )

        self.assertTrue(
            interface_flux_items,
            msg=(
                f"No interface FLUID_FLUX_VECTOR values found at "
                f"time index {time_index} (time={interface_time:.12g})"
            ),
        )
        self.assertTrue(
            soil_flux_items,
            msg=(
                f"No soil FLUID_FLUX_VECTOR values found at "
                f"time index {time_index} (time={soil_time:.12g})"
            ),
        )

        # The interface output can contain both line-interface and bulk-element
        # FLUID_FLUX_VECTOR blocks at the same time. Select the pair sharing
        # the largest set of element IDs (the bulk domain overlap).
        soil_flux_item = max(soil_flux_items, key=lambda item: len(item["values"]))
        soil_flux_by_element = self._element_id_to_integration_values_map(
            soil_flux_item
        )
        soil_element_ids = set(soil_flux_by_element.keys())

        def _overlap_with_soil(item):
            interface_element_ids = {value["element"] for value in item["values"]}
            return len(interface_element_ids.intersection(soil_element_ids))

        interface_flux_item = max(
            interface_flux_items,
            key=lambda item, overlap_with_soil=_overlap_with_soil: (
                overlap_with_soil(item),
                len(item["values"]),
            ),
        )
        interface_flux_by_element = self._element_id_to_integration_values_map(
            interface_flux_item
        )

        common_element_ids = sorted(
            set(interface_flux_by_element.keys()).intersection(soil_element_ids)
        )
        self.assertTrue(
            common_element_ids,
            msg=(
                f"No overlapping flux elements found at time index {time_index} "
                f"(time={interface_time:.12g})"
            ),
        )

        flux_component_differences = []
        for element_id in common_element_ids:
            interface_integration_values = interface_flux_by_element[element_id]
            soil_integration_values = soil_flux_by_element[element_id]

            self.assertEqual(
                len(interface_integration_values), len(soil_integration_values)
            )

            for integration_point_index, (
                interface_vector,
                soil_vector,
            ) in enumerate(zip(interface_integration_values, soil_integration_values)):
                self.assertEqual(len(interface_vector), len(soil_vector))
                for component_index, (
                    interface_component,
                    soil_component,
                ) in enumerate(zip(interface_vector, soil_vector)):
                    flux_component_differences.append(
                        (
                            abs(interface_component - soil_component),
                            element_id,
                            integration_point_index,
                            component_index,
                            interface_component,
                            soil_component,
                        )
                    )

        top_flux_differences = sorted(
            flux_component_differences,
            key=lambda item: item[0],
            reverse=True,
        )[:5]
        max_flux_component_difference = (
            top_flux_differences[0][0] if top_flux_differences else 0.0
        )
        self.assertLess(
            max_flux_component_difference,
            flux_tolerance,
            msg=(
                f"FLUID_FLUX_VECTOR mismatch at time index {time_index} "
                f"(time={interface_time:.12g}). Top component differences: "
                f"{self._format_flux_difference_summary(top_flux_differences)}"
            ),
        )

    def _assert_nodal_match(
        self,
        interface_output,
        soil_output,
        interface_times,
        soil_times,
        interface_node_ids,
        soil_node_ids,
        displacement_tolerance,
        pressure_tolerance,
        flux_tolerance,
    ):
        self.assertEqual(
            len(interface_times),
            len(soil_times),
            msg=(
                "Different number of output time steps found: "
                f"interface={len(interface_times)}, soil={len(soil_times)}"
            ),
        )

        for time_index, (interface_time, soil_time) in enumerate(
            zip(interface_times, soil_times)
        ):
            self.assertTrue(
                math.isclose(interface_time, soil_time),
                msg=(
                    f"Time step mismatch at index {time_index}: "
                    f"interface time={interface_time:.12g}, soil time={soil_time:.12g}"
                ),
            )

            self._assert_displacement_match_at_time(
                interface_output,
                soil_output,
                interface_time,
                soil_time,
                interface_node_ids,
                soil_node_ids,
                time_index,
                displacement_tolerance,
            )
            self._assert_pressure_match_at_time(
                interface_output,
                soil_output,
                interface_time,
                soil_time,
                interface_node_ids,
                soil_node_ids,
                time_index,
                pressure_tolerance,
            )
            self._assert_flux_match_at_time(
                interface_output,
                soil_output,
                interface_time,
                soil_time,
                time_index,
                flux_tolerance,
            )

    def _assert_non_zero_displacement_jump(
        self,
        output_data,
        time_steps,
        side_a_node_ids,
        side_b_node_ids,
        min_jump,
    ):
        for time_step in time_steps:
            interface_side_a = GiDOutputFileReader.nodal_values_at_time(
                "TOTAL_DISPLACEMENT", time_step, output_data, side_a_node_ids
            )
            interface_side_b = GiDOutputFileReader.nodal_values_at_time(
                "TOTAL_DISPLACEMENT", time_step, output_data, side_b_node_ids
            )
            max_interface_displacement_jump_y = max(
                abs(a[1] - b[1]) for a, b in zip(interface_side_a, interface_side_b)
            )
            self.assertGreater(max_interface_displacement_jump_y, min_jump)

    def _run_interface_vs_soil_scenario(
        self,
        interface_case_name,
        soil_case_name,
        all_interface_to_soil_min_pairs,
        shared_unique_min_pairs,
        displacement_tolerance,
        pressure_tolerance,
        flux_tolerance,
        side_pair_min_pairs=None,
        min_jump=1e-10,
    ):

        interface_path = test_helper.get_file_path(
            os.path.join("UPw_interface", interface_case_name)
        )
        soil_path = test_helper.get_file_path(
            os.path.join("UPw_interface", soil_case_name)
        )

        interface_output = self._run_two_stage_interface_case(interface_path)
        soil_output = self._run_two_stage_soil_case(soil_path)

        interface_times = self._all_times(interface_output)
        soil_times = self._all_times(soil_output)

        interface_node_ids, soil_node_ids = self._all_interface_to_soil_node_ids(
            interface_path,
            soil_path,
            min_pairs=all_interface_to_soil_min_pairs,
        )

        shared_unique_interface_node_ids, shared_unique_soil_node_ids = (
            self._shared_unique_node_ids(
                interface_path,
                soil_path,
                min_pairs=shared_unique_min_pairs,
            )
        )

        self._assert_nodal_match(
            interface_output,
            soil_output,
            interface_times,
            soil_times,
            interface_node_ids,
            soil_node_ids,
            displacement_tolerance=displacement_tolerance,
            pressure_tolerance=pressure_tolerance,
            flux_tolerance=flux_tolerance,
        )

        self._assert_nodal_match(
            interface_output,
            soil_output,
            interface_times,
            soil_times,
            shared_unique_interface_node_ids,
            shared_unique_soil_node_ids,
            displacement_tolerance=displacement_tolerance,
            pressure_tolerance=pressure_tolerance,
            flux_tolerance=flux_tolerance,
        )

        if side_pair_min_pairs is not None:
            side_a_node_ids, side_b_node_ids = self._interface_node_ids_by_side(
                interface_path,
                min_pairs=side_pair_min_pairs,
            )

            self._assert_non_zero_displacement_jump(
                interface_output,
                interface_times,
                side_a_node_ids=side_a_node_ids,
                side_b_node_ids=side_b_node_ids,
                min_jump=min_jump,
            )

    def test_vertical_interface_matches_column_on_shared_unique_nodes(self):
        self._run_interface_vs_soil_scenario(
            interface_case_name="column_vertical_interface",
            soil_case_name="column",
            all_interface_to_soil_min_pairs=3,
            shared_unique_min_pairs=3,
            displacement_tolerance=1e-6,
            pressure_tolerance=1.71,
            flux_tolerance=7e-6,
        )

    def test_horizontal_interface_matches_column_on_shared_nodes(self):
        self._run_interface_vs_soil_scenario(
            interface_case_name="column_horizontal_interface",
            soil_case_name="column",
            all_interface_to_soil_min_pairs=2,
            shared_unique_min_pairs=4,
            displacement_tolerance=1e-6,
            pressure_tolerance=0.11,
            flux_tolerance=5e-6,
            side_pair_min_pairs=2,
        )

    def test_horizontal_interface_diff_order_matches_column_diff_order(self):
        self._run_interface_vs_soil_scenario(
            interface_case_name="column_horizontal_interface_diff_order_elements",
            soil_case_name="column_diff_order_elements",
            all_interface_to_soil_min_pairs=3,
            shared_unique_min_pairs=10,
            displacement_tolerance=1e-6,
            pressure_tolerance=0.11,
            flux_tolerance=5e-6,
            side_pair_min_pairs=3,
        )


if __name__ == "__main__":
    KratosUnittest.main()
