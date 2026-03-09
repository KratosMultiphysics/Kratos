import os
import re
import shutil
import tempfile
from pathlib import Path

import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_helper
import KratosMultiphysics.GeoMechanicsApplication.run_geo_settlement as run_geo_settlement
from KratosMultiphysics.GeoMechanicsApplication.gid_output_file_reader import GiDOutputFileReader

class KratosGeoMechanicsUPwInterfaceTests(KratosUnittest.TestCase):
    """
    This class contains benchmark tests which check if interfaces between soil and structural elements are
    calculated correctly
    """

    def setUp(self):
        self.n_stages = 2
        self.project_parameters_filenames = [
            os.path.join(f"ProjectParameters_stage{i+1}.json") for i in range(self.n_stages)
        ]
        self.interface_project_parameters_filenames = [
            os.path.join("..", "common", f"ProjectParameters_interface_stage{i+1}.json")
            for i in range(self.n_stages)
        ]
        self.soil_project_parameters_filenames = [
            os.path.join("..", "common", f"ProjectParameters_stage{i+1}.json")
            for i in range(self.n_stages)
        ]
        self.output_reader = GiDOutputFileReader()

    def tearDown(self):
        pass

    def _run_two_stage_case(self, case_path, project_parameters_filenames=None):
        project_parameters_filenames = (
            project_parameters_filenames or self.project_parameters_filenames
        )
        status = run_geo_settlement.run_stages(case_path, project_parameters_filenames)
        self.assertEqual(status, 0)
        return self._read_output(case_path, "stage2")

    def _run_two_stage_interface_case(self, case_path):
        return self._run_two_stage_case(
            case_path, self.interface_project_parameters_filenames
        )

    def _run_two_stage_soil_case(self, case_path):
        return self._run_two_stage_case(
            case_path, self.soil_project_parameters_filenames
        )

    def _read_output(self, case_path, stage_name):
        output_file_name = os.path.join(case_path, f"{stage_name}.post.res")
        return self.output_reader.read_output_from(output_file_name)

    @staticmethod
    def _last_time(output_data):
        times = GiDOutputFileReader.get_time_steps_from_first_valid_result(output_data)
        if not times:
            raise RuntimeError("No valid output time steps found")
        return times[-1]

    @staticmethod
    def _max_component_difference(vectors_a, vectors_b, component_index):
        return max(
            abs(vector_a[component_index] - vector_b[component_index])
            for vector_a, vector_b in zip(vectors_a, vectors_b)
        )

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
        first_nodes = KratosGeoMechanicsNewInterfaceTests._read_mdpa_nodes(first_mdpa_file_path)
        second_nodes = KratosGeoMechanicsNewInterfaceTests._read_mdpa_nodes(second_mdpa_file_path)

        first_coordinate_to_node_ids_map = KratosGeoMechanicsNewInterfaceTests._coordinate_to_node_ids_map(first_nodes)
        second_coordinate_to_node_ids_map = KratosGeoMechanicsNewInterfaceTests._coordinate_to_node_ids_map(second_nodes)

        shared_coordinates = sorted(
            set(first_coordinate_to_node_ids_map.keys()).intersection(second_coordinate_to_node_ids_map.keys())
        )

        node_id_pairs = []
        for coordinates in shared_coordinates:
            first_node_ids = first_coordinate_to_node_ids_map[coordinates]
            second_node_ids = second_coordinate_to_node_ids_map[coordinates]

            if len(first_node_ids) == 1 and len(second_node_ids) == 1:
                node_id_pairs.append((first_node_ids[0], second_node_ids[0]))

        return node_id_pairs

    @staticmethod
    def _interface_unique_node_id_pairs(mdpa_file_path):
        node_id_pairs = []
        seen_pairs = set()
        reading_interface_elements = False

        with open(mdpa_file_path, "r", encoding="utf-8") as mdpa_file:
            for line in mdpa_file:
                stripped_line = line.strip()

                if stripped_line.startswith("Begin Elements"):
                    element_name = stripped_line[len("Begin Elements") :].strip()
                    reading_interface_elements = "Interface" in element_name
                    continue

                if stripped_line == "End Elements":
                    reading_interface_elements = False
                    continue

                if not reading_interface_elements or not stripped_line:
                    continue

                columns = stripped_line.split()
                if len(columns) < 4:
                    continue

                element_node_ids = [int(value) for value in columns[2:]]
                if len(element_node_ids) < 2 or len(element_node_ids) % 2 != 0:
                    continue

                half = len(element_node_ids) // 2
                for node_id_a, node_id_b in zip(element_node_ids[:half], element_node_ids[half:]):
                    pair = (node_id_a, node_id_b)
                    if pair not in seen_pairs:
                        seen_pairs.add(pair)
                        node_id_pairs.append(pair)

        return node_id_pairs

    def _run_interface_and_soil_cases(self, interface_case_name, soil_case_name):
        interface_path = test_helper.get_file_path(
            os.path.join('UPw_interface', interface_case_name)
        )
        soil_path = test_helper.get_file_path(
            os.path.join('UPw_interface', soil_case_name)
        )

        interface_output = self._run_two_stage_interface_case(interface_path)
        soil_output = self._run_two_stage_soil_case(soil_path)

        interface_time = self._last_time(interface_output)
        soil_time = self._last_time(soil_output)

        return interface_path, soil_path, interface_output, soil_output, interface_time, soil_time

    def _shared_unique_node_ids(self, interface_path, soil_path, min_pairs):
        node_id_pairs = self._shared_unique_node_id_pairs(
            os.path.join(interface_path, "column.mdpa"),
            os.path.join(soil_path, "column.mdpa"),
        )
        self.assertGreaterEqual(len(node_id_pairs), min_pairs)

        interface_node_ids = [pair[0] for pair in node_id_pairs]
        soil_node_ids = [pair[1] for pair in node_id_pairs]

        return interface_node_ids, soil_node_ids

    def _all_interface_to_soil_node_ids(self, interface_path, soil_path, min_pairs):
        interface_mdpa_file_path = os.path.join(interface_path, "column.mdpa")
        soil_mdpa_file_path = os.path.join(soil_path, "column.mdpa")

        interface_node_id_pairs = self._interface_unique_node_id_pairs(interface_mdpa_file_path)
        self.assertGreaterEqual(len(interface_node_id_pairs), min_pairs)

        interface_node_ids = list(dict.fromkeys(
            node_id
            for interface_node_id_pair in interface_node_id_pairs
            for node_id in interface_node_id_pair
        ))

        interface_nodes = self._read_mdpa_nodes(interface_mdpa_file_path)
        soil_nodes = self._read_mdpa_nodes(soil_mdpa_file_path)
        soil_coordinate_to_node_ids_map = self._coordinate_to_node_ids_map(soil_nodes)

        soil_node_ids = []
        for interface_node_id in interface_node_ids:
            self.assertIn(interface_node_id, interface_nodes)

            coordinate_key = tuple(round(value, 12) for value in interface_nodes[interface_node_id])
            matching_soil_node_ids = soil_coordinate_to_node_ids_map.get(coordinate_key, [])
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

    def _assert_nodal_match(
        self,
        interface_output,
        soil_output,
        interface_time,
        soil_time,
        interface_node_ids,
        soil_node_ids,
        displacement_tolerance,
        pressure_tolerance,
    ):
        interface_displacements = GiDOutputFileReader.nodal_values_at_time(
            "TOTAL_DISPLACEMENT", interface_time, interface_output, interface_node_ids
        )
        soil_displacements = GiDOutputFileReader.nodal_values_at_time(
            "TOTAL_DISPLACEMENT", soil_time, soil_output, soil_node_ids
        )
        max_displacement_difference_y = self._max_component_difference(
            interface_displacements, soil_displacements, component_index=1
        )
        self.assertLess(max_displacement_difference_y, displacement_tolerance)

        interface_pressures = GiDOutputFileReader.nodal_values_at_time(
            "WATER_PRESSURE", interface_time, interface_output, interface_node_ids
        )
        soil_pressures = GiDOutputFileReader.nodal_values_at_time(
            "WATER_PRESSURE", soil_time, soil_output, soil_node_ids
        )
        max_pressure_difference = max(
            abs(interface_pressure - soil_pressure)
            for interface_pressure, soil_pressure in zip(interface_pressures, soil_pressures)
        )
        self.assertLess(max_pressure_difference, pressure_tolerance)

    def _assert_non_zero_displacement_jump(
        self,
        output_data,
        time_step,
        side_a_node_ids,
        side_b_node_ids,
        min_jump,
    ):
        interface_side_a = GiDOutputFileReader.nodal_values_at_time(
            "TOTAL_DISPLACEMENT", time_step, output_data, side_a_node_ids
        )
        interface_side_b = GiDOutputFileReader.nodal_values_at_time(
            "TOTAL_DISPLACEMENT", time_step, output_data, side_b_node_ids
        )
        max_interface_displacement_jump_y = self._max_component_difference(
            interface_side_a, interface_side_b, component_index=1
        )
        self.assertGreater(max_interface_displacement_jump_y, min_jump)

    def test_horizontal_interface(self):
        file_path = test_helper.get_file_path(os.path.join('UPw_interface', 'column'))
        output_data = self._run_two_stage_soil_case(file_path)
        self.assertTrue(output_data.get("results"))

    def test_horizontal_interface_diff_order(self):
        file_path = test_helper.get_file_path(os.path.join('UPw_interface', 'column_diff_order_elements'))
        output_data = self._run_two_stage_soil_case(file_path)
        self.assertTrue(output_data.get("results"))

    def test_vertical_interface(self):
        file_path = test_helper.get_file_path(os.path.join('UPw_interface', 'column_vertical_interface'))
        output_data = self._run_two_stage_interface_case(file_path)
        self.assertTrue(output_data.get("results"))

    def test_vertical_interface_matches_column_on_shared_unique_nodes(self):
        interface_path, soil_path, interface_output, soil_output, interface_time, soil_time = self._run_interface_and_soil_cases(
            'column_vertical_interface', 'column'
        )
        interface_node_ids, soil_node_ids = self._all_interface_to_soil_node_ids(
            interface_path, soil_path, min_pairs=3
        )
        shared_unique_interface_node_ids, shared_unique_soil_node_ids = self._shared_unique_node_ids(
            interface_path, soil_path, min_pairs=3
        )

        self._assert_nodal_match(
            interface_output,
            soil_output,
            interface_time,
            soil_time,
            interface_node_ids,
            soil_node_ids,
            displacement_tolerance=1e-6,
            pressure_tolerance=0.1,
        )

        self._assert_nodal_match(
            interface_output,
            soil_output,
            interface_time,
            soil_time,
            shared_unique_interface_node_ids,
            shared_unique_soil_node_ids,
            displacement_tolerance=1e-6,
            pressure_tolerance=0.1,
        )

    def test_horizontal_interface_matches_column_on_shared_nodes(self):
        interface_path, soil_path, interface_output, soil_output, interface_time, soil_time = self._run_interface_and_soil_cases(
            'column_horizontal_interface', 'column'
        )
        interface_node_ids, soil_node_ids = self._all_interface_to_soil_node_ids(
            interface_path, soil_path, min_pairs=2
        )
        shared_unique_interface_node_ids, shared_unique_soil_node_ids = self._shared_unique_node_ids(
            interface_path, soil_path, min_pairs=4
        )
        side_a_node_ids, side_b_node_ids = self._interface_node_ids_by_side(
            interface_path, min_pairs=2
        )

        self._assert_nodal_match(
            interface_output,
            soil_output,
            interface_time,
            soil_time,
            interface_node_ids,
            soil_node_ids,
            displacement_tolerance=5e-7,
            pressure_tolerance=0.1,
        )

        self._assert_nodal_match(
            interface_output,
            soil_output,
            interface_time,
            soil_time,
            shared_unique_interface_node_ids,
            shared_unique_soil_node_ids,
            displacement_tolerance=5e-7,
            pressure_tolerance=0.1,
        )

        self._assert_non_zero_displacement_jump(
            interface_output,
            interface_time,
            side_a_node_ids=side_a_node_ids,
            side_b_node_ids=side_b_node_ids,
            min_jump=1e-10,
        )

    def test_horizontal_interface_diff_order_matches_column_diff_order(self):
        interface_path, soil_path, interface_output, soil_output, interface_time, soil_time = self._run_interface_and_soil_cases(
            'column_horizontal_interface_diff_order_elements', 'column_diff_order_elements'
        )

        interface_node_ids, soil_node_ids = self._all_interface_to_soil_node_ids(
            interface_path, soil_path, min_pairs=3
        )
        shared_unique_interface_node_ids, shared_unique_soil_node_ids = self._shared_unique_node_ids(
            interface_path, soil_path, min_pairs=10
        )
        side_a_node_ids, side_b_node_ids = self._interface_node_ids_by_side(
            interface_path, min_pairs=3
        )

        self._assert_nodal_match(
            interface_output,
            soil_output,
            interface_time,
            soil_time,
            interface_node_ids,
            soil_node_ids,
            displacement_tolerance=5e-7,
            pressure_tolerance=0.1,
        )

        self._assert_nodal_match(
            interface_output,
            soil_output,
            interface_time,
            soil_time,
            shared_unique_interface_node_ids,
            shared_unique_soil_node_ids,
            displacement_tolerance=5e-7,
            pressure_tolerance=0.1,
        )

        self._assert_non_zero_displacement_jump(
            interface_output,
            interface_time,
            side_a_node_ids=side_a_node_ids,
            side_b_node_ids=side_b_node_ids,
            min_jump=1e-10,
        )

    def test_top_pressure_table_changes_interface_water_pressure(self):
        file_path = test_helper.get_file_path(os.path.join('UPw_interface', 'column_horizontal_interface'))

        baseline_output = self._run_two_stage_interface_case(file_path)
        baseline_time = self._last_time(baseline_output)
        baseline_top_pressure = GiDOutputFileReader.nodal_values_at_time(
            "WATER_PRESSURE", baseline_time, baseline_output, [5]
        )[0]

        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_case = Path(tmp_dir) / "column_horizontal_interface"
            shutil.copytree(file_path, tmp_case)

            source_upw_interface_dir = Path(file_path).parent
            source_common_dir = source_upw_interface_dir / "common"
            for stage_index in range(1, self.n_stages + 1):
                shutil.copyfile(
                    source_common_dir / f"ProjectParameters_interface_stage{stage_index}.json",
                    tmp_case / f"ProjectParameters_stage{stage_index}.json",
                )

            tmp_common_dir = tmp_case.parent / "common"
            tmp_common_dir.mkdir(exist_ok=True)
            shutil.copyfile(
                source_common_dir / "TwoMaterialParameters.json",
                tmp_common_dir / "TwoMaterialParameters.json",
            )

            mdpa_path = tmp_case / "column.mdpa"
            mdpa_text = mdpa_path.read_text(encoding="utf-8")
            mdpa_text_updated = re.sub(
                r"(Begin Table 2 TIME WATER_PRESSURE\s+0\.0000000000\s+0\.0000000000\s+)1\.0000000000\s+1000\.0000000000",
                r"\g<1>1.0000000000 0.0000000000",
                mdpa_text,
                flags=re.MULTILINE,
            )
            self.assertNotEqual(mdpa_text, mdpa_text_updated)
            mdpa_path.write_text(mdpa_text_updated, encoding="utf-8")

            no_pressure_output = self._run_two_stage_case(str(tmp_case))

        no_pressure_time = self._last_time(no_pressure_output)
        no_pressure_top_pressure = GiDOutputFileReader.nodal_values_at_time(
            "WATER_PRESSURE", no_pressure_time, no_pressure_output, [5]
        )[0]

        self.assertGreater(baseline_top_pressure - no_pressure_top_pressure, 900.0)
        self.assertAlmostEqual(no_pressure_top_pressure, 0.0, delta=1.0)

if __name__ == '__main__':
    KratosUnittest.main()
    
