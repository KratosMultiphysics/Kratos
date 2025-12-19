from KratosMultiphysics.GeoMechanicsApplication.gid_output_file_reader import GiDOutputFileReader
from KratosMultiphysics.GeoMechanicsApplication import unit_conversions
import KratosMultiphysics.KratosUnittest as KratosUnittest
import os
import pathlib

import test_helper

if test_helper.want_test_plots():
    import KratosMultiphysics.GeoMechanicsApplication.geo_plot_utilities as plot_utils


def _extract_x_and_y_from_line(line, index_of_x=0, index_of_y=1):
    words = line.split()
    return (float(words[index_of_x]), float(words[index_of_y]))


def extract_time_and_settlement_from_line(line):
    return _extract_x_and_y_from_line(line, index_of_x=1, index_of_y=2)


def extract_stress_and_y_from_line(line):
    return _extract_x_and_y_from_line(line, index_of_x=2, index_of_y=1)


def make_compression_positive_and_convert_Pa_to_kPa(kratos_stress_components):
    return [
        -1.0 * unit_conversions.Pa_to_kPa(stress) for stress in kratos_stress_components
    ]


def get_nodal_vertical_effective_stress_at_time(time_in_s, output_data, node_ids=None):
    return make_compression_positive_and_convert_Pa_to_kPa(
        [
            stress_vector[1]
            for stress_vector in GiDOutputFileReader.nodal_values_at_time(
                "CAUCHY_STRESS_TENSOR", time_in_s, output_data, node_ids=node_ids
            )
        ]
    )


def get_nodal_vertical_total_stress_at_time(time_in_s, output_data, node_ids=None):
    return make_compression_positive_and_convert_Pa_to_kPa(
        [
            stress_vector[1]
            for stress_vector in GiDOutputFileReader.nodal_values_at_time(
                "TOTAL_STRESS_TENSOR", time_in_s, output_data, node_ids=node_ids
            )
        ]
    )


def get_nodal_water_pressures_at_time(time_in_s, output_data, node_ids=None):
    return make_compression_positive_and_convert_Pa_to_kPa(
        GiDOutputFileReader.nodal_values_at_time(
            "WATER_PRESSURE", time_in_s, output_data, node_ids=node_ids
        )
    )


def extract_nodal_settlement_over_time(output_data, node_id):
    result = []
    for item in output_data["results"]["TOTAL_DISPLACEMENT"]:
        if item["location"] != "OnNodes":
            continue

        settlement = None
        for value_item in item["values"]:
            if value_item["node"] == node_id:
                settlement = -1.0 * value_item["value"][1]
                break
        assert settlement is not None

        result.append((unit_conversions.seconds_to_days(item["time"]), settlement))

    return result


def make_settlement_history_plot(
    stage_outputs, node_ids, path_to_ref_data_points, figure_filename
):
    data_series_collection = []

    data_series_collection.append(
        plot_utils.DataSeries(
            test_helper.get_data_points_from_file(
                path_to_ref_data_points, extract_time_and_settlement_from_line
            ),
            "D-Settlement",
            marker="1",
        )
    )

    data_points_by_node = {node_id: [] for node_id in node_ids}
    for output_data in stage_outputs:
        for node_id in node_ids:
            data_points_by_node[node_id].extend(
                extract_nodal_settlement_over_time(output_data, node_id)
            )
    for node_id in node_ids:
        data_series_collection.append(
            plot_utils.DataSeries(
                data_points_by_node[node_id],
                f"node {node_id} [Kratos]",
                line_style=":",
                marker="+",
            )
        )

    plot_utils.make_settlement_history_plot(data_series_collection, figure_filename)


def get_ref_y_coordinates(post_msh_file_path, node_ids):
    coordinates = test_helper.read_coordinates_from_post_msh_file(
        post_msh_file_path, node_ids=node_ids
    )

    # The top edge of the soil column in the Kratos model is at y = 50.0. In the corresponding D-Settlement model, the
    # same edge is located at y = 0.0. The following lambda transforms the y coordinate of the Kratos model to the
    # corresponding one of the D-Settlement model.
    to_ref_y_coordinate = lambda y: y - 50.0
    return [to_ref_y_coordinate(coord[1]) for coord in coordinates]


def make_stress_over_y_plot(
    output_data,
    time_in_s,
    y_coordinates,
    node_ids_over_depth,
    ref_data,
    plot_file_path,
    want_water_pressure_plot=True,
):
    data_series_collection = []

    # Extract reference data points from files
    for item in ref_data:
        data_points = test_helper.get_data_points_from_file(
            item["file_path"], extract_stress_and_y_from_line
        )
        data_series_collection.append(
            plot_utils.DataSeries(data_points, item["label"], marker="1")
        )

    # Extract data points from the Kratos analysis results
    if want_water_pressure_plot:
        water_pressures = get_nodal_water_pressures_at_time(
            time_in_s, output_data, node_ids=node_ids_over_depth
        )
        data_series_collection.append(
            plot_utils.DataSeries(
                zip(water_pressures, y_coordinates, strict=True),
                r"$p_{\mathrm{w}}$ [Kratos]",
                line_style=":",
                marker="+",
            )
        )

    effective_vertical_stresses = get_nodal_vertical_effective_stress_at_time(
        time_in_s, output_data, node_ids=node_ids_over_depth
    )
    data_series_collection.append(
        plot_utils.DataSeries(
            zip(effective_vertical_stresses, y_coordinates, strict=True),
            r"$\sigma_{\mathrm{eff, yy}}$ [Kratos]",
            line_style=":",
            marker="+",
        )
    )

    total_vertical_stresses = get_nodal_vertical_total_stress_at_time(
        time_in_s, output_data, node_ids=node_ids_over_depth
    )
    data_series_collection.append(
        plot_utils.DataSeries(
            zip(total_vertical_stresses, y_coordinates, strict=True),
            r"$\sigma_{\mathrm{tot, yy}}$ [Kratos]",
            line_style=":",
            marker="+",
        )
    )

    plot_utils.make_stress_over_y_plot(data_series_collection, plot_file_path)


class KratosGeoMechanicsDSettlementValidationTests(KratosUnittest.TestCase):
    def test_settlement_dry_column(self):
        """
        This test validates the settlement of a dry column under uniform load.
        The test runs multiple stages of a settlement simulation and checks the
        settlement values at specific times against expected results.
        The expected settlement values are based on an analytical solution.
        """
        test_name = "dry_column_uniform_load"
        test_root = "dsettlement"
        project_path = test_helper.get_file_path(os.path.join(test_root, test_name, "coarse_mesh"))
        ref_path = test_helper.get_file_path(os.path.join(test_root, test_name))

        import KratosMultiphysics.GeoMechanicsApplication.run_geo_settlement as run_geo_settlement

        n_stages = 5
        project_parameters_filenames = [
            os.path.join("..", "common", f"ProjectParameters_stage{i+1}.json") for i in range(n_stages)
        ]
        status = run_geo_settlement.run_stages(
            project_path, project_parameters_filenames
        )
        self.assertEqual(status, 0)

        project_path = pathlib.Path(project_path)
        ref_path = pathlib.Path(ref_path)

        reader = GiDOutputFileReader()

        output_stage_3 = reader.read_output_from(project_path / "stage3.post.res")
        output_stage_5 = reader.read_output_from(project_path / "stage5.post.res")

        if test_helper.want_test_plots():
            output_stage_4 = reader.read_output_from(project_path / "stage4.post.res")
            top_node_ids = [2, 3, 104]
            make_settlement_history_plot(
                (output_stage_3, output_stage_4, output_stage_5),
                top_node_ids,
                ref_path / "ref_settlement_data.txt",
                project_path / "test_case_1_settlement_plot.svg",
            )

            left_side_corner_node_ids = [3] + list(range(105, 154)) + [4]
            ref_y_coordinates = get_ref_y_coordinates(
                project_path / "stage1.post.msh", left_side_corner_node_ids
            )

            # Make a stress plot at the start of the analysis
            ref_data = [
                {
                    "file_path": ref_path
                    / "ref_effective_vertical_stresses_after_0.1_days.txt",
                    "label": r"$\sigma_{\mathrm{eff, yy}}$ [D-Settlement]",
                },
                {
                    "file_path": ref_path
                    / "ref_total_vertical_stresses_after_0.1_days.txt",
                    "label": r"$\sigma_{\mathrm{tot, yy}}$ [D-Settlement]",
                },
            ]
            make_stress_over_y_plot(
                output_stage_3,
                unit_conversions.days_to_seconds(0.1) + 1.0,
                ref_y_coordinates,
                left_side_corner_node_ids,
                ref_data,
                project_path / "test_case_1_stress_plot_after_0.1_days.svg",
                want_water_pressure_plot=False,
            )

            # Make a stress plot after 100.1 days
            ref_data = [
                {
                    "file_path": ref_path
                    / "ref_effective_vertical_stresses_after_100.1_days.txt",
                    "label": r"$\sigma_{\mathrm{eff, yy}}$ [D-Settlement]",
                },
                {
                    "file_path": ref_path
                    / "ref_total_vertical_stresses_after_100.1_days.txt",
                    "label": r"$\sigma_{\mathrm{tot, yy}}$ [D-Settlement]",
                },
            ]
            make_stress_over_y_plot(
                output_stage_5,
                unit_conversions.days_to_seconds(100.1) + 1.0,
                ref_y_coordinates,
                left_side_corner_node_ids,
                ref_data,
                project_path / "test_case_1_stress_plot_after_100.1_days.svg",
                want_water_pressure_plot=False,
            )

            # Make a stress plot at the end of the fifth stage (when consolidation is supposed to be finished)
            ref_data = [
                {
                    "file_path": ref_path
                    / "ref_effective_vertical_stresses_after_10000_days.txt",
                    "label": r"$\sigma_{\mathrm{eff, yy}}$ [D-Settlement]",
                },
                {
                    "file_path": ref_path
                    / "ref_total_vertical_stresses_after_10000_days.txt",
                    "label": r"$\sigma_{\mathrm{tot, yy}}$ [D-Settlement]",
                },
            ]
            make_stress_over_y_plot(
                output_stage_5,
                unit_conversions.days_to_seconds(10000),
                ref_y_coordinates,
                left_side_corner_node_ids,
                ref_data,
                project_path / "test_case_1_stress_plot_after_10000_days.svg",
                want_water_pressure_plot=False,
            )

        # Check some results
        top_middle_node_id = 104
        actual_settlement_after_one_hundred_days = reader.nodal_values_at_time(
            "TOTAL_DISPLACEMENT", 8640000, output_stage_3, [top_middle_node_id]
        )[0][1]
        expected_settlement_after_one_hundred_days = -3.22

        # Assert the value to be within 1% of the analytical solution
        self.assertTrue(
            abs(
                (
                    expected_settlement_after_one_hundred_days
                    - actual_settlement_after_one_hundred_days
                )
                / expected_settlement_after_one_hundred_days
            )
            < 0.01
        )

        actual_settlement_after_ten_thousand_days = reader.nodal_values_at_time(
            "TOTAL_DISPLACEMENT", 864000000, output_stage_5, [top_middle_node_id]
        )[0][1]

        expected_settlement_after_ten_thousand_days = -8.01

        # Assert the value to be within 1% of the analytical solution
        self.assertTrue(
            abs(
                (
                    expected_settlement_after_ten_thousand_days
                    - actual_settlement_after_ten_thousand_days
                )
                / actual_settlement_after_ten_thousand_days
            )
            < 0.01
        )

    def test_settlement_dry_column_fine_mesh(self):
        """
        This test validates the settlement of a dry column under uniform load.
        The test runs multiple stages of a settlement simulation and checks the
        settlement values at specific times against expected results.
        The expected settlement values are based on an analytical solution.
        """
        test_name = "dry_column_uniform_load"
        test_root = "dsettlement"
        project_path = test_helper.get_file_path(os.path.join(test_root, test_name, "fine_mesh"))
        ref_path = test_helper.get_file_path(os.path.join(test_root, test_name))

        import KratosMultiphysics.GeoMechanicsApplication.run_geo_settlement as run_geo_settlement

        n_stages = 5
        project_parameters_filenames = [
            f"../common/ProjectParameters_stage{i+1}.json" for i in range(n_stages)
        ]
        status = run_geo_settlement.run_stages(
            project_path, project_parameters_filenames
        )
        self.assertEqual(status, 0)

        project_path = pathlib.Path(project_path)
        ref_path = pathlib.Path(ref_path)

        reader = GiDOutputFileReader()

        output_stage_3 = reader.read_output_from(project_path / "stage3.post.res")
        output_stage_5 = reader.read_output_from(project_path / "stage5.post.res")

        if test_helper.want_test_plots():
            output_stage_4 = reader.read_output_from(project_path / "stage4.post.res")
            top_node_ids = [2, 3, 1008]
            make_settlement_history_plot(
                (output_stage_3, output_stage_4, output_stage_5),
                top_node_ids,
                ref_path / "ref_settlement_data.txt",
                project_path / "test_case_1_settlement_plot.svg",
                )

            left_side_corner_node_ids = [3] + list(range(1023, 2021)) + [4]

            # Taking every 10th node, to improve clarity of the plots for this model with a finer mesh.
            left_side_corner_node_ids = left_side_corner_node_ids[0::10]

            ref_y_coordinates = get_ref_y_coordinates(
                project_path / "stage1.post.msh", left_side_corner_node_ids
            )

            # Make a stress plot at the start of the analysis
            ref_data = [
                {
                    "file_path": ref_path
                                 / "ref_effective_vertical_stresses_after_0.1_days.txt",
                    "label": r"$\sigma_{\mathrm{eff, yy}}$ [D-Settlement]",
                },
                {
                    "file_path": ref_path
                                 / "ref_total_vertical_stresses_after_0.1_days.txt",
                    "label": r"$\sigma_{\mathrm{tot, yy}}$ [D-Settlement]",
                },
            ]
            make_stress_over_y_plot(
                output_stage_3,
                unit_conversions.days_to_seconds(0.1) + 1.0,
                ref_y_coordinates,
                left_side_corner_node_ids,
                ref_data,
                project_path / "test_case_1_stress_plot_after_0.1_days.svg",
                want_water_pressure_plot=False,
                )

            # Make a stress plot after 100.1 days
            ref_data = [
                {
                    "file_path": ref_path
                                 / "ref_effective_vertical_stresses_after_100.1_days.txt",
                    "label": r"$\sigma_{\mathrm{eff, yy}}$ [D-Settlement]",
                },
                {
                    "file_path": ref_path
                                 / "ref_total_vertical_stresses_after_100.1_days.txt",
                    "label": r"$\sigma_{\mathrm{tot, yy}}$ [D-Settlement]",
                },
            ]
            make_stress_over_y_plot(
                output_stage_5,
                unit_conversions.days_to_seconds(100.1) + 1.0,
                ref_y_coordinates,
                left_side_corner_node_ids,
                ref_data,
                project_path / "test_case_1_stress_plot_after_100.1_days.svg",
                want_water_pressure_plot=False,
                )

            # Make a stress plot at the end of the fifth stage (when consolidation is supposed to be finished)
            ref_data = [
                {
                    "file_path": ref_path
                                 / "ref_effective_vertical_stresses_after_10000_days.txt",
                    "label": r"$\sigma_{\mathrm{eff, yy}}$ [D-Settlement]",
                },
                {
                    "file_path": ref_path
                                 / "ref_total_vertical_stresses_after_10000_days.txt",
                    "label": r"$\sigma_{\mathrm{tot, yy}}$ [D-Settlement]",
                },
            ]
            make_stress_over_y_plot(
                output_stage_5,
                unit_conversions.days_to_seconds(10000),
                ref_y_coordinates,
                left_side_corner_node_ids,
                ref_data,
                project_path / "test_case_1_stress_plot_after_10000_days.svg",
                want_water_pressure_plot=False,
                )

        # Check some results
        top_middle_node_id = 1008
        actual_settlement_after_one_hundred_days = reader.nodal_values_at_time(
            "TOTAL_DISPLACEMENT", 8640000, output_stage_3, [top_middle_node_id]
        )[0][1]
        expected_settlement_after_one_hundred_days = -3.22

        # Assert the value to be within 1% of the analytical solution
        self.assertTrue(
            abs(
                (
                        expected_settlement_after_one_hundred_days
                        - actual_settlement_after_one_hundred_days
                )
                / expected_settlement_after_one_hundred_days
            )
            < 0.01
        )

        actual_settlement_after_ten_thousand_days = reader.nodal_values_at_time(
            "TOTAL_DISPLACEMENT", 864000000, output_stage_5, [top_middle_node_id]
        )[0][1]

        expected_settlement_after_ten_thousand_days = -8.01

        # Assert the value to be within 1% of the analytical solution
        self.assertTrue(
            abs(
                (
                        expected_settlement_after_ten_thousand_days
                        - actual_settlement_after_ten_thousand_days
                )
                / actual_settlement_after_ten_thousand_days
            )
            < 0.01
        )

    def test_settlement_consolidation_coarse_mesh(self):
        """
        This test validates the settlement of a fully saturated column under uniform load.
        The test runs multiple stages of a settlement simulation and checks the
        pore pressure and stresses at specific time but several locations against expected results.
        The expected settlement values are based on regression tests.
        """
        test_name = "coarse_mesh"
        test_root = "dsettlement/consolidation_uniform_load"
        project_path = test_helper.get_file_path(os.path.join(test_root, test_name))

        import KratosMultiphysics.GeoMechanicsApplication.run_geo_settlement as run_geo_settlement

        n_stages = 3
        project_parameters_filenames = [
            f"ProjectParameters_stage{i+1}.json" for i in range(n_stages)
        ]
        status = run_geo_settlement.run_stages(
            project_path, project_parameters_filenames
        )
        self.assertEqual(status, 0)

        reader = GiDOutputFileReader()

        output_data = reader.read_output_from(
            os.path.join(project_path, "stage3.post.res")
        )

        nodes = [3, 15, 16, 17, 18, 4]
        expected_pressure = [0.0, -18524.3, -24096.0, -33456.2, -42816.6, -50000.0]
        expected_total_stress = [
            -5332.28,
            -24858.1,
            -52279.1,
            -70030.9,
            -89823.8,
            -109627.0,
        ]

        for i in range(6):
            actual_pressure = reader.nodal_values_at_time(
                "WATER_PRESSURE", 8726.4, output_data, [nodes[i]]
            )[0]
            self.assertAlmostEqual(actual_pressure, expected_pressure[i], 4)

            actual_total_stress = reader.nodal_values_at_time(
                "TOTAL_STRESS_TENSOR", 8726.4, output_data, [nodes[i]]
            )[0][1]
            self.assertAlmostEqual(actual_total_stress, expected_total_stress[i], 4)

    def test_settlement_consolidation_fine_mesh(self):
        """
        This test validates the settlement of a fully saturated column under uniform load.
        The test runs multiple stages of a settlement simulation and checks the
        pore pressure and stresses at specific time but several locations against expected results.
        The expected settlement values are based on regression tests.
        """
        test_name = "fine_mesh"
        test_root = "dsettlement/consolidation_uniform_load"
        project_path = test_helper.get_file_path(os.path.join(test_root, test_name))

        import KratosMultiphysics.GeoMechanicsApplication.run_geo_settlement as run_geo_settlement

        n_stages = 3
        project_parameters_filenames = [
            f"ProjectParameters_stage{i+1}.json" for i in range(n_stages)
        ]
        status = run_geo_settlement.run_stages(
            project_path, project_parameters_filenames
        )
        self.assertEqual(status, 0)

        reader = GiDOutputFileReader()

        output_data = reader.read_output_from(
            os.path.join(project_path, "stage3.post.res")
        )

        nodes = [87, 154, 266, 390, 516, 641]
        expected_pressure = [0.0, -16382.0, -24416.0, -33477.7, -42701.9, -50000.0]
        expected_total_stress = [
            -10324.4,
            -30053.2,
            -50021.0,
            -70014.2,
            -89998.3,
            -110028.0,
        ]

        for i in range(6):
            actual_pressure = reader.nodal_values_at_time(
                "WATER_PRESSURE", 8726.4, output_data, [nodes[i]]
            )[0]
            self.assertAlmostEqual(actual_pressure, expected_pressure[i], 4)

            actual_total_stress = reader.nodal_values_at_time(
                "TOTAL_STRESS_TENSOR", 8726.4, output_data, [nodes[i]]
            )[0][1]
            self.assertAlmostEqual(actual_total_stress, expected_total_stress[i], 4)

    def test_settlement_fully_saturated_column(self):
        """
        This test validates the settlement of a fully saturated column under uniform load.
        The test runs multiple stages of a settlement simulation and checks the
        settlement values at specific times against expected results.
        The expected settlement values are based on an analytical solution.
        """
        test_name = "high_permeability"
        test_root = os.path.join("dsettlement", "fully_saturated_column_uniform_load")
        project_path = test_helper.get_file_path(os.path.join(test_root, test_name))

        import KratosMultiphysics.GeoMechanicsApplication.run_geo_settlement as run_geo_settlement

        n_stages = 5
        project_parameters_filenames = [
            os.path.join("..", "common", f"ProjectParameters_stage{i+1}.json") for i in range(n_stages)
        ]
        status = run_geo_settlement.run_stages(
            project_path, project_parameters_filenames
        )
        self.assertEqual(status, 0)

        project_path = pathlib.Path(project_path)

        reader = GiDOutputFileReader()

        output_stage_3 = reader.read_output_from(project_path / "stage3.post.res")
        output_stage_5 = reader.read_output_from(project_path / "stage5.post.res")

        if test_helper.want_test_plots():
            output_stage_2 = reader.read_output_from(project_path / "stage2.post.res")
            output_stage_4 = reader.read_output_from(project_path / "stage4.post.res")
            top_node_ids = [2, 3, 104]
            make_settlement_history_plot(
                (output_stage_3, output_stage_4, output_stage_5),
                top_node_ids,
                project_path / "ref_settlement_data.txt",
                project_path / "test_case_2_settlement_plot.svg",
            )

            left_side_corner_node_ids = [3] + list(range(105, 154)) + [4]
            ref_y_coordinates = get_ref_y_coordinates(
                project_path / "stage1.post.msh", left_side_corner_node_ids
            )

            # Make a stress plot at the start of the analysis
            ref_data = [
                {
                    "file_path": project_path / "ref_water_pressures_after_0_days.txt",
                    "label": r"$p_{\mathrm{w}}$ [D-Settlement]",
                },
                {
                    "file_path": project_path
                    / "ref_effective_vertical_stresses_after_0_days.txt",
                    "label": r"$\sigma_{\mathrm{eff, yy}}$ [D-Settlement]",
                },
                {
                    "file_path": project_path
                    / "ref_total_vertical_stresses_after_0_days.txt",
                    "label": r"$\sigma_{\mathrm{tot, yy}}$ [D-Settlement]",
                },
            ]
            make_stress_over_y_plot(
                output_stage_2,
                unit_conversions.days_to_seconds(0) + 1.0,
                ref_y_coordinates,
                left_side_corner_node_ids,
                ref_data,
                project_path / "test_case_2_stress_plot_after_0_days.svg",
            )

            # Make a stress plot after 100 days of consolidation have passed
            ref_data = [
                {
                    "file_path": project_path
                    / "ref_water_pressures_after_100_days.txt",
                    "label": r"$p_{\mathrm{w}}$ [D-Settlement]",
                },
                {
                    "file_path": project_path
                    / "ref_effective_vertical_stresses_after_100_days.txt",
                    "label": r"$\sigma_{\mathrm{eff, yy}}$ [D-Settlement]",
                },
                {
                    "file_path": project_path
                    / "ref_total_vertical_stresses_after_100_days.txt",
                    "label": r"$\sigma_{\mathrm{tot, yy}}$ [D-Settlement]",
                },
            ]
            make_stress_over_y_plot(
                output_stage_3,
                unit_conversions.days_to_seconds(100),
                ref_y_coordinates,
                left_side_corner_node_ids,
                ref_data,
                project_path / "test_case_2_stress_plot_after_100_days.svg",
            )

            # Make a stress plot after applying the surface load
            ref_data = [
                {
                    "file_path": project_path
                    / "ref_water_pressures_after_100.1_days.txt",
                    "label": r"$p_{\mathrm{w}}$ [D-Settlement]",
                },
                {
                    "file_path": project_path
                    / "ref_effective_vertical_stresses_after_100.1_days.txt",
                    "label": r"$\sigma_{\mathrm{eff, yy}}$ [D-Settlement]",
                },
                {
                    "file_path": project_path
                    / "ref_total_vertical_stresses_after_100.1_days.txt",
                    "label": r"$\sigma_{\mathrm{tot, yy}}$ [D-Settlement]",
                },
            ]
            make_stress_over_y_plot(
                output_stage_5,
                unit_conversions.days_to_seconds(100.1) + 1.0,
                ref_y_coordinates,
                left_side_corner_node_ids,
                ref_data,
                project_path / "test_case_2_stress_plot_after_100.1_days.svg",
            )

            # Make a stress plot at the end of the fifth stage (when consolidation is supposed to be finished)
            ref_data = [
                {
                    "file_path": project_path
                    / "ref_water_pressures_after_10000_days.txt",
                    "label": r"$p_{\mathrm{w}}$ [D-Settlement]",
                },
                {
                    "file_path": project_path
                    / "ref_effective_vertical_stresses_after_10000_days.txt",
                    "label": r"$\sigma_{\mathrm{eff, yy}}$ [D-Settlement]",
                },
                {
                    "file_path": project_path
                    / "ref_total_vertical_stresses_after_10000_days.txt",
                    "label": r"$\sigma_{\mathrm{tot, yy}}$ [D-Settlement]",
                },
            ]
            make_stress_over_y_plot(
                output_stage_5,
                unit_conversions.days_to_seconds(10000),
                ref_y_coordinates,
                left_side_corner_node_ids,
                ref_data,
                project_path / "test_case_2_stress_plot_after_10000_days.svg",
            )

        # Check some results
        actual_settlement_after_one_hundred_days = reader.nodal_values_at_time(
            "TOTAL_DISPLACEMENT", unit_conversions.days_to_seconds(100), output_stage_3, [104]
        )[0][1]
        self.assertAlmostEqual(actual_settlement_after_one_hundred_days, -1.70997, 4)

        actual_settlement_after_ten_thousand_days = reader.nodal_values_at_time(
            "TOTAL_DISPLACEMENT", unit_conversions.days_to_seconds(10000), output_stage_5, [104]
        )[0][1]
        self.assertTrue(
            abs((-8.63753 - actual_settlement_after_ten_thousand_days) / -8.63753) < 0.01
        )

    def test_settlement_phreatic_line_below_surface(self):
        """
        This test validates the settlement of a soil column where the phreatic line
        is 10 m below the soil surface.
        """
        test_name = "phreatic_line_below_soil_surface"
        test_root = "dsettlement"
        project_path = test_helper.get_file_path(os.path.join(test_root, test_name))

        import KratosMultiphysics.GeoMechanicsApplication.run_geo_settlement as run_geo_settlement

        n_stages = 5
        project_parameters_filenames = [
            f"ProjectParameters_stage{i+1}.json" for i in range(n_stages)
        ]
        status = run_geo_settlement.run_stages(
            project_path, project_parameters_filenames
        )
        self.assertEqual(status, 0)

        reader = GiDOutputFileReader()
        top_node_ids = [2, 3, 104]
        project_path = pathlib.Path(project_path)

        output_stage_3 = reader.read_output_from(project_path / "stage3.post.res")
        output_stage_5 = reader.read_output_from(project_path / "stage5.post.res")

        if test_helper.want_test_plots():
            output_stage_4 = reader.read_output_from(project_path / "stage4.post.res")
            make_settlement_history_plot(
                (output_stage_3, output_stage_4, output_stage_5),
                top_node_ids,
                project_path / "ref_settlement_data.txt",
                project_path / "test_case_3_settlement_plot.svg",
            )

            left_side_corner_node_ids = [3] + list(range(105, 154)) + [4]
            ref_y_coordinates = get_ref_y_coordinates(
                project_path / "stage1.post.msh", left_side_corner_node_ids
            )

            # Make a stress plot at the start of the analysis
            output_stage_2 = reader.read_output_from(project_path / "stage2.post.res")
            ref_data = [
                {
                    "file_path": project_path / "ref_water_pressures_after_0_days.txt",
                    "label": r"$p_{\mathrm{w}}$ [D-Settlement]",
                },
                {
                    "file_path": project_path
                    / "ref_effective_vertical_stresses_after_0_days.txt",
                    "label": r"$\sigma_{\mathrm{eff, yy}}$ [D-Settlement]",
                },
                {
                    "file_path": project_path
                    / "ref_total_vertical_stresses_after_0_days.txt",
                    "label": r"$\sigma_{\mathrm{tot, yy}}$ [D-Settlement]",
                },
            ]
            make_stress_over_y_plot(
                output_stage_2,
                unit_conversions.days_to_seconds(0) + 1.0,
                ref_y_coordinates,
                left_side_corner_node_ids,
                ref_data,
                project_path / "test_case_3_stress_plot_after_0_days.svg",
            )

            # Make a stress plot at the end of the third stage
            ref_data = [
                {
                    "file_path": project_path
                    / "ref_water_pressures_after_100_days.txt",
                    "label": r"$p_{\mathrm{w}}$ [D-Settlement]",
                },
                {
                    "file_path": project_path
                    / "ref_effective_vertical_stresses_after_100_days.txt",
                    "label": r"$\sigma_{\mathrm{eff, yy}}$ [D-Settlement]",
                },
                {
                    "file_path": project_path
                    / "ref_total_vertical_stresses_after_100_days.txt",
                    "label": r"$\sigma_{\mathrm{tot, yy}}$ [D-Settlement]",
                },
            ]
            make_stress_over_y_plot(
                output_stage_3,
                unit_conversions.days_to_seconds(100),
                ref_y_coordinates,
                left_side_corner_node_ids,
                ref_data,
                project_path / "test_case_3_stress_plot_after_100_days.svg",
            )

            # Make a stress plot at the start of the fifth stage
            ref_data = [
                {
                    "file_path": project_path
                    / "ref_water_pressures_after_100.1_days.txt",
                    "label": r"$p_{\mathrm{w}}$ [D-Settlement]",
                },
                {
                    "file_path": project_path
                    / "ref_effective_vertical_stresses_after_100.1_days.txt",
                    "label": r"$\sigma_{\mathrm{eff, yy}}$ [D-Settlement]",
                },
                {
                    "file_path": project_path
                    / "ref_total_vertical_stresses_after_100.1_days.txt",
                    "label": r"$\sigma_{\mathrm{tot, yy}}$ [D-Settlement]",
                },
            ]
            make_stress_over_y_plot(
                output_stage_5,
                unit_conversions.days_to_seconds(100.1) + 1.0,
                ref_y_coordinates,
                left_side_corner_node_ids,
                ref_data,
                project_path / "test_case_3_stress_plot_after_100.1_days.svg",
            )

            # Make a stress plot at the end of the fifth stage (when consolidation is supposed to be finished)
            ref_data = [
                {
                    "file_path": project_path
                    / "ref_water_pressures_after_10000_days.txt",
                    "label": r"$p_{\mathrm{w}}$ [D-Settlement]",
                },
                {
                    "file_path": project_path
                    / "ref_effective_vertical_stresses_after_10000_days.txt",
                    "label": r"$\sigma_{\mathrm{eff, yy}}$ [D-Settlement]",
                },
                {
                    "file_path": project_path
                    / "ref_total_vertical_stresses_after_10000_days.txt",
                    "label": r"$\sigma_{\mathrm{tot, yy}}$ [D-Settlement]",
                },
            ]
            make_stress_over_y_plot(
                output_stage_5,
                unit_conversions.days_to_seconds(10000),
                ref_y_coordinates,
                left_side_corner_node_ids,
                ref_data,
                project_path / "test_case_3_stress_plot_after_10000_days.svg",
            )

        # Check some settlement values
        check_data = [
            {
                "output_data": output_stage_3,
                "time_in_s": unit_conversions.days_to_seconds(0.1) + 1.0,
                "expected_total_u_y": 0.0,  # analytical value
                "delta": 0.02,
            },
            {
                "output_data": output_stage_3,
                "time_in_s": 129601,
                "expected_total_u_y": -0.057,  # regression value
                "delta": 0.02,
            },
            {
                "output_data": output_stage_3,
                "time_in_s": 1097281,
                "expected_total_u_y": -0.46,  # regression value
                "delta": 0.02,
            },
            {
                "output_data": output_stage_3,
                "time_in_s": unit_conversions.days_to_seconds(100),
                "expected_total_u_y": -1.75,  # analytical value
                "delta": 0.06,
            },
            {
                "output_data": output_stage_5,
                "time_in_s": 17487361,
                "expected_total_u_y": -3.69,  # regression value
                "delta": 0.02,
            },
            {
                "output_data": output_stage_5,
                "time_in_s": 79418881,
                "expected_total_u_y": -5.43,  # regression value
                "delta": 0.02,
            },
            {
                "output_data": output_stage_5,
                "time_in_s": unit_conversions.days_to_seconds(10000),
                "expected_total_u_y": -7.90,  # analytical value
                "delta": 0.12,
            },
        ]
        for item in check_data:
            actual_total_displacement_of_top_edge = reader.nodal_values_at_time(
                "TOTAL_DISPLACEMENT",
                item["time_in_s"],
                item["output_data"],
                top_node_ids,
            )
            for total_displacement_vector, node_id in zip(
                actual_total_displacement_of_top_edge, top_node_ids, strict=True
            ):
                self.assertAlmostEqual(
                    total_displacement_vector[1],
                    item["expected_total_u_y"],
                    places=None,
                    delta=item["delta"],
                    msg=f"total vertical displacement at node {node_id} at time {item['time_in_s']} [s]",
                )

    def test_settlement_fully_saturated_column_low_permeability(self):
        """
        This test validates the settlement of a fully saturated column under uniform load.
        The test runs multiple stages of a settlement simulation and checks the
        settlement values at specific times against expected results.
        The expected settlement values are partly based on an analytical solution, partly
        on regression values.
        """
        test_name = "low_permeability"
        test_root = os.path.join("dsettlement", "fully_saturated_column_uniform_load")
        project_path = test_helper.get_file_path(os.path.join(test_root, test_name))

        import KratosMultiphysics.GeoMechanicsApplication.run_geo_settlement as run_geo_settlement

        n_stages = 5
        project_parameters_filenames = [
            os.path.join("..", "common", f"ProjectParameters_stage{i+1}.json") for i in range(n_stages)
        ]
        status = run_geo_settlement.run_stages(
            project_path, project_parameters_filenames
        )
        self.assertEqual(status, 0)

        project_path = pathlib.Path(project_path)

        reader = GiDOutputFileReader()

        output_stage_3 = reader.read_output_from(project_path / "stage3.post.res")
        output_stage_5 = reader.read_output_from(project_path / "stage5.post.res")

        if test_helper.want_test_plots():
            output_stage_2 = reader.read_output_from(project_path / "stage2.post.res")
            output_stage_4 = reader.read_output_from(project_path / "stage4.post.res")
            top_node_ids = [2, 3, 104]
            make_settlement_history_plot(
                (output_stage_3, output_stage_4, output_stage_5),
                top_node_ids,
                project_path / "ref_settlement_data.txt",
                project_path / "test_case_4_settlement_plot.svg",
            )

            left_side_corner_node_ids = [3] + list(range(105, 154)) + [4]
            ref_y_coordinates = get_ref_y_coordinates(
                project_path / "stage1.post.msh", left_side_corner_node_ids
            )

            # Make a stress plot at the start of the analysis
            ref_data = [
                {
                    "file_path": project_path / "ref_water_pressures_after_0_days.txt",
                    "label": r"$p_{\mathrm{w}}$ [D-Settlement]",
                },
                {
                    "file_path": project_path
                    / "ref_effective_vertical_stresses_after_0_days.txt",
                    "label": r"$\sigma_{\mathrm{eff, yy}}$ [D-Settlement]",
                },
                {
                    "file_path": project_path
                    / "ref_total_vertical_stresses_after_0_days.txt",
                    "label": r"$\sigma_{\mathrm{tot, yy}}$ [D-Settlement]",
                },
            ]
            make_stress_over_y_plot(
                output_stage_2,
                unit_conversions.days_to_seconds(0) + 1.0,
                ref_y_coordinates,
                left_side_corner_node_ids,
                ref_data,
                project_path / "test_case_4_stress_plot_after_0_days.svg",
            )

            # Make a stress plot after 100 days of consolidation have passed
            ref_data = [
                {
                    "file_path": project_path
                    / "ref_water_pressures_after_100_days.txt",
                    "label": r"$p_{\mathrm{w}}$ [D-Settlement]",
                },
                {
                    "file_path": project_path
                    / "ref_effective_vertical_stresses_after_100_days.txt",
                    "label": r"$\sigma_{\mathrm{eff, yy}}$ [D-Settlement]",
                },
                {
                    "file_path": project_path
                    / "ref_total_vertical_stresses_after_100_days.txt",
                    "label": r"$\sigma_{\mathrm{tot, yy}}$ [D-Settlement]",
                },
            ]
            make_stress_over_y_plot(
                output_stage_3,
                unit_conversions.days_to_seconds(100),
                ref_y_coordinates,
                left_side_corner_node_ids,
                ref_data,
                project_path / "test_case_4_stress_plot_after_100_days.svg",
            )

            # Make a stress plot after applying the surface load
            ref_data = [
                {
                    "file_path": project_path
                    / "ref_water_pressures_after_100.1_days.txt",
                    "label": r"$p_{\mathrm{w}}$ [D-Settlement]",
                },
                {
                    "file_path": project_path
                    / "ref_effective_vertical_stresses_after_100.1_days.txt",
                    "label": r"$\sigma_{\mathrm{eff, yy}}$ [D-Settlement]",
                },
                {
                    "file_path": project_path
                    / "ref_total_vertical_stresses_after_100.1_days.txt",
                    "label": r"$\sigma_{\mathrm{tot, yy}}$ [D-Settlement]",
                },
            ]
            make_stress_over_y_plot(
                output_stage_5,
                unit_conversions.days_to_seconds(100.1) + 1.0,
                ref_y_coordinates,
                left_side_corner_node_ids,
                ref_data,
                project_path / "test_case_4_stress_plot_after_100.1_days.svg",
            )

            # Make a stress plot at the end of the fifth stage (when consolidation is supposed to be finished)
            ref_data = [
                {
                    "file_path": project_path
                    / "ref_water_pressures_after_10000_days.txt",
                    "label": r"$p_{\mathrm{w}}$ [D-Settlement]",
                },
                {
                    "file_path": project_path
                    / "ref_effective_vertical_stresses_after_10000_days.txt",
                    "label": r"$\sigma_{\mathrm{eff, yy}}$ [D-Settlement]",
                },
                {
                    "file_path": project_path
                    / "ref_total_vertical_stresses_after_10000_days.txt",
                    "label": r"$\sigma_{\mathrm{tot, yy}}$ [D-Settlement]",
                },
            ]
            make_stress_over_y_plot(
                output_stage_5,
                unit_conversions.days_to_seconds(10000),
                ref_y_coordinates,
                left_side_corner_node_ids,
                ref_data,
                project_path / "test_case_4_stress_plot_after_10000_days.svg",
            )

        # Check some results
        actual_settlement_after_one_hundred_days = reader.nodal_values_at_time(
            "TOTAL_DISPLACEMENT", unit_conversions.days_to_seconds(100), output_stage_3, [104]
        )[0][1]
        self.assertAlmostEqual(
            actual_settlement_after_one_hundred_days, -0.495277, 4
        )  # Regression value

        actual_settlement_after_ten_thousand_days = reader.nodal_values_at_time(
            "TOTAL_DISPLACEMENT", unit_conversions.days_to_seconds(10000), output_stage_5, [104]
        )[0][1]
        # Assert the value to be within 1% of the analytical solution
        self.assertTrue(
            abs((-8.48 - actual_settlement_after_ten_thousand_days) / -8.48) < 0.01
        )


if __name__ == "__main__":
    KratosUnittest.main()
