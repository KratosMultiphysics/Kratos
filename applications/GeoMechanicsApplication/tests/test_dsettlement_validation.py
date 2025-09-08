from KratosMultiphysics.GeoMechanicsApplication import run_multiple_stages
from KratosMultiphysics.GeoMechanicsApplication import unit_conversions
import KratosMultiphysics.KratosUnittest as KratosUnittest
import os
import pathlib

import test_helper

if test_helper.want_test_plots():
    import KratosMultiphysics.GeoMechanicsApplication.geo_plot_utilities as plot_utils


def get_data_points_from_file(file_path, line_to_data_point):
    result = []
    with open(file_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            result.append(line_to_data_point(line))

    return result


def _extract_x_and_y_from_line(line, index_of_x=0, index_of_y=1):
    words = line.split()
    return (float(words[index_of_x]), float(words[index_of_y]))


def extract_time_and_settlement_from_line(line):
    return _extract_x_and_y_from_line(line, index_of_x=1, index_of_y=2)


def extract_stress_and_y_from_line(line):
    return _extract_x_and_y_from_line(line, index_of_x=2, index_of_y=1)


def extract_nodal_settlement_over_time(output_data, node_id):
    result = []
    for item in output_data["results"]["TOTAL_DISPLACEMENT"]:
        if item["location"] != "OnNodes":
            continue

        total_y_displacement = None
        for value_item in item["values"]:
            if value_item["node"] == node_id:
                total_y_displacement = -1.0 * value_item["value"][1]
                break
        assert total_y_displacement is not None

        result.append((unit_conversions.seconds_to_days(item["time"]), total_y_displacement))

    return result


def get_nodal_vertical_stress_component_at_time(stress_item_name, time_in_seconds, output_data, node_ids=None):
    stress_vectors = test_helper.GiDOutputFileReader.nodal_values_at_time(stress_item_name, time_in_seconds, output_data, node_ids=node_ids)
    # Invert the sign of the vertical stress component such that compression becomes positive.
    return [-1.0 * unit_conversions.Pa_to_kPa(vector[1]) for vector in stress_vectors]


def get_nodal_vertical_effective_stress_at_time(time_in_seconds, output_data, node_ids=None):
    return get_nodal_vertical_stress_component_at_time("CAUCHY_STRESS_TENSOR", time_in_seconds, output_data, node_ids=node_ids)


def get_nodal_water_pressures_at_time(time_in_seconds, output_data, node_ids=None):
    water_pressures = test_helper.GiDOutputFileReader.nodal_values_at_time("WATER_PRESSURE", time_in_seconds, output_data, node_ids=node_ids)
    # Invert the sign of the water pressure such that compression becomes positive.
    return [-1.0 * unit_conversions.Pa_to_kPa(value) for value in water_pressures]


def shift_y_of_kratos_model(y):
    # The top edge of the soil column in the Kratos model is at y = 50.0. In the corresponding "old" D-Settlement model,
    # the same edge is located at y = 0.0. This function transforms the y coordinate of the Kratos model to the
    # corresponding one of the "old" D-Settlement model.
    return y - 50.0


def make_settlement_plot(stage_outputs, node_ids, path_to_ref_data_points, figure_filename):
    data_points_by_node = {node_id : [] for node_id in node_ids}
    for output_data in stage_outputs:
        for node_id in node_ids:
            data_points_by_node[node_id].extend(extract_nodal_settlement_over_time(output_data, node_id))

    data_series_collection = []
    data_series_collection.append(plot_utils.DataSeries(get_data_points_from_file(path_to_ref_data_points, extract_time_and_settlement_from_line), 'ref', marker='+'))
    for node_id in node_ids:
        data_series_collection.append(plot_utils.DataSeries(data_points_by_node[node_id], f'node {node_id}', line_style=':', marker='+'))

    plot_utils.plot_settlement_results(data_series_collection, figure_filename)


class StressPlotDataFilePaths:
    def __init__(self):
        self.path_to_water_pressure_data = None
        self.path_to_vertical_effective_stress_data = None


def make_stress_over_depth_plot(output_data, time_in_sec, post_msh_file_path, node_ids_over_depth, ref_data, plot_file_path):
    data_series_collection = []

    # Extract reference data points from files
    if ref_data.path_to_water_pressure_data:
        data_points = get_data_points_from_file(ref_data.path_to_water_pressure_data, extract_stress_and_y_from_line)
        data_series_collection.append(plot_utils.DataSeries(data_points, 'ref P_w', marker='+'))

    if ref_data.path_to_vertical_effective_stress_data:
        data_points = get_data_points_from_file(ref_data.path_to_vertical_effective_stress_data, extract_stress_and_y_from_line)
        data_series_collection.append(plot_utils.DataSeries(data_points, 'ref sigma_yy;eff', marker='+'))

    # Extract data points from the Kratos analysis results
    coordinates = test_helper.read_coordinates_from_post_msh_file(post_msh_file_path, node_ids=node_ids_over_depth)
    y_coordinates = [shift_y_of_kratos_model(coord[1]) for coord in coordinates]

    water_pressures = get_nodal_water_pressures_at_time(time_in_sec, output_data, node_ids=node_ids_over_depth)
    data_series_collection.append(plot_utils.DataSeries(zip(water_pressures, y_coordinates, strict=True), 'P_w [Kratos]', line_style=':', marker='+'))

    effective_vertical_stresses = get_nodal_vertical_effective_stress_at_time(time_in_sec, output_data, node_ids=node_ids_over_depth)
    data_series_collection.append(plot_utils.DataSeries(zip(effective_vertical_stresses, y_coordinates, strict=True), 'sigma_yy;eff [Kratos]', line_style=':', marker='+'))

    plot_utils.make_stress_plot(data_series_collection, plot_file_path)


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

        project_path = pathlib.Path(project_path)

        reader = test_helper.GiDOutputFileReader()

        output_stage_3 = reader.read_output_from(project_path / "stage3.post.res")
        top_middle_node_id = 104
        actual_settlement_after_one_hundred_days = reader.nodal_values_at_time(
            "TOTAL_DISPLACEMENT", 8640000, output_stage_3, [top_middle_node_id]
        )[0][1]
        expected_settlement_after_one_hundred_days = -3.22

        # Assert the value to be within 1% of the analytical solution
        self.assertTrue(
            abs((
                expected_settlement_after_one_hundred_days
                - actual_settlement_after_one_hundred_days
            )
            / expected_settlement_after_one_hundred_days)
            < 0.01
        )

        output_stage_5 = reader.read_output_from(project_path / "stage5.post.res")
        actual_settlement_after_ten_thousand_days = reader.nodal_values_at_time(
            "TOTAL_DISPLACEMENT", 864000000, output_stage_5, [top_middle_node_id]
        )[0][1]

        expected_settlement_after_ten_thousand_days = -8.01

        # Assert the value to be within 1% of the analytical solution
        self.assertTrue(
            abs((
                expected_settlement_after_ten_thousand_days
                - actual_settlement_after_ten_thousand_days
            )
            / actual_settlement_after_ten_thousand_days)
            < 0.01
        )

        if test_helper.want_test_plots():
            output_stage_4 = reader.read_output_from(project_path / "stage4.post.res")
            top_node_ids = [2, 3, 104]
            make_settlement_plot((output_stage_3, output_stage_4, output_stage_5), top_node_ids, project_path / "ref_settlement_data.txt", project_path / "test_case_1_settlement_plot.svg")


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

        reader = test_helper.GiDOutputFileReader()

        output_data = reader.read_output_from(
            os.path.join(project_path, "stage3.post.res")
        )
        
        nodes = [3, 15, 16, 17, 18, 4]
        expected_pressure = [0.0, -18524.3, -24096.0, -33456.2, -42816.6, -50000.0]
        expected_total_stress = [-5332.28, -24858.1, -52279.1, -70030.9, -89823.8, -109627.0]
        
        for i in range(6):
            actual_pressure = reader.nodal_values_at_time("WATER_PRESSURE", 8726.4, output_data, [nodes[i]])[0]
            self.assertAlmostEqual(actual_pressure, expected_pressure[i], 4)
        
            actual_total_stress = reader.nodal_values_at_time(
                "TOTAL_STRESS_TENSOR", 8726.4, output_data, [nodes[i]])[0][1]
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

        reader = test_helper.GiDOutputFileReader()

        output_data = reader.read_output_from(
            os.path.join(project_path, "stage3.post.res")
        )
        
        nodes = [87, 154, 266, 390, 516, 641]
        expected_pressure = [0.0, -16382.0, -24416.0, -33477.7, -42701.9, -50000.0]
        expected_total_stress = [-10324.4, -30053.2, -50021.0, -70014.2, -89998.3, -110028.0]
        
        for i in range(6):
            actual_pressure = reader.nodal_values_at_time("WATER_PRESSURE", 8726.4, output_data, [nodes[i]])[0]
            self.assertAlmostEqual(actual_pressure, expected_pressure[i], 4)
        
            actual_total_stress = reader.nodal_values_at_time(
                "TOTAL_STRESS_TENSOR", 8726.4, output_data, [nodes[i]])[0][1]
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
            f"ProjectParameters_stage{i+1}.json" for i in range(n_stages)
        ]
        status = run_geo_settlement.run_stages(
            project_path, project_parameters_filenames
        )
        self.assertEqual(status, 0)

        project_path = pathlib.Path(project_path)

        reader = test_helper.GiDOutputFileReader()

        output_stage_2 = reader.read_output_from(project_path / "stage2.post.res")
        actual_settlement_after_one_hundred_days = reader.nodal_values_at_time(
            "TOTAL_DISPLACEMENT", 8640000, output_stage_2, [104]
        )[0][1]
        self.assertAlmostEqual(actual_settlement_after_one_hundred_days, -1.71094, 4)

        output_stage_5 = reader.read_output_from(project_path / "stage5.post.res")

        actual_settlement_after_ten_thousand_days = reader.nodal_values_at_time(
            "TOTAL_DISPLACEMENT", 864000000, output_stage_5, [104]
        )[0][1]
        self.assertAlmostEqual(actual_settlement_after_ten_thousand_days, -8.63753, 4)

        if test_helper.want_test_plots():
            output_stage_3 = reader.read_output_from(project_path / "stage3.post.res")
            output_stage_4 = reader.read_output_from(project_path / "stage4.post.res")
            top_node_ids = [2, 3, 104]
            make_settlement_plot((output_stage_2, output_stage_3, output_stage_4, output_stage_5), top_node_ids, project_path / "ref_settlement_data.txt", project_path / "test_case_2_settlement_plot.svg")

            ref_data = StressPlotDataFilePaths()
            left_side_corner_node_ids = [3] + list(range(105, 154)) + [4]

            # Make a stress plot after 100 days of consolidation have passed
            ref_data.path_to_water_pressure_data = project_path / "ref_water_pressures_after_100_days.txt"
            ref_data.path_to_vertical_effective_stress_data = project_path / "ref_effective_vertical_stresses_after_100_days.txt"
            make_stress_over_depth_plot(output_stage_2, unit_conversions.days_to_seconds(100), project_path / "stage2.post.msh", left_side_corner_node_ids, ref_data, project_path / "test_case_2_stress_plot_after_100_days.svg")

            # Make a stress plot after applying the surface load
            ref_data.path_to_water_pressure_data = project_path / "ref_water_pressures_after_100.1_days.txt"
            ref_data.path_to_vertical_effective_stress_data = project_path / "ref_effective_vertical_stresses_after_100.1_days.txt"
            make_stress_over_depth_plot(output_stage_4, unit_conversions.days_to_seconds(100.1001), project_path / "stage4.post.msh", left_side_corner_node_ids, ref_data, project_path / "test_case_2_stress_plot_after_100.1_days.svg")

            # Make a stress plot at the end of the fifth stage (when consolidation is supposed to be finished)
            ref_data.path_to_water_pressure_data = project_path / "ref_water_pressures_after_10000_days.txt"
            ref_data.path_to_vertical_effective_stress_data = project_path / "ref_effective_vertical_stresses_after_10000_days.txt"
            make_stress_over_depth_plot(output_stage_5, unit_conversions.days_to_seconds(10000), project_path / "stage5.post.msh", left_side_corner_node_ids, ref_data, project_path / "test_case_2_stress_plot_after_10000_days.svg")


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

        reader = test_helper.GiDOutputFileReader()
        top_node_ids = [2, 3, 104]
        project_path = pathlib.Path(project_path)

        # Check some settlement values
        output_stage_3 = reader.read_output_from(project_path / "stage3.post.res")
        output_stage_5 = reader.read_output_from(project_path / "stage5.post.res")

        check_data = [{"output_data": output_stage_3, "time_in_s": unit_conversions.days_to_seconds(0.1) + 1.0, "expected_total_u_y": 0.0, "delta": 0.01},
                      {"output_data": output_stage_3, "time_in_s": unit_conversions.days_to_seconds(100), "expected_total_u_y": -1.75, "delta": 0.06},
                      {"output_data": output_stage_5, "time_in_s": unit_conversions.days_to_seconds(10000), "expected_total_u_y": -7.90, "delta": 0.12}]
        for item in check_data:
            actual_total_displacement_of_top_edge = reader.nodal_values_at_time(
                "TOTAL_DISPLACEMENT", item["time_in_s"], item["output_data"], top_node_ids
            )
            for total_displacement_vector, node_id in zip(actual_total_displacement_of_top_edge, top_node_ids, strict=True):
                self.assertAlmostEqual(total_displacement_vector[1], item["expected_total_u_y"], places=None, delta=item["delta"], msg=f"total vertical displacement at node {node_id} at time {item["time_in_s"]} [s]")

        if test_helper.want_test_plots():
            left_side_corner_node_ids = [3] + list(range(105, 154)) + [4]
            output_stage_4 = reader.read_output_from(project_path / "stage4.post.res")
            make_settlement_plot((output_stage_3, output_stage_4, output_stage_5), top_node_ids, project_path / "ref_settlement_data.txt", project_path / "test_case_3_settlement_plot.svg")

            ref_data = StressPlotDataFilePaths()

            # Make a stress plot at the end of the third stage
            ref_data.path_to_water_pressure_data = project_path / "ref_water_pressures_after_100_days.txt"
            ref_data.path_to_vertical_effective_stress_data = project_path / "ref_effective_vertical_stresses_after_100_days.txt"
            make_stress_over_depth_plot(output_stage_3, unit_conversions.days_to_seconds(100), project_path / "stage3.post.msh", left_side_corner_node_ids, ref_data, project_path / "test_case_3_stress_plot_after_100_days.svg")

            # Make a stress plot at the start of the fifth stage
            ref_data.path_to_water_pressure_data = project_path / "ref_water_pressures_after_100.1_days.txt"
            ref_data.path_to_vertical_effective_stress_data = project_path / "ref_effective_vertical_stresses_after_100.1_days.txt"
            make_stress_over_depth_plot(output_stage_5, unit_conversions.days_to_seconds(100.1) + 1.0, project_path / "stage5.post.msh", left_side_corner_node_ids, ref_data, project_path / "test_case_3_stress_plot_after_100.1_days.svg")

            # Make a stress plot at the end of the fifth stage (when consolidation is supposed to be finished)
            ref_data.path_to_water_pressure_data = project_path / "ref_water_pressures_after_10000_days.txt"
            ref_data.path_to_vertical_effective_stress_data = project_path / "ref_effective_vertical_stresses_after_10000_days.txt"
            make_stress_over_depth_plot(output_stage_5, unit_conversions.days_to_seconds(10000), project_path / "stage5.post.msh", left_side_corner_node_ids, ref_data, project_path / "test_case_3_stress_plot_after_10000_days.svg")

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
            f"ProjectParameters_stage{i+1}.json" for i in range(n_stages)
        ]
        status = run_geo_settlement.run_stages(
            project_path, project_parameters_filenames
        )
        self.assertEqual(status, 0)

        project_path = pathlib.Path(project_path)

        reader = test_helper.GiDOutputFileReader()

        output_stage_2 = reader.read_output_from(project_path / "stage2.post.res")
        actual_settlement_after_one_hundred_days = reader.nodal_values_at_time(
            "TOTAL_DISPLACEMENT", 8640000, output_stage_2, [104]
        )[0][1]
        self.assertAlmostEqual(actual_settlement_after_one_hundred_days, -0.496382, 4) # Regression value

        output_stage_5 = reader.read_output_from(project_path / "stage5.post.res")
        actual_settlement_after_ten_thousand_days = reader.nodal_values_at_time(
            "TOTAL_DISPLACEMENT", 864000000, output_stage_5, [104]
        )[0][1]

        # Assert the value to be within 1% of the analytical solution
        self.assertTrue(
            abs((
                    -8.48
                    - actual_settlement_after_ten_thousand_days
            )
            / -8.48)
            < 0.01
        )

        if test_helper.want_test_plots():
            output_stage_3 = reader.read_output_from(project_path / "stage3.post.res")
            output_stage_4 = reader.read_output_from(project_path / "stage4.post.res")
            top_node_ids = [2, 3, 104]
            make_settlement_plot((output_stage_2, output_stage_3, output_stage_4, output_stage_5), top_node_ids, project_path / "ref_settlement_data.txt", project_path / "test_case_4_settlement_plot.svg")

            ref_data = StressPlotDataFilePaths()
            left_side_corner_node_ids = [3] + list(range(105, 154)) + [4]

            # Make a stress plot after 100 days of consolidation have passed
            ref_data.path_to_water_pressure_data = project_path / "ref_water_pressures_after_100_days.txt"
            ref_data.path_to_vertical_effective_stress_data = project_path / "ref_effective_vertical_stresses_after_100_days.txt"
            make_stress_over_depth_plot(output_stage_2, unit_conversions.days_to_seconds(100), project_path / "stage2.post.msh", left_side_corner_node_ids, ref_data, project_path / "test_case_4_stress_plot_after_100_days.svg")

            # Make a stress plot after applying the surface load
            ref_data.path_to_water_pressure_data = project_path / "ref_water_pressures_after_100.1_days.txt"
            ref_data.path_to_vertical_effective_stress_data = project_path / "ref_effective_vertical_stresses_after_100.1_days.txt"
            make_stress_over_depth_plot(output_stage_4, unit_conversions.days_to_seconds(100.1001), project_path / "stage4.post.msh", left_side_corner_node_ids, ref_data, project_path / "test_case_4_stress_plot_after_100.1_days.svg")

            # Make a stress plot at the end of the fifth stage (when consolidation is supposed to be finished)
            ref_data.path_to_water_pressure_data = project_path / "ref_water_pressures_after_10000_days.txt"
            ref_data.path_to_vertical_effective_stress_data = project_path / "ref_effective_vertical_stresses_after_10000_days.txt"
            make_stress_over_depth_plot(output_stage_5, unit_conversions.days_to_seconds(10000), project_path / "stage5.post.msh", left_side_corner_node_ids, ref_data, project_path / "test_case_4_stress_plot_after_10000_days.svg")

      
if __name__ == "__main__":
    KratosUnittest.main()
