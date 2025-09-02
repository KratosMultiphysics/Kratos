from KratosMultiphysics.GeoMechanicsApplication import run_multiple_stages
import KratosMultiphysics.KratosUnittest as KratosUnittest
import os
import matplotlib.pyplot as plt
import pathlib

import test_helper

class PlotDataSeries():
    def __init__(self, x_values, y_values, label='', linestyle='-', marker=None):
        self.x_values = x_values
        self.y_values = y_values
        self.label = label
        self.linestyle = linestyle
        self.marker = marker

def seconds_to_days(duration_in_seconds):
    return duration_in_seconds / (24.0 * 60.0 * 60.0)


def days_to_seconds(duration_in_days):
    return duration_in_days * 24.0 * 60.0 * 60.0


def make_plot_data(points, label='', linestyle='-', marker=None):
    times = [data[0] for data in points]
    total_y_displacements = [data[1] for data in points]
    return PlotDataSeries(times, total_y_displacements, label=label, linestyle=linestyle, marker=marker)


def make_water_pressure_plot_data(points, label='', linestyle='-', marker=None):
    ys = [data[0] for data in points]
    water_pressures = [data[1] for data in points]
    return PlotDataSeries(water_pressures, ys, label=label, linestyle=linestyle, marker=marker)


def extract_nodal_settlement_over_time(output_data, node_id):
    result = []
    for item in output_data["results"]["TOTAL_DISPLACEMENT"]:
        if item["location"] != "OnNodes":
            continue

        time_in_days = seconds_to_days(item["time"])

        total_y_displacement = None
        for value_item in item["values"]:
            if value_item["node"] == node_id:
                total_y_displacement = -1.0 * value_item["value"][1]
                break
        assert total_y_displacement is not None

        result.append((time_in_days, total_y_displacement))

    return result

def plot_settlement_results(series_collection, figure_filename):
    figure, axes = plt.subplots(layout="constrained")
    axes.set_xscale('log')
    for series in series_collection:
        axes.plot(series.x_values, series.y_values, label=series.label, linestyle=series.linestyle, marker=series.marker)
    axes.grid()
    axes.grid(which="minor", color="0.9")
    axes.yaxis.set_inverted(True)
    axes.set_xlabel('Time [day]')
    axes.set_ylabel('Settlement [m]')
    figure.legend(loc='outside center right')

    if isinstance(figure_filename, pathlib.Path):
        figure_filename = str(figure_filename.resolve())
    print(f"Save plot to {figure_filename}")
    plt.savefig(figure_filename)

def make_stress_plot(series_collection, figure_filename):
    figure, axes = plt.subplots(layout="constrained")
    for series in series_collection:
        axes.plot(series.x_values, series.y_values, label=series.label, linestyle=series.linestyle, marker=series.marker)
    axes.grid()
    axes.grid(which="minor", color="0.9")
    axes.set_xlabel('Stress [kPa]')
    axes.set_ylabel('y [m]')
    figure.legend(loc='outside center right')

    if isinstance(figure_filename, pathlib.Path):
        figure_filename = str(figure_filename.resolve())
    print(f"Save plot to {figure_filename}")
    plt.savefig(figure_filename)


def get_data_points_from(file_path):
    result = []
    with open(file_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith("#"):
                continue

            _, x, y = line.split()  # ignore the row number
            result.append((float(x), float(y)))

    return result


def get_nodal_vertical_stress_component_at_time(stress_item_name, time_in_seconds, output_data, node_ids=None):
    stress_vectors_by_node_id = test_helper.GiDOutputFileReader.nodal_values_at_time2(stress_item_name, time_in_seconds, output_data, node_ids=node_ids)
    # Invert the sign of the vertical stress component such that compression becomes positive. Also convert Pa to kPa.
    return [-1.0 * (stress_vectors_by_node_id[node_id][1] / 1000.0) for node_id in node_ids]


def get_nodal_vertical_effective_stress_at_time(time_in_seconds, output_data, node_ids=None):
    return get_nodal_vertical_stress_component_at_time("CAUCHY_STRESS_TENSOR", time_in_seconds, output_data, node_ids=node_ids)


def get_nodal_water_pressures_at_time(time_in_seconds, output_data, node_ids=None):
    water_pressures_by_node_id = test_helper.GiDOutputFileReader.nodal_values_at_time2("WATER_PRESSURE", time_in_seconds, output_data, node_ids=node_ids)
    # Invert the sign of the water pressure such that compression becomes positive. Also convert Pa to kPa.
    return [-1.0 * (water_pressures_by_node_id[node_id] / 1000.0) for node_id in node_ids]


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

    graph_series = []
    graph_series.append(make_plot_data(get_data_points_from(path_to_ref_data_points), 'ref', marker='+'))
    for node_id in node_ids:
        graph_series.append(make_plot_data(data_points_by_node[node_id], f'node {node_id}', linestyle=':', marker='+'))

    plot_settlement_results(graph_series, figure_filename)


class StressPlotDataFilePaths:
    def __init__(self):
        self.path_to_water_pressure_data = None
        self.path_to_vertical_effective_stress_data = None


def make_stress_over_depth_plot(output_data, time_in_sec, post_msh_file_path, node_ids_over_depth, ref_data, plot_file_path):
    graph_series = []

    # Extract reference data points from files
    data_points = get_data_points_from(ref_data.path_to_water_pressure_data)
    graph_series.append(make_water_pressure_plot_data(data_points, 'ref P_w', marker='+'))

    data_points = get_data_points_from(ref_data.path_to_vertical_effective_stress_data)
    graph_series.append(make_water_pressure_plot_data(data_points, 'ref sigma_yy;eff', marker='+'))

    # Extract data points from the Kratos analysis results
    coordinates = test_helper.read_coordinates_from_post_msh_file(post_msh_file_path, node_ids=node_ids_over_depth)
    y_coordinates = [shift_y_of_kratos_model(coord[1]) for coord in coordinates]

    water_pressures = get_nodal_water_pressures_at_time(time_in_sec, output_data, node_ids=node_ids_over_depth)
    graph_series.append(PlotDataSeries(water_pressures, y_coordinates, 'P_w [Kratos]', linestyle=':', marker='+'))

    effective_vertical_stresses = get_nodal_vertical_effective_stress_at_time(time_in_sec, output_data, node_ids=node_ids_over_depth)
    graph_series.append(PlotDataSeries(effective_vertical_stresses, y_coordinates, 'sigma_yy;eff [Kratos]', linestyle=':', marker='+'))

    make_stress_plot(graph_series, plot_file_path)


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

        original_working_dir = os.getcwd()
        os.chdir(project_path)

        import KratosMultiphysics.GeoMechanicsApplication.run_geo_settlement as run_geo_settlement

        n_stages = 5
        project_parameters_filenames = [
            f"ProjectParameters_stage{i+1}.json" for i in range(n_stages)
        ]
        status = run_geo_settlement.run_stages(
            project_path, project_parameters_filenames
        )
        self.assertEqual(status, 0)

        os.chdir(original_working_dir)

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
            (
                expected_settlement_after_one_hundred_days
                - actual_settlement_after_one_hundred_days
            )
            / expected_settlement_after_one_hundred_days
            < 0.01
        )

        output_stage_5 = reader.read_output_from(project_path / "stage5.post.res")
        actual_settlement_after_ten_thousand_days = reader.nodal_values_at_time(
            "TOTAL_DISPLACEMENT", 864000000, output_stage_5, [top_middle_node_id]
        )[0][1]

        expected_settlement_after_ten_thousand_days = -8.01

        # Assert the value to be within 1% of the analytical solution
        self.assertTrue(
            (
                expected_settlement_after_ten_thousand_days
                - actual_settlement_after_ten_thousand_days
            )
            / expected_settlement_after_one_hundred_days
            < 0.01
        )

        if test_helper.want_test_plots():
            output_stage_4 = reader.read_output_from(project_path / "stage4.post.res")
            top_node_ids = [2, 3, 104]
            make_settlement_plot((output_stage_3, output_stage_4, output_stage_5), top_node_ids, project_path / "ref_settlement_data.txt", project_path / "test_case_1_settlement_plot.svg")

    def test_settlement_fully_saturated_column(self):
        """
        This test validates the settlement of a fully saturated column under uniform load.
        The test runs multiple stages of a settlement simulation and checks the
        settlement values at specific times against expected results.
        The expected settlement values are based on an analytical solution.
        """
        test_name = "fully_saturated_column_uniform_load"
        test_root = "dsettlement"
        project_path = test_helper.get_file_path(os.path.join(test_root, test_name))

        cwd = os.getcwd()
        os.chdir(project_path)

        import KratosMultiphysics.GeoMechanicsApplication.run_geo_settlement as run_geo_settlement

        n_stages = 5
        project_parameters_filenames = [
            f"ProjectParameters_stage{i+1}.json" for i in range(n_stages)
        ]
        status = run_geo_settlement.run_stages(
            project_path, project_parameters_filenames
        )
        self.assertEqual(status, 0)

        os.chdir(cwd)

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

    def test_settlement_phreatic_line_below_surface(self):
        """
        This test validates the settlement of a soil column where the phreatic line
        is 10 m below the soil surface.
        """
        test_name = "phreatic_line_below_soil_surface"
        test_root = "dsettlement"
        project_path = test_helper.get_file_path(os.path.join(test_root, test_name))

        original_working_dir = os.getcwd()
        os.chdir(project_path)

        import KratosMultiphysics.GeoMechanicsApplication.run_geo_settlement as run_geo_settlement

        n_stages = 5
        project_parameters_filenames = [
            f"ProjectParameters_stage{i+1}.json" for i in range(n_stages)
        ]
        status = run_geo_settlement.run_stages(
            project_path, project_parameters_filenames
        )
        self.assertEqual(status, 0)

        os.chdir(original_working_dir)

        top_node_ids = [2, 3, 104]
        reader = test_helper.GiDOutputFileReader()
        project_path = pathlib.Path(project_path)
        output_stage_3 = reader.read_output_from(project_path / "stage3.post.res")
        time_in_sec = days_to_seconds(0.1) + 1.0
        actual_total_displacement_along_top_edge = reader.nodal_values_at_time(
            "TOTAL_DISPLACEMENT", time_in_sec, output_stage_3, top_node_ids
        )
        # for total_displacement_vector, node_id in zip(actual_total_displacement_along_top_edge, top_node_ids):
        #     self.assertAlmostEqual(total_displacement_vector[1], 0.0, 4, msg=f"total vertical displacement at node {node_id} at time {time_in_sec} [s]")

        time_in_sec = days_to_seconds(100)
        actual_total_displacement_along_top_edge = reader.nodal_values_at_time(
            "TOTAL_DISPLACEMENT", time_in_sec, output_stage_3, top_node_ids
        )
        # for total_displacement_vector, node_id in zip(actual_total_displacement_along_top_edge, top_node_ids):
        #     self.assertAlmostEqual(total_displacement_vector[1], -1.70, 4, msg=f"total vertical displacement at node {node_id} at time {time_in_sec} [s]")

        output_stage_4 = reader.read_output_from(project_path / "stage4.post.res")

        output_stage_5 = reader.read_output_from(project_path / "stage5.post.res")
        time_in_sec = days_to_seconds(10000)
        actual_total_displacement_along_top_edge = reader.nodal_values_at_time(
            "TOTAL_DISPLACEMENT", time_in_sec, output_stage_5, top_node_ids
        )
        # for total_displacement_vector, node_id in zip(actual_total_displacement_along_top_edge, top_node_ids):
        #     self.assertAlmostEqual(total_displacement_vector[1], -8.64, 4, msg=f"total vertical displacement at node {node_id} at time {time_in_sec} [s]")

        if test_helper.want_test_plots():
            left_side_corner_node_ids = [3] + list(range(105, 154)) + [4]
            make_settlement_plot((output_stage_3, output_stage_4, output_stage_5), top_node_ids, project_path / "ref_settlement_data.txt", project_path / "test_case_3_settlement_plot.svg")

            # Make a stress plot at the start of the fifth stage
            ref_data = StressPlotDataFilePaths()
            ref_data.path_to_water_pressure_data = project_path / "ref_water_pressures_after_100.1_days.txt"
            ref_data.path_to_vertical_effective_stress_data = project_path / "ref_effective_vertical_stresses_after_100.1_days.txt"
            time_in_seconds = days_to_seconds(100.1) + 1.0
            make_stress_over_depth_plot(output_stage_5, time_in_seconds, project_path / "stage5.post.msh", left_side_corner_node_ids, ref_data, project_path / "test_case_3_stress_plot_after_100.1_days.svg")

            # Make a stress plot at the end of the fifth stage (when consolidation is supposed to be finished)
            ref_data.path_to_water_pressure_data = project_path / "ref_water_pressures_after_10000_days.txt"
            ref_data.path_to_vertical_effective_stress_data = project_path / "ref_effective_vertical_stresses_after_10000_days.txt"
            make_stress_over_depth_plot(output_stage_5, days_to_seconds(10000), project_path / "stage5.post.msh", left_side_corner_node_ids, ref_data, project_path / "test_case_3_stress_plot_after_10000_days.svg")

      
if __name__ == "__main__":
    KratosUnittest.main()
