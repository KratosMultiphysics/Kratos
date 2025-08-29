from KratosMultiphysics.GeoMechanicsApplication import run_multiple_stages
import KratosMultiphysics.KratosUnittest as KratosUnittest
import os
import matplotlib.pyplot as plt

import test_helper

class PlotDataSeries():
    def __init__(self, x_values, y_values, label='', linestyle='-', marker=None):
        self.x_values = x_values
        self.y_values = y_values
        self.label = label
        self.linestyle = linestyle
        self.marker = marker

def sec_to_day(time_in_sec):
    return time_in_sec / (24.0 * 60.0 * 60.0)


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

        time_in_days = sec_to_day(item["time"])

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

    #plt.show()
    png_file_path = figure_filename
    print(f"Save plot to {png_file_path}")
    plt.savefig(png_file_path)

def plot_water_pressure_results(series_collection, figure_filename):
    figure, axes = plt.subplots(layout="constrained")
    for series in series_collection:
        axes.plot(series.x_values, series.y_values, label=series.label, linestyle=series.linestyle, marker=series.marker)
    axes.grid()
    axes.grid(which="minor", color="0.9")
    axes.set_xlabel('Water pressure [kPa]')
    axes.set_ylabel('Depth [m]')
    figure.legend(loc='outside center right')

    #plt.show()
    png_file_path = figure_filename
    print(f"Save plot to {png_file_path}")
    plt.savefig(png_file_path)


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

        reader = test_helper.GiDOutputFileReader()

        output_data = reader.read_output_from(
            os.path.join(project_path, "stage3.post.res")
        )
        top_middle_node_id = 104
        actual_settlement_after_one_hundred_days = reader.nodal_values_at_time(
            "TOTAL_DISPLACEMENT", 8640000, output_data, [top_middle_node_id]
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

        output_data = reader.read_output_from(
            os.path.join(project_path, "stage5.post.res")
        )
        actual_settlement_after_ten_thousand_days = reader.nodal_values_at_time(
            "TOTAL_DISPLACEMENT", 864000000, output_data, [top_middle_node_id]
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

        os.chdir(original_working_dir)

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

        reader = test_helper.GiDOutputFileReader()

        output_data = reader.read_output_from(
            os.path.join(project_path, "stage2.post.res")
        )
        actual_settlement_after_one_hundred_days = reader.nodal_values_at_time(
            "TOTAL_DISPLACEMENT", 8640000, output_data, [104]
        )[0][1]
        self.assertAlmostEqual(actual_settlement_after_one_hundred_days, -1.71094, 4)

        output_data = reader.read_output_from(
            os.path.join(project_path, "stage5.post.res")
        )
        actual_settlement_after_ten_thousand_days = reader.nodal_values_at_time(
            "TOTAL_DISPLACEMENT", 864000000, output_data, [104]
        )[0][1]
        self.assertAlmostEqual(actual_settlement_after_ten_thousand_days, -8.63753, 4)

        os.chdir(cwd)

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
        output_data = reader.read_output_from(os.path.join(project_path, "stage3.post.res"))
        points_for_node_2 = extract_nodal_settlement_over_time(output_data, 2)
        points_for_node_3 = extract_nodal_settlement_over_time(output_data, 3)
        points_for_node_104 = extract_nodal_settlement_over_time(output_data, 104)

        time_in_sec = 0.1 * 86400.0 + 1.0
        actual_total_displacement_along_top_edge = reader.nodal_values_at_time(
            "TOTAL_DISPLACEMENT", time_in_sec, output_data, top_node_ids
        )
        # for total_displacement_vector, node_id in zip(actual_total_displacement_along_top_edge, top_node_ids):
        #     self.assertAlmostEqual(total_displacement_vector[1], 0.0, 4, msg=f"total vertical displacement at node {node_id} at time {time_in_sec} [s]")

        time_in_sec = 100.0 * 86400.0
        actual_total_displacement_along_top_edge = reader.nodal_values_at_time(
            "TOTAL_DISPLACEMENT", time_in_sec, output_data, top_node_ids
        )
        # for total_displacement_vector, node_id in zip(actual_total_displacement_along_top_edge, top_node_ids):
        #     self.assertAlmostEqual(total_displacement_vector[1], -1.70, 4, msg=f"total vertical displacement at node {node_id} at time {time_in_sec} [s]")

        output_data = reader.read_output_from(os.path.join(project_path, "stage5.post.res"))
        points_for_node_2.extend(extract_nodal_settlement_over_time(output_data, 2))
        points_for_node_3.extend(extract_nodal_settlement_over_time(output_data, 3))
        points_for_node_104.extend(extract_nodal_settlement_over_time(output_data, 104))

        time_in_sec = 10000.0 * 86400.0
        actual_total_displacement_along_top_edge = reader.nodal_values_at_time(
            "TOTAL_DISPLACEMENT", time_in_sec, output_data, top_node_ids
        )
        # for total_displacement_vector, node_id in zip(actual_total_displacement_along_top_edge, top_node_ids):
        #     self.assertAlmostEqual(total_displacement_vector[1], -8.64, 4, msg=f"total vertical displacement at node {node_id} at time {time_in_sec} [s]")

        # Provided by Wijtze Pieter (time in days, settlement in meters)
        reference_points = [
(0.10,	0.00),
(0.10,	0.00),
(0.20,	0.00),
(0.33,	0.00),
(0.49,	0.00),
(0.69,	0.01),
(0.94,	0.01),
(1.26,	0.02),
(1.67,	0.03),
(2.18,	0.05),
(2.82,	0.08),
(3.63,	0.11),
(4.65,	0.14),
(5.94,	0.19),
(7.56,	0.25),
(9.60,	0.31),
(12.18,	0.39),
(15.43,	0.48),
(19.53,	0.59),
(24.70,	0.71),
(31.21,	0.84),
(39.42,	0.99),
(49.77,	1.15),
(62.82,	1.32),
(79.26,	1.51),
(100.00,	1.70),
(100.10,	2.01),
(100.20,	2.13),
(100.33,	2.22),
(100.49,	2.31),
(100.69,	2.40),
(100.94,	2.48),
(101.00,	2.50),
(101.26,	2.56),
(101.66,	2.63),
(102.16,	2.70),
(102.79,	2.76),
(103.59,	2.82),
(104.60,	2.88),
(105.87,	2.94),
(107.47,	3.01),
(109.48,	3.09),
(112.01,	3.17),
(115.20,	3.25),
(119.22,	3.35),
(124.29,	3.46),
(130.66,	3.57),
(138.69,	3.70),
(148.81,	3.84),
(161.55,	3.99),
(177.60,	4.14),
(197.81,	4.31),
(223.26,	4.49),
(255.32,	4.68),
(295.71,	4.87),
(346.57,	5.07),
(410.63,	5.28),
(491.31,	5.49),
(592.93,	5.70),
(720.93,	5.92),
(882.14,	6.14),
(1085.18,	6.36),
(1340.92,	6.59),
(1663.02,	6.81),
(2068.71,	7.04),
(2579.67,	7.27),
(3223.24,	7.49),
(4033.82,	7.72),
(5054.74,	7.95),
(6340.61,	8.18),
(7960.16,	8.41),
(10000.00,	8.64),
            ]

        graph_series = []
        graph_series.append(make_plot_data(reference_points, 'ref', marker='+'))
        graph_series.append(make_plot_data(points_for_node_2, 'node 2', linestyle='-', marker='+'))
        graph_series.append(make_plot_data(points_for_node_3, 'node 3', linestyle='--', marker='+'))
        graph_series.append(make_plot_data(points_for_node_104, 'node 104', linestyle=':', marker='+'))

        plot_settlement_results(graph_series, "test_case_3_settlement_plot.png")

        # water pressure values from "old" D-Settlement [in kPa] at t = 10,000 days
        reference_points = [
(0.00,	0.00),
(-0.10,	0.00),
(-0.20,	0.00),
(-0.30,	0.00),
(-0.40,	0.00),
(-0.50,	0.00),
(-0.60,	0.00),
(-0.70,	0.00),
(-0.80,	0.00),
(-0.90,	0.00),
(-1.00,	0.00),
(-1.50,	0.00),
(-2.00,	0.00),
(-2.50,	0.00),
(-3.00,	0.00),
(-3.50,	0.00),
(-4.00,	0.00),
(-4.50,	0.00),
(-5.00,	0.00),
(-5.50,	0.00),
(-6.00,	0.00),
(-6.50,	0.00),
(-7.00,	0.00),
(-7.50,	0.00),
(-8.00,	0.00),
(-8.50,	0.00),
(-9.00,	0.00),
(-9.50,	0.00),
(-10.00,	0.00),
(-10.57,	5.67),
(-11.13,	11.34),
(-11.70,	17.01),
(-12.20,	22.01),
(-12.70,	27.01),
(-13.20,	32.01),
(-13.70,	37.01),
(-14.20,	42.01),
(-14.70,	47.01),
(-15.20,	52.01),
(-15.70,	57.01),
(-16.20,	62.01),
(-16.70,	67.01),
(-17.20,	72.01),
(-17.70,	77.01),
(-18.20,	82.01),
(-18.70,	87.01),
(-19.20,	92.01),
(-19.70,	97.01),
(-20.20,	102.01),
(-20.70,	107.01),
(-21.20,	112.01),
(-21.70,	117.01),
(-22.20,	122.01),
(-22.70,	127.01),
(-23.20,	132.01),
(-23.70,	137.01),
(-24.13,	141.34),
(-24.57,	145.67),
(-25.00,	150.01),
(-25.50,	155.01),
(-26.00,	160.01),
(-26.50,	165.01),
(-27.00,	170.01),
(-27.50,	175.01),
(-28.00,	180.01),
(-28.50,	185.01),
(-29.00,	190.01),
(-29.50,	195.01),
(-30.00,	200.01),
(-30.53,	205.34),
(-31.07,	210.67),
(-31.60,	216.01),
(-32.10,	221.01),
(-32.60,	226.01),
(-33.10,	231.01),
(-33.60,	236.01),
(-34.07,	240.67),
(-34.53,	245.34),
(-35.00,	250.01),
(-35.53,	255.34),
(-36.07,	260.67),
(-36.60,	266.01),
(-37.10,	271.01),
(-37.60,	276.01),
(-38.10,	281.01),
(-38.60,	286.00),
(-39.07,	290.67),
(-39.53,	295.34),
(-40.00,	300.00),
(-40.50,	305.00),
(-41.00,	310.00),
(-41.50,	315.00),
(-42.00,	320.00),
(-42.50,	325.00),
(-43.00,	330.00),
(-43.50,	335.00),
(-44.00,	340.00),
(-44.50,	345.00),
(-45.00,	350.00),
(-45.50,	355.00),
(-46.00,	360.00),
(-46.50,	365.00),
(-47.00,	370.00),
(-47.50,	375.00),
(-48.00,	380.00),
(-48.50,	385.00),
(-49.00,	390.00),
(-49.50,	395.00),
(-50.00,	400.00),
        ]

        graph_series = []
        graph_series.append(make_water_pressure_plot_data(reference_points, 'ref P_w', marker='+'))

        # Effective vertical stresses from "old" D-Settlement [in kPa] at t = 10,000 days
        reference_points = [
(0.00,	20.00),
(-0.10,	20.60),
(-0.20,	21.20),
(-0.30,	21.80),
(-0.40,	22.40),
(-0.50,	23.00),
(-0.60,	23.60),
(-0.70,	24.20),
(-0.80,	24.80),
(-0.90,	25.40),
(-1.00,	26.00),
(-1.50,	29.00),
(-2.00,	32.00),
(-2.50,	35.00),
(-3.00,	38.00),
(-3.50,	41.00),
(-4.00,	44.00),
(-4.50,	47.00),
(-5.00,	50.00),
(-5.50,	53.00),
(-6.00,	56.00),
(-6.50,	59.00),
(-7.00,	62.00),
(-7.50,	65.00),
(-8.00,	68.00),
(-8.50,	71.00),
(-9.00,	74.00),
(-9.50,	77.00),
(-10.00,	80.00),
(-10.57,	83.40),
(-11.13,	86.80),
(-11.70,	90.19),
(-12.20,	93.19),
(-12.70,	96.19),
(-13.20,	99.19),
(-13.70,	102.19),
(-14.20,	105.19),
(-14.70,	108.19),
(-15.20,	111.19),
(-15.70,	114.19),
(-16.20,	117.19),
(-16.70,	120.19),
(-17.20,	123.19),
(-17.70,	126.19),
(-18.20,	129.19),
(-18.70,	132.19),
(-19.20,	135.19),
(-19.70,	138.19),
(-20.20,	141.19),
(-20.70,	144.19),
(-21.20,	147.19),
(-21.70,	150.19),
(-22.20,	153.19),
(-22.70,	156.19),
(-23.20,	159.19),
(-23.70,	162.19),
(-24.13,	164.79),
(-24.57,	167.39),
(-25.00,	169.99),
(-25.50,	172.99),
(-26.00,	175.99),
(-26.50,	178.99),
(-27.00,	181.99),
(-27.50,	184.99),
(-28.00,	187.99),
(-28.50,	190.99),
(-29.00,	193.99),
(-29.50,	196.99),
(-30.00,	199.99),
(-30.53,	203.19),
(-31.07,	206.39),
(-31.60,	209.59),
(-32.10,	212.59),
(-32.60,	215.59),
(-33.10,	218.59),
(-33.60,	221.59),
(-34.07,	224.39),
(-34.53,	227.19),
(-35.00,	229.99),
(-35.53,	233.19),
(-36.07,	236.39),
(-36.60,	239.59),
(-37.10,	242.59),
(-37.60,	245.59),
(-38.10,	248.59),
(-38.60,	251.60),
(-39.07,	254.40),
(-39.53,	257.20),
(-40.00,	260.00),
(-40.50,	263.00),
(-41.00,	266.00),
(-41.50,	269.00),
(-42.00,	272.00),
(-42.50,	275.00),
(-43.00,	278.00),
(-43.50,	281.00),
(-44.00,	284.00),
(-44.50,	287.00),
(-45.00,	290.00),
(-45.50,	293.00),
(-46.00,	296.00),
(-46.50,	299.00),
(-47.00,	302.00),
(-47.50,	305.00),
(-48.00,	308.00),
(-48.50,	311.00),
(-49.00,	314.00),
(-49.50,	317.00),
(-50.00,	320.00),
        ]

        graph_series.append(make_water_pressure_plot_data(reference_points, 'ref eff sigma_yy', marker='+'))

        plot_water_pressure_results(graph_series, "test_case_3_stress_plot_after_10000_days.png")

      
if __name__ == "__main__":
    KratosUnittest.main()
