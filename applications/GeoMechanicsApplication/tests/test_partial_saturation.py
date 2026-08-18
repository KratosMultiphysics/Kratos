import csv
import os
from dataclasses import dataclass

import KratosMultiphysics as Kratos
import KratosMultiphysics.GeoMechanicsApplication.geomechanics_analysis as analysis
import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_helper
from KratosMultiphysics.GeoMechanicsApplication.gid_output_file_reader import (
    GiDOutputFileReader,
)
from KratosMultiphysics.GeoMechanicsApplication.unit_conversions import Pa_to_kPa

if test_helper.want_test_plots():
    import KratosMultiphysics.GeoMechanicsApplication.geo_plot_utilities as plot_utils


class KratosGeoMechanicsPartialSaturation(KratosUnittest.TestCase):
    """
    This class contains benchmark tests which are checked with the analytical solution
    """

    def __test_saturated_below_phreatic_level_pw(self, test_name):
        n_stages = 2
        # get the parameter file names for all stages
        file_path = test_helper.get_file_path(
            os.path.join("test_partially_saturated", test_name)
        )
        parameter_file_names = [
            os.path.join(file_path, f"ProjectParameters_stage{i+1}.json")
            for i in range(n_stages)
        ]

        # set stage parameters
        parameters_stages = []
        initial_directory = os.getcwd()
        os.chdir(file_path)
        for parameter_file_name in parameter_file_names:
            with open(parameter_file_name, "r") as parameter_file:
                parameters_stages.append(Kratos.Parameters(parameter_file.read()))

        model = Kratos.Model()

        # run stages and get water pressure/displacement results per stage
        stage_water_pressure = []
        coords = []
        for stage_parameters in parameters_stages:
            stage = analysis.GeoMechanicsAnalysis(model, stage_parameters)
            stage.Run()
            stage_water_pressure.append(test_helper.get_water_pressure(stage))
            coords.append(test_helper.get_nodal_coordinates(stage))

        os.chdir(initial_directory)
        # get y coords of all the nodes
        y_coords = [coord[1] for coord in coords[0]]

        # calculate water pressure analytical solution for all stages and calculate the error
        rel_p_stage = [
            self.__compute_hydrostatic_water_pressure(y_coord, -2.0)
            for y_coord in y_coords
        ]

        errors_stage = [
            actual_pressure - expected_pressure
            for actual_pressure, expected_pressure in zip(
                stage_water_pressure[1], rel_p_stage
            )
        ]
        rmse_stages = (
            sum([error**2 for error in errors_stage]) / len(errors_stage)
        ) ** 0.5

        # assert if average error in all stages is below accuracy
        accuracy = 1.0e-3
        self.assertLess(rmse_stages, accuracy)

    # def test_saturated_below_phreatic_level_pw_triangle3N(self):
    #     self.__test_saturated_below_phreatic_level_pw(
    #         "test_saturated_below_phreatic_level_pw_triangle3N"
    #     )
    #
    # def test_saturated_below_phreatic_level_pw_triangle6N(self):
    #     self.__test_saturated_below_phreatic_level_pw(
    #         "test_saturated_below_phreatic_level_pw_triangle6N"
    #     )
    #
    # def test_saturated_below_phreatic_level_upw_difforder_triangle6n(self):
    #     self.__test_saturated_below_phreatic_level_pw(
    #         "test_saturated_below_phreatic_level_upw_difforder_triangle6n"
    #     )
    #
    # def test_saturated_below_phreatic_level_upw_smallstrain_triangle3n(self):
    #     self.__test_saturated_below_phreatic_level_pw(
    #         "test_saturated_below_phreatic_level_upw_smallstrain_triangle3n"
    #     )
    #
    # def test_saturated_below_phreatic_level_upw_smallstrain_triangle6n(self):
    #     self.__test_saturated_below_phreatic_level_pw(
    #         "test_saturated_below_phreatic_level_upw_smallstrain_triangle6n"
    #     )

    # def test_climbing_falling_phreatic_level_upw_smallstrain_quad4n(self):
    #     # only waterpressures below phreatic level are checked with an analytical solution.
    #     # values above phreatic level give suction of an unchecked amount.
    #     file_path = test_helper.get_file_path(
    #         os.path.join(
    #             "test_partially_saturated",
    #             "test_rising_falling_phreatic_level_pw_quad4N",
    #         )
    #     )
    #     simulation = test_helper.run_kratos(file_path)
    #
    #     reader = GiDOutputFileReader()
    #     output_data = reader.read_output_from(
    #         os.path.join(file_path, "rising_falling_phreatic_level_pw_quad4n.post.res")
    #     )
    #     coords = test_helper.get_nodal_coordinates(simulation)
    #     times = [1.0, 5.0, 9.0, 13.0, 17.0, 21.0, 25.0, 29.0]
    #     water_levels = [-4.0, -3.0, -2.0, -1.0, -2.0, -3.0, -4.0, -5.0]
    #     for time, water_level in zip(times, water_levels):
    #         water_pressures = reader.nodal_values_at_time(
    #             "WATER_PRESSURE", time, output_data
    #         )
    #         negative_water_pressures = [
    #             min([water_pressure, 0.0]) for water_pressure in water_pressures
    #         ]
    #         analytical_water_pressures = [
    #             self.__compute_hydrostatic_water_pressure(coord[1], water_level)
    #             for coord in coords
    #         ]
    #         self.assertVectorAlmostEqual(
    #             negative_water_pressures,
    #             analytical_water_pressures,
    #             places=None,
    #             msg=f"water pressures at time {time}",
    #             delta=10.0,
    #         )

    def __compute_hydrostatic_water_pressure(self, y_coord, phreatic_level):
        water_weight = -10000.0
        result = water_weight * (phreatic_level - y_coord)
        return min([result, 0.0])

    def test_infiltration_from_top_boundary(self):
        file_path = test_helper.get_file_path(
            os.path.join(
                "test_partially_saturated", "infiltration", "test_infiltration_pw"
            )
        )
        simulation = test_helper.run_kratos(file_path)

        reader = GiDOutputFileReader()
        output_data = reader.read_output_from(
            os.path.join(file_path, "run1sim5_map_hydro.post.res")
        )

        depth_by_id_for_left_boundary_nodes = self._calculate_depth_for_boundary_nodes(
            simulation
        )

        @dataclass
        class ExpectedResult:
            node_id: int
            value: float

        expected_results_at_times = {
            12000.0: [
                ExpectedResult(node_id=1, value=0.0),
                ExpectedResult(node_id=17, value=6243.59),
                ExpectedResult(node_id=26, value=16400.0),
            ],
            72000.0: [
                ExpectedResult(node_id=55, value=0.0),
                ExpectedResult(node_id=61, value=5013.57),
                ExpectedResult(node_id=70, value=7535.47),
            ],
            96000.0: [
                ExpectedResult(node_id=87, value=-293.249),
                ExpectedResult(node_id=91, value=970.378),
                ExpectedResult(node_id=100, value=1672.53),
            ],
            192000.0: [
                # This is the hydrostatic line from p = 0 at depth = 0m
                # and p = -20000 Pa at depth = 2m, as seen in the plot
                ExpectedResult(node_id=1, value=0.0),
                ExpectedResult(node_id=58, value=-10000.0),
                ExpectedResult(node_id=3, value=-20000.0),
            ],
        }

        if test_helper.want_test_plots():
            self.create_pressure_depth_plots(
                depth_by_id_for_left_boundary_nodes,
                expected_results_at_times,
                file_path,
                output_data,
                [12000.0, 24000.0, 36000.0, 48000.0, 72000.0, 96000.0, 192000.0],
            )

        variable_name = "WATER_PRESSURE"
        self._validate_outputs_against_expected_results(
            reader, output_data, expected_results_at_times, variable_name
        )

    def _validate_outputs_against_expected_results(
        self, reader, output_data, expected_results_at_times, variable_name
    ):
        for time, expected_results in expected_results_at_times.items():
            water_pressures = reader.nodal_values_at_time(
                variable_name,
                time,
                output_data,
                [result.node_id for result in expected_results],
            )
            expected_values = [result.value for result in expected_results]
            self.assertVectorAlmostEqual(
                water_pressures, expected_values, places=None, delta=10.0
            )

    # def test_infiltration_from_top_boundary_O13(self):
    #     file_path = test_helper.get_file_path(
    #         os.path.join("test_partially_saturated", "test_infiltration_pw_caseO13")
    #     )
    #     simulation = test_helper.run_kratos(file_path)
    #
    #     reader = GiDOutputFileReader()
    #     output_data = reader.read_output_from(
    #         os.path.join(file_path, "run1sim5_map_hydro.post.res")
    #     )
    #
    #     depth_by_id_for_left_boundary_nodes = self._calculate_depth_for_boundary_nodes(
    #         simulation
    #     )
    #
    #     @dataclass
    #     class ExpectedResult:
    #         node_id: int
    #         value: float
    #
    #     expected_results_at_times = {
    #         60: [],
    #         3600.0: [],
    #         7200.0: [],
    #         10800.0: [],
    #         14400.0: [],
    #     }
    #
    #     if test_helper.want_test_plots():
    #         self.create_pressure_depth_plots(
    #             depth_by_id_for_left_boundary_nodes,
    #             expected_results_at_times,
    #             file_path,
    #             output_data,
    #             plot_times=[5220.0, 10440.0, 14400.0],
    #         )
    #         self.create_saturation_depth_plots(
    #             depth_by_id_for_left_boundary_nodes,
    #             expected_results_at_times,
    #             file_path,
    #             output_data,
    #             plot_times=expected_results_at_times.keys(),
    #         )

    def test_infiltration_from_top_boundary_O06(self):
        file_path = test_helper.get_file_path(
            os.path.join(
                "test_partially_saturated",
                "infiltration",
                "test_infiltration_pw_caseO06",
            )
        )
        simulation = test_helper.run_kratos(file_path)

        reader = GiDOutputFileReader()
        output_data = reader.read_output_from(
            os.path.join(file_path, "run1sim5_map_hydro.post.res")
        )

        depth_by_id_for_left_boundary_nodes = self._calculate_depth_for_boundary_nodes(
            simulation
        )

        @dataclass
        class ExpectedResult:
            node_id: int
            value: float

        expected_results_at_times = {
            60.0: [],
            3600.0: [],
            7200.0: [],
            10800.0: [],
            14400.0: [],
        }

        if test_helper.want_test_plots():
            self.create_pressure_depth_plots(
                depth_by_id_for_left_boundary_nodes,
                expected_results_at_times,
                file_path,
                output_data,
                plot_times=expected_results_at_times.keys(),
                test_name="O6",
                time_strings=["3622", "7235", "10,81*10\^3", "14,40*10\^3"],
            )
            self.create_saturation_depth_plots(
                depth_by_id_for_left_boundary_nodes,
                expected_results_at_times,
                file_path,
                output_data,
                plot_times=expected_results_at_times.keys(),
                test_name="O6",
            )

    def test_infiltration_from_top_boundary_O10(self):
        file_path = test_helper.get_file_path(
            os.path.join(
                "test_partially_saturated",
                "infiltration",
                "test_infiltration_pw_caseO10",
            )
        )
        simulation = test_helper.run_kratos(file_path)

        reader = GiDOutputFileReader()
        output_data = reader.read_output_from(
            os.path.join(file_path, "run1sim5_map_hydro.post.res")
        )

        depth_by_id_for_left_boundary_nodes = self._calculate_depth_for_boundary_nodes(
            simulation
        )

        @dataclass
        class ExpectedResult:
            node_id: int
            value: float

        expected_results_at_times = {
            60.0: [],
            3600.0: [],
            7200.0: [],
            10800.0: [],
            14400.0: [],
        }

        if test_helper.want_test_plots():
            self.create_pressure_depth_plots(
                depth_by_id_for_left_boundary_nodes,
                expected_results_at_times,
                file_path,
                output_data,
                plot_times=expected_results_at_times.keys(),
                test_name="O10",
                time_strings=["3604", "7218", "10,83*10\^3", "14,40*10\^3"],
            )
            self.create_saturation_depth_plots(
                depth_by_id_for_left_boundary_nodes,
                expected_results_at_times,
                file_path,
                output_data,
                plot_times=expected_results_at_times.keys(),
                test_name="O10",
            )
    def test_infiltration_from_top_boundary_B10(self):
        file_path = test_helper.get_file_path(
            os.path.join(
                "test_partially_saturated",
                "infiltration",
                "test_infiltration_pw_caseB10",
            )
        )
        simulation = test_helper.run_kratos(file_path)

        reader = GiDOutputFileReader()
        output_data = reader.read_output_from(
            os.path.join(file_path, "run1sim5_map_hydro.post.res")
        )

        depth_by_id_for_left_boundary_nodes = self._calculate_depth_for_boundary_nodes(
            simulation
        )

        @dataclass
        class ExpectedResult:
            node_id: int
            value: float

        expected_results_at_times = {
            60.0: [],
            3600.0: [],
            7200.0: [],
            10800.0: [],
            14400.0: [],
        }

        if test_helper.want_test_plots():
            self.create_pressure_depth_plots(
                depth_by_id_for_left_boundary_nodes,
                expected_results_at_times,
                file_path,
                output_data,
                plot_times=expected_results_at_times.keys(),
                test_name="B10",
                time_strings=["3617", "7231", "10,81*10\^3", "14,40*10\^3"],
            )
            self.create_saturation_depth_plots(
                depth_by_id_for_left_boundary_nodes,
                expected_results_at_times,
                file_path,
                output_data,
                plot_times=expected_results_at_times.keys(),
                test_name="B10",
            )
    def test_infiltration_from_top_boundary_B06(self):
        file_path = test_helper.get_file_path(
            os.path.join(
                "test_partially_saturated",
                "infiltration",
                "test_infiltration_pw_caseB06",
            )
        )
        simulation = test_helper.run_kratos(file_path)

        reader = GiDOutputFileReader()
        output_data = reader.read_output_from(
            os.path.join(file_path, "run1sim5_map_hydro.post.res")
        )

        depth_by_id_for_left_boundary_nodes = self._calculate_depth_for_boundary_nodes(
            simulation
        )

        @dataclass
        class ExpectedResult:
            node_id: int
            value: float

        expected_results_at_times = {
            60.0: [],
            3600.0: [],
            7200.0: [],
            10800.0: [],
            14400.0: [],
        }

        if test_helper.want_test_plots():
            self.create_pressure_depth_plots(
                depth_by_id_for_left_boundary_nodes,
                expected_results_at_times,
                file_path,
                output_data,
                plot_times=expected_results_at_times.keys(),
                test_name="B6",
                time_strings=["3608", "7241", "10,82*10\^3", "14,40*10\^3"],
            )
            self.create_saturation_depth_plots(
                depth_by_id_for_left_boundary_nodes,
                expected_results_at_times,
                file_path,
                output_data,
                plot_times=expected_results_at_times.keys(),
            )

    def test_infiltration_from_top_boundary_B4(self):
        file_path = test_helper.get_file_path(
            os.path.join(
                "test_partially_saturated",
                "infiltration",
                "test_infiltration_pw_caseB4",
            )
        )
        simulation = test_helper.run_kratos(file_path)

        reader = GiDOutputFileReader()
        output_data = reader.read_output_from(
            os.path.join(file_path, "run1sim5_map_hydro.post.res")
        )

        depth_by_id_for_left_boundary_nodes = self._calculate_depth_for_boundary_nodes(
            simulation
        )

        @dataclass
        class ExpectedResult:
            node_id: int
            value: float

        expected_results_at_times = {
            60.0: [],
            3600.0: [],
            7200.0: [],
            10800.0: [],
            14400.0: [],
        }

        if test_helper.want_test_plots():
            self.create_pressure_depth_plots(
                depth_by_id_for_left_boundary_nodes,
                expected_results_at_times,
                file_path,
                output_data,
                plot_times=expected_results_at_times.keys(),
                test_name="B4",
                time_strings=["3640", "7338", "10,90*10\^3", "14,40*10\^3"],
            )
            self.create_saturation_depth_plots(
                depth_by_id_for_left_boundary_nodes,
                expected_results_at_times,
                file_path,
                output_data,
                plot_times=expected_results_at_times.keys(),
            )

    def _calculate_depth_for_boundary_nodes(self, simulation):
        return {
            node.Id: -1.0 * node.Y
            for node in simulation.model.GetModelPart(
                "PorousDomain.porous_computational_model_part"
            ).Nodes
            if node.X == 0.0
        }

    def test_infiltration_from_top_boundary_B09(self):
        file_path = test_helper.get_file_path(
            os.path.join(
                "test_partially_saturated",
                "infiltration",
                "test_infiltration_pw_caseB09",
            )
        )
        simulation = test_helper.run_kratos(file_path)

        output_data = GiDOutputFileReader().read_output_from(
            os.path.join(file_path, "run1sim5_map_hydro.post.res")
        )

        depth_by_id_for_left_boundary_nodes = self._calculate_depth_for_boundary_nodes(
            simulation
        )

        @dataclass
        class ExpectedResult:
            node_id: int
            value: float

        expected_results_at_times = {
            150.0: [],
            3600.0: [],
            7200.0: [],
            10800.0: [],
            14400.0: [],
        }

        if test_helper.want_test_plots():
            self.create_pressure_depth_plots(
                depth_by_id_for_left_boundary_nodes,
                expected_results_at_times,
                file_path,
                output_data,
                plot_times=expected_results_at_times.keys(),
                test_name="B9",
                time_strings=["3626", "7354", "10,86*10\^3", "14,40*10\^3"],
            )
            self.create_saturation_depth_plots(
                depth_by_id_for_left_boundary_nodes,
                expected_results_at_times,
                file_path,
                output_data,
                plot_times=expected_results_at_times.keys(),
            )

    def create_pressure_depth_plots(
        self,
        depth_by_id_for_left_boundary_nodes,
        expected_results_at_times,
        file_path,
        output_data,
        plot_times,
        test_name=None,
        time_strings=None,
    ):
        data_series_collection = []
        colors = ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple", "tab:brown", "tab:pink", "tab:gray"]
        for time, color in zip(plot_times, colors):
            water_pressures = GiDOutputFileReader.nodal_values_at_time(
                "WATER_PRESSURE",
                time,
                output_data,
                depth_by_id_for_left_boundary_nodes.keys(),
            )
            water_pressures = [Pa_to_kPa(pressure) for pressure in water_pressures]
            sorted_depth, sorted_pressures = zip(
                *sorted(
                    zip(
                        depth_by_id_for_left_boundary_nodes.values(),
                        water_pressures,
                    )
                )
            )
            data = zip(sorted_pressures, sorted_depth)
            data_series_collection.append(
                plot_utils.DataSeries(
                    data,
                    label=f"Time = {int(time)}s",
                    line_style="-",
                    marker="",
                    color=color,
                )
            )

        file_name = "expected_water_pressures_internal_reference.csv"
        expected_pressures_file = os.path.join(file_path, file_name)
        if os.path.exists(expected_pressures_file):
            for time, color in zip([0.0, 3600.0, 7200.0, 10800.0, 14400.0], colors):
                with open(
                    expected_pressures_file,
                    newline="",
                ) as csv_file:
                    reader = csv.DictReader(csv_file, skipinitialspace=True)
                    data_points = [
                        (-1.0 * float(row[f"t={int(time)}"]), -1.0 * float(row["d"]))
                        for row in reader
                    ]
                    data_series_collection.append(
                        plot_utils.DataSeries(
                            data_points,
                            label=f"DG-Flow Reference t = {time}",
                            line_style="--",
                        )
                    )

        file_name = "expected_water_pressures_hydrus_reference.csv"
        expected_pressures_file = os.path.join(file_path, file_name)
        if os.path.exists(expected_pressures_file):
            for time, color in zip([0.0, 3600.0, 7200.0, 10800.0, 14400.0], colors):
                with open(
                    expected_pressures_file,
                    newline="",
                ) as csv_file:
                    reader = csv.DictReader(csv_file, skipinitialspace=True)
                    data_points = [
                        (-1.0 * float(row[f"t={int(time)}"]), float(row["d"]))
                        for row in reader
                    ]
                    data_series_collection.append(
                        plot_utils.DataSeries(
                            data_points,
                            label=f"Hydrus Reference t = {time}",
                            line_style=":",
                            color=color,
                        )
                    )

        asserted_data_points = []
        for time, expected_results in expected_results_at_times.items():
            for expected_result in expected_results:
                water_pressure = Pa_to_kPa(expected_result.value)

                asserted_data_points.append(
                    (
                        water_pressure,
                        depth_by_id_for_left_boundary_nodes[expected_result.node_id],
                    )
                )
        data_series_collection.append(
            plot_utils.DataSeries(
                asserted_data_points,
                label=f"Asserted pressures",
                line_style="",
                marker="x",
                color="r",
            )
        )
        expected_water_pressures = os.path.join(
            file_path, "expected_water_pressures.csv"
        )

        if os.path.exists(expected_water_pressures) and test_name:
            for time, color in zip(time_strings, colors[1:]):
                with open(
                    expected_water_pressures,
                    newline="",
                ) as csv_file:
                    reader = csv.DictReader(csv_file, skipinitialspace=True)
                    data_points = [
                        (
                            float(row[f"{test_name}_{time} s "]),
                            -1.0 * float(row["Y coordinate [m]"]),
                        )
                        for row in reader
                    ]
                    data_series_collection.append(
                        plot_utils.DataSeries(
                            data_points,
                            label=f"External Reference t = {time}",
                            line_style="-.",
                            color=color,
                        )
                    )
        plot_utils._make_plot(
            data_series_collection,
            os.path.join(file_path, "infiltration_from_top_boundary.svg"),
            xlabel="water pressure [kPa]",
            ylabel="depth [m]",
            yaxis_inverted=True,
        )

    def create_saturation_depth_plots(
        self,
        depth_by_id_for_left_boundary_nodes,
        expected_results_at_times,
        file_path,
        output_data,
        plot_times,
        test_name=None,
    ):
        data_series_collection = []
        colors = ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple"]

        for time, color in zip(plot_times, colors):
            saturations = GiDOutputFileReader.nodal_values_at_time(
                "EFFECTIVE_SATURATION",
                time,
                output_data,
                depth_by_id_for_left_boundary_nodes.keys(),
            )
            saturations = [saturation * 100 for saturation in saturations]
            sorted_depth, sorted_pressures = zip(
                *sorted(
                    zip(
                        depth_by_id_for_left_boundary_nodes.values(),
                        saturations,
                    )
                )
            )
            data = zip(sorted_pressures, sorted_depth)
            data_series_collection.append(
                plot_utils.DataSeries(
                    data, label=f"Time = {int(time)}s", line_style="-", marker=""
                )
            )
        asserted_data_points = []
        for time, expected_results in expected_results_at_times.items():
            for expected_result in expected_results:
                water_pressure = Pa_to_kPa(expected_result.value)

                asserted_data_points.append(
                    (
                        water_pressure,
                        depth_by_id_for_left_boundary_nodes[expected_result.node_id],
                    )
                )
        data_series_collection.append(
            plot_utils.DataSeries(
                asserted_data_points,
                label=f"Asserted pressures",
                line_style="",
                marker="x",
                color="r",
            )
        )
        expected_saturations_file = os.path.join(file_path, "expected_saturations.csv")

        if os.path.exists(expected_saturations_file) and test_name:
            for time, color in zip([3782.0, 7382.0, 10982.0, 14400.0], colors[1:]):
                with open(
                    expected_saturations_file,
                    newline="",
                ) as csv_file:
                    reader = csv.DictReader(csv_file, skipinitialspace=True)
                    data_points = [
                        (
                            float(row[f"{test_name}_{int(time)}s"]),
                            -1.0 * float(row["Y"]),
                        )
                        for row in reader
                    ]
                    data_series_collection.append(
                        plot_utils.DataSeries(
                            data_points,
                            label=f"External Reference t = {time}",
                            line_style="--",
                            color=color,
                        )
                    )
        plot_utils._make_plot(
            data_series_collection,
            os.path.join(file_path, "saturation_depth_plots.svg"),
            xlabel="saturation [%]",
            ylabel="depth [m]",
            yaxis_inverted=True,
        )


if __name__ == "__main__":
    KratosUnittest.main()
