import os

import KratosMultiphysics as Kratos
import KratosMultiphysics.GeoMechanicsApplication.geomechanics_analysis as analysis
from KratosMultiphysics.GeoMechanicsApplication.gid_output_file_reader import (
    GiDOutputFileReader,
)
from dataclasses import dataclass
import KratosMultiphysics.KratosUnittest as KratosUnittest

import test_helper

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

    def test_saturated_below_phreatic_level_pw_triangle3N(self):
        self.__test_saturated_below_phreatic_level_pw('test_saturated_below_phreatic_level_pw_triangle3N')

    def test_saturated_below_phreatic_level_pw_triangle6N(self):
        self.__test_saturated_below_phreatic_level_pw('test_saturated_below_phreatic_level_pw_triangle6N')

    def test_saturated_below_phreatic_level_upw_difforder_triangle6n(self):
        self.__test_saturated_below_phreatic_level_pw('test_saturated_below_phreatic_level_upw_difforder_triangle6n')

    def test_saturated_below_phreatic_level_upw_smallstrain_triangle3n(self):
        self.__test_saturated_below_phreatic_level_pw('test_saturated_below_phreatic_level_upw_smallstrain_triangle3n')

    def test_saturated_below_phreatic_level_upw_smallstrain_triangle6n(self):
        self.__test_saturated_below_phreatic_level_pw('test_saturated_below_phreatic_level_upw_smallstrain_triangle6n')

    def test_climbing_falling_phreatic_level_upw_smallstrain_quad4n(self):
        # only waterpressures below phreatic level are checked with an analytical solution.
        # values above phreatic level give suction of an unchecked amount.
        file_path = test_helper.get_file_path(os.path.join('test_partially_saturated', 'test_rising_falling_phreatic_level_pw_quad4N'))
        simulation = test_helper.run_kratos(file_path)

        reader = GiDOutputFileReader()
        output_data = reader.read_output_from(os.path.join(file_path, 'rising_falling_phreatic_level_pw_quad4n.post.res'))
        coords = test_helper.get_nodal_coordinates(simulation)
        times = [1.0, 5.0, 9.0, 13.0, 17.0, 21.0, 25.0, 29.0]
        water_levels = [-4.0, -3.0, -2.0, -1.0, -2.0, -3.0, -4.0, -5.0]
        for time, water_level in zip(times, water_levels):
            water_pressures = reader.nodal_values_at_time('WATER_PRESSURE', time, output_data)
            negative_water_pressures = [min([water_pressure, 0.0]) for water_pressure in water_pressures]
            analytical_water_pressures = [self.__compute_hydrostatic_water_pressure(coord[1], water_level) for coord in coords]
            self.assertVectorAlmostEqual(negative_water_pressures, analytical_water_pressures, places=None, msg=f"water pressures at time {time}", delta=10.0)

    def __compute_hydrostatic_water_pressure(self, y_coord, phreatic_level):
        water_weight = -10000.0
        result = water_weight * (phreatic_level - y_coord)
        return min([result, 0.0])

    def test_infiltration_from_top_boundary(self):
        file_path = test_helper.get_file_path(
            os.path.join("test_partially_saturated", "test_infiltration_pw")
        )
        simulation = test_helper.run_kratos(file_path)

        reader = GiDOutputFileReader()
        output_data = reader.read_output_from(
            os.path.join(file_path, "run1sim5_map_hydro.post.res")
        )

        node_ids_left_boundary = []
        depth_boundary_nodes = []
        for node in simulation.model.GetModelPart(
            "PorousDomain.porous_computational_model_part"
        ).Nodes:
            if node.X == 0.0:
                node_ids_left_boundary.append(node.Id)
                depth_boundary_nodes.append(-1.0 * node.Y)

        @dataclass
        class ExpectedResult:
            node_id: int
            value: float

        expected_results_at_times = {
            12000: [
                ExpectedResult(node_id=1, value=0.0),
                ExpectedResult(node_id=17, value=6243.59),
                ExpectedResult(node_id=26, value=16400),
            ],
            72000: [
                ExpectedResult(node_id=55, value=0.0),
                ExpectedResult(node_id=61, value=5013.57),
                ExpectedResult(node_id=70, value=7535.47),
            ],
            96000: [
                ExpectedResult(node_id=87, value=-293.249),
                ExpectedResult(node_id=91, value=970.378),
                ExpectedResult(node_id=100, value=1672.53),
            ],
            192000: [
                # This is the hydrostatic line from p = 0 at depth = 0m
                # and p = -20000 Pa at depth = 2m, as seen in the plot
                ExpectedResult(node_id=1, value=0.0),
                ExpectedResult(node_id=58, value=-10000),
                ExpectedResult(node_id=3, value=-20000),
            ],
        }

        for time, expected_results in expected_results_at_times.items():
            water_pressures = reader.nodal_values_at_time(
                "WATER_PRESSURE",
                time,
                output_data,
                [result.node_id for result in expected_results],
            )
            expected_water_pressures = [
                result.value for result in expected_results
            ]
            self.assertVectorAlmostEqual(
                water_pressures, expected_water_pressures, places=None, delta=10.0
            )

        if test_helper.want_test_plots():
            plot_times = [12000, 24000, 36000, 48000, 72000, 96000, 192000]
            data_series_collection = []
            for time in plot_times:
                water_pressures = reader.nodal_values_at_time(
                    "WATER_PRESSURE", time, output_data, node_ids_left_boundary
                )
                sorted_y, sorted_data = zip(
                    *sorted(zip(depth_boundary_nodes, water_pressures))
                )
                list = zip(sorted_data, sorted_y)
                data_series_collection.append(
                    plot_utils.DataSeries(
                        list, label=f"Time = {time}", line_style="-", marker=""
                    )
                )
            asserted_data_points = []
            for time, expected_results in expected_results_at_times.items():
                for expected_result in expected_results:
                    water_pressure = expected_result.value

                    asserted_data_points.append((
                        water_pressure,
                        -1.0 * simulation.model.GetModelPart(
                            "PorousDomain.porous_computational_model_part"
                        ).GetNode(expected_result.node_id).Y
                    ))
            data_series_collection.append(
                plot_utils.DataSeries(
                    asserted_data_points, label=f"Asserted pressures", line_style="", marker="x", color="r"
                )
            )
            plot_utils._make_plot(
                data_series_collection,
                os.path.join(file_path, "infiltration_from_top_boundary.svg"),
                xlabel="water pressure [Pa]",
                ylabel="depth [m]",
                yaxis_inverted=True,
            )

    # def test_infiltration_from_top_boundary_upw(self):
    #     file_path = test_helper.get_file_path(
    #         os.path.join("test_partially_saturated", "test_infilteration_upw")
    #     )
    #     simulation = test_helper.run_kratos(file_path)


if __name__ == "__main__":
    KratosUnittest.main()
