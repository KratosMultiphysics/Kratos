import os

import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_helper

if test_helper.want_test_plots():
    import KratosMultiphysics.GeoMechanicsApplication.geo_plot_utilities as plot_utils

from dataclasses import dataclass

from KratosMultiphysics.GeoMechanicsApplication.gid_output_file_reader import (
    GiDOutputFileReader,
)


def _extract_x_and_y_from_line(line, index_of_x=0, index_of_y=1, x_transform=None):
    line = line.strip().lstrip("\ufeff").replace("ï»¿", "")
    words = [
        word.strip() for word in (line.split(",") if "," in line else line.split())
    ]

    x_ = float(words[index_of_x])
    if x_transform:
        x_ = x_transform(x_)
    y_ = float(words[index_of_y])

    return x_, y_


def extract_saturation_and_y_from_line(line):
    return _extract_x_and_y_from_line(
        line, index_of_x=1, index_of_y=0, x_transform=lambda x: x
    )


def extract_x_and_y_from_line(line):
    return _extract_x_and_y_from_line(line, index_of_x=0, index_of_y=1)


class KratosGeoMechanicsMuskatTests(KratosUnittest.TestCase):
    """
    To be detailed out
    """

    def _assert_fully_saturated_flow(
        self, parent_name, test_name, expected_results_for_variable
    ):
        file_path = test_helper.get_file_path(os.path.join(parent_name, test_name))
        simulation = test_helper.run_kratos(file_path)

        reader = GiDOutputFileReader()
        output_data = reader.read_output_from(
            os.path.join(file_path, "not_Muskat_yet.post.res")
        )
        if test_helper.want_test_plots():
            data_series_collection = []
            y_coord_by_id_for_right_boundary_nodes = {}
            for node in simulation.model.GetModelPart(
                "PorousDomain.porous_computational_model_part"
            ).Nodes:
                if abs(node.X - 1.52) < 0.01:
                    y_coord_by_id_for_right_boundary_nodes[node.Id] = node.Y
            saturations = GiDOutputFileReader.nodal_values_at_time(
                "EFFECTIVE_SATURATION",
                1.0,
                output_data,
                y_coord_by_id_for_right_boundary_nodes.keys(),
            )
            print(y_coord_by_id_for_right_boundary_nodes.keys())
            saturations = [saturation * 100 for saturation in saturations]
            sorted_depth, sorted_data = zip(
                *sorted(
                    zip(
                        y_coord_by_id_for_right_boundary_nodes.values(),
                        saturations,
                    )
                )
            )
            data = zip(sorted_data, sorted_depth)
            data_series_collection.append(
                plot_utils.DataSeries(data, label="Kratos", line_style="-", marker="")
            )

            data_points = test_helper.get_data_points_from_file(
                os.path.join(file_path, "expected_saturation_at_x_1_52.csv"),
                extract_saturation_and_y_from_line,
            )
            data_series_collection.append(
                plot_utils.DataSeries(
                    data_points,
                    label="Commercial FE package",
                )
            )

            asserted_data_points = []
            for expected_result in expected_results_for_variable["EFFECTIVE_SATURATION"]:
                effective_saturation = expected_result.value * 100

                asserted_data_points.append(
                    (
                        effective_saturation,
                        y_coord_by_id_for_right_boundary_nodes[expected_result.node_id],
                    )
                )
            data_series_collection.append(
                plot_utils.DataSeries(
                    asserted_data_points,
                    label=f"Asserted effective saturation",
                    line_style="",
                    marker="x",
                )
            )
            plot_utils._make_plot(
                data_series_collection,
                os.path.join(file_path, "saturation_on_right_boundary.svg"),
                xlabel="Degree of Saturation [%]",
                ylabel="Y [m]",
                yaxis_inverted=False,
            )

        xs = []
        ys = []
        pressures = []
        for node in simulation.model.GetModelPart(
            "PorousDomain.porous_computational_model_part"
        ).Nodes:
            xs.append(node.X)
            ys.append(node.Y)
            pressures.append(
                GiDOutputFileReader.nodal_values_at_time(
                    "WATER_PRESSURE",
                    1.0,
                    output_data,
                    [node.Id],
                )[0]
            )

        if test_helper.want_test_plots():
            import matplotlib.pyplot as plt
            import numpy as np
            from matplotlib.tri import Triangulation

            tri = Triangulation(xs, ys)
            plt.figure(figsize=(10, 8))
            contour_zero = plt.tricontour(tri, pressures, levels=[0])

            # Extract isoline coordinates
            isoline_coords = []
            if hasattr(contour_zero, "allsegs"):
                for level_segs in contour_zero.allsegs:
                    for seg in level_segs:
                        isoline_coords.append(seg)

            data_series_collection = [
                plot_utils.DataSeries(
                    isoline_coords[0],
                    label="p = 0 line Kratos",
                    line_style="-",
                    marker="",
                )
            ]
            data_points = sorted(
                test_helper.get_data_points_from_file(
                    os.path.join(file_path, "expected_phreatic_line.csv"),
                    extract_x_and_y_from_line,
                ),
                key=lambda point: point[1],
            )

            data_series_collection.append(
                plot_utils.DataSeries(
                    data_points,
                    label="p = 0 line commercial FE",
                )
            )

            plot_utils._make_plot(
                data_series_collection,
                os.path.join(file_path, "phreatic_line.svg"),
                xlabel="X [m]",
                ylabel="Y [m]",
                yaxis_inverted=False,
            )

        for variable_name, expected_results in expected_results_for_variable.items():
            actual_results = GiDOutputFileReader.nodal_values_at_time(
                variable_name,
                1.0,
                output_data,
                [result.node_id for result in expected_results],
            )
            for idx, expected_result in enumerate(expected_results):
                self.assertAlmostEqual(expected_result.value, actual_results[idx], 6)

    def test_van_genuchten_hydrostatic(self):

        @dataclass
        class ExpectedResult:
            node_id: int
            value: float

        expected_results_for_variables = {
            "WATER_PRESSURE": [
                ExpectedResult(node_id=198, value=8325.15),
                ExpectedResult(node_id=244, value=-3031.34),
                ExpectedResult(node_id=183, value=-20243.4),
            ],
            "EFFECTIVE_SATURATION": [
                ExpectedResult(node_id=273, value=1.0),
                ExpectedResult(node_id=890, value=0.462247),
                ExpectedResult(node_id=1440, value=0.924914),
            ],
        }
        self._assert_fully_saturated_flow(
            "Muskat", "van_genuchten_hydrostatic", expected_results_for_variables
        )


if __name__ == "__main__":
    KratosUnittest.main()
