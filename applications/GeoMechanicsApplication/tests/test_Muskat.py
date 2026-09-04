import csv
import os
import math

import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_helper
import KratosGeoUnittest

if test_helper.want_test_plots():
    import KratosMultiphysics.GeoMechanicsApplication.geo_plot_utilities as plot_utils

from dataclasses import dataclass

from KratosMultiphysics.GeoMechanicsApplication.gid_output_file_reader import (
    GiDOutputFileReader,
)
from KratosMultiphysics.GeoMechanicsApplication.unit_conversions import (
    fraction_to_percentage,
)


class KratosGeoMechanicsMuskatTests(KratosGeoUnittest.TestCase):
    """
    A test suite representing (semi)-Muskat test cases. Instead of a seepage boundary, a fixed head is used on the
    right boundary. For more information on the test case, see the README.md file in the corresponding test folder.
    """

    def _assert_muskat_results(
        self, parent_name, test_name, expected_results_for_variable
    ):
        file_path = test_helper.get_file_path(os.path.join(parent_name, test_name))
        simulation = test_helper.run_kratos(file_path)

        reader = GiDOutputFileReader()
        output_data = reader.read_output_from(
            os.path.join(file_path, "Muskat_Dirichlet.post.res")
        )
        if test_helper.want_test_plots():
            self._create_effective_saturation_plot(
                expected_results_for_variable, file_path, output_data, simulation
            )
            self._create_phreatic_line_plot(file_path, output_data, simulation)
            self._create_fluid_flux_plot(
                expected_results_for_variable, file_path, output_data, simulation
            )

        for variable_name, expected_results in expected_results_for_variable.items():
            actual_results = GiDOutputFileReader.nodal_values_at_time(
                variable_name,
                1.0,
                output_data,
                [result.node_id for result in expected_results],
            )
            self.assert_nodal_values_at_time(
                output_data,
                variable_name,
                [result.node_id for result in expected_results],
                [expected.value for expected in expected_results],
                1.0,
            )

    def _create_effective_saturation_plot(
        self, expected_results_for_variable, file_path, output_data, simulation
    ):
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
        saturations = [fraction_to_percentage(saturation) for saturation in saturations]
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

        data_points = []
        with open(
            os.path.join(file_path, "expected_saturation_at_x_1_52.csv"),
            newline="",
        ) as csv_file:
            reader = csv.DictReader(csv_file, skipinitialspace=True)
            data_points = [
                (float(row["Effective Saturation"]), float(row["Y"])) for row in reader
            ]
        data_series_collection.append(
            plot_utils.DataSeries(
                data_points,
                label="Commercial FE package",
            )
        )

        asserted_data_points = []
        for expected_result in expected_results_for_variable["EFFECTIVE_SATURATION"]:
            effective_saturation = fraction_to_percentage(expected_result.value)

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
            os.path.join(file_path, "saturation_for_x_1_52.svg"),
            xlabel="Effective Saturation [%]",
            ylabel="Y [m]",
            yaxis_inverted=False,
        )

    def _create_fluid_flux_plot(
        self, expected_results_for_variable, file_path, output_data, simulation
    ):
        data_series_collection = []
        y_coord_by_id_for_right_boundary_nodes = {}
        for node in simulation.model.GetModelPart(
            "PorousDomain.porous_computational_model_part"
        ).Nodes:
            if abs(node.X - 1.52) < 0.01:
                y_coord_by_id_for_right_boundary_nodes[node.Id] = node.Y
        fluxes = GiDOutputFileReader.nodal_values_at_time(
            "FLUID_FLUX_VECTOR",
            1.0,
            output_data,
            y_coord_by_id_for_right_boundary_nodes.keys(),
        )
        flux_length = [
            math.sqrt(flux[0] ** 2 + flux[1] ** 2 + flux[2]**2) * 86400
            for flux in fluxes
        ]
        sorted_depth, sorted_data = zip(
            *sorted(
                zip(
                    y_coord_by_id_for_right_boundary_nodes.values(),
                    flux_length,
                )
            )
        )
        data = zip(sorted_data, sorted_depth)
        data_series_collection.append(
            plot_utils.DataSeries(data, label="Kratos", line_style="-", marker="")
        )

        data_points = []
        with open(
            os.path.join(file_path, "expected_flux_at_x_1_52.csv"),
            newline="",
        ) as csv_file:
            reader = csv.DictReader(csv_file, skipinitialspace=True)
            data_points = [(float(row["Flux"]), float(row["Y"])) for row in reader]
        data_series_collection.append(
            plot_utils.DataSeries(
                data_points,
                label="Commercial FE package",
            )
        )

        plot_utils._make_plot(
            data_series_collection,
            os.path.join(file_path, "fluxes_for_x_1_52.svg"),
            xlabel="|Flux| [m/d]",
            ylabel="Y [m]",
            yaxis_inverted=False,
        )

    def _create_phreatic_line_plot(self, file_path, output_data, simulation):
        import matplotlib.pyplot as plt
        from matplotlib.tri import Triangulation

        model_part = simulation.model.GetModelPart(
            "PorousDomain.porous_computational_model_part"
        )
        nodes = list(model_part.Nodes)
        node_ids = [node.Id for node in nodes]
        xs = [node.X for node in nodes]
        ys = [node.Y for node in nodes]
        pressures = GiDOutputFileReader.nodal_values_at_time(
            "WATER_PRESSURE",
            1.0,
            output_data,
            node_ids,
        )

        # Although this triangulation (Delaunay triangulation) is not identical to
        # the mesh used in the simulation, it was decided, it is sufficient for
        # visualizing the p = 0 line.
        tri = Triangulation(xs, ys)
        plt.figure(figsize=(10, 8))
        contour_zero = plt.tricontour(tri, pressures, levels=[0.0])

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
        with open(
            os.path.join(file_path, "expected_phreatic_line.csv"), newline=""
        ) as csv_file:
            reader = csv.DictReader(csv_file, skipinitialspace=True)
            data_points = sorted(
                [(float(row["X"]), float(row["Y"])) for row in reader],
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

    @dataclass
    class ExpectedResult:
        node_id: int
        value: float | list[float]

    def test_muskat_van_genuchten_hydrostatic(self):

        expected_results_for_variables = {
            "WATER_PRESSURE": [
                self.ExpectedResult(node_id=198, value=7939.92),
                self.ExpectedResult(node_id=244, value=-2753.0),
                self.ExpectedResult(node_id=183, value=-20549.3),
            ],
            "EFFECTIVE_SATURATION": [
                self.ExpectedResult(node_id=273, value=1.0),
                self.ExpectedResult(node_id=890, value=0.485245),
                self.ExpectedResult(node_id=1440, value=0.991294),
            ],
            "FLUID_FLUX_VECTOR": [
                self.ExpectedResult(node_id=273, value=[8.1665e-06, -1.08412e-06, 0.0]),
                self.ExpectedResult(node_id=890, value=[4.6101e-09, 4.7292e-11, 0.0]),
                self.ExpectedResult(
                    node_id=1440, value=[2.37856e-06, -1.57026e-06, 0.0]
                ),
            ],
        }
        self._assert_muskat_results(
            "Muskat", "van_genuchten_hydrostatic", expected_results_for_variables
        )

    def test_muskat_saturated_hydrostatic(self):
        expected_results_for_variables = {
            "WATER_PRESSURE": [
                self.ExpectedResult(node_id=198, value=14346.4),
                self.ExpectedResult(node_id=244, value=-37.4452),
                self.ExpectedResult(node_id=183, value=-20335.1),
            ],
            "EFFECTIVE_SATURATION": [
                self.ExpectedResult(node_id=273, value=1.0),
                self.ExpectedResult(node_id=890, value=1.0),
                self.ExpectedResult(node_id=1440, value=1.0),
            ],
            "FLUID_FLUX_VECTOR": [
                self.ExpectedResult(
                    node_id=273, value=[6.00133e-06, -1.34793e-11, 0.0]
                ),
                self.ExpectedResult(
                    node_id=890, value=[6.00131e-06, -8.38038e-12, 0.0]
                ),
                self.ExpectedResult(
                    node_id=1440, value=[6.00133e-06, -1.35382e-11, 0.0]
                ),
            ],
        }
        self._assert_muskat_results(
            "Muskat", "saturated_hydrostatic", expected_results_for_variables
        )


if __name__ == "__main__":
    KratosUnittest.main()
