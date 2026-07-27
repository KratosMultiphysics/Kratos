import os

import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_helper

if test_helper.want_test_plots():
    import KratosMultiphysics.GeoMechanicsApplication.geo_plot_utilities as plot_utils

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
        line, index_of_x=1, index_of_y=0, x_transform=lambda x : x / 100.0
    )
class KratosGeoMechanicsMuskatTests(KratosUnittest.TestCase):
    """
    To be detailed out
    """




    def _assert_fully_saturated_flow(
        self, parent_name, test_name, node_number, expected_value
    ):
        file_path = test_helper.get_file_path(os.path.join(parent_name, test_name))
        simulation = test_helper.run_kratos(file_path)

        if test_helper.want_test_plots():
            reader = GiDOutputFileReader()
            output_data = reader.read_output_from(
                os.path.join(file_path, "not_Muskat_yet.post.res")
            )
            data_series_collection = []
            y_coord_by_id_for_right_boundary_nodes = {}
            for node in simulation.model.GetModelPart(
                    "PorousDomain.porous_computational_model_part"
            ).Nodes:
                if abs(node.X - 1.62) < 0.01:
                    y_coord_by_id_for_right_boundary_nodes[node.Id] = node.Y
            water_pressures = GiDOutputFileReader.nodal_values_at_time(
                "DEGREE_OF_SATURATION",
                1.0,
                output_data,
                y_coord_by_id_for_right_boundary_nodes.keys(),
            )
            sorted_depth, sorted_data = zip(
                *sorted(
                    zip(
                        y_coord_by_id_for_right_boundary_nodes.values(),
                        water_pressures,
                    )
                )
            )
            data = zip(sorted_data, sorted_depth)
            data_series_collection.append(
                plot_utils.DataSeries(
                    data, label="Kratos", line_style="-", marker=""
                )
            )

            data_points = test_helper.get_data_points_from_file(
                os.path.join(file_path,"expected_results.csv"), extract_saturation_and_y_from_line
            )
            data_series_collection.append(
                plot_utils.DataSeries(
                    data_points,
                    label="Commercial FE package",
                )
            )
            plot_utils._make_plot(
                data_series_collection,
                os.path.join(file_path, "saturation_on_right_boundary.svg"),
                xlabel="Degree of Saturation [-]",
                ylabel="Y [m]",
                yaxis_inverted=False,
            )



        water_pressure = [
            water_pressure
            for water_pressure in test_helper.get_water_pressure(simulation)
        ]
        self.assertAlmostEqual(expected_value, water_pressure[node_number - 1], 6)

    # To be changed, node 100 is a random pick and may be on the prescribed boundary. pick a better result, with computed result
    # def test_fully_saturated_hydrostatic(self):
    #     self._assert_fully_saturated_flow(
    #         "Muskat", "fully_saturated_hydrostatic", 100, 22800.0
    #     )
    #
    # def test_fully_saturated_hydrostatic_cutoff(self):
    #     self._assert_fully_saturated_flow(
    #         "Muskat", "fully_saturated_hydrostatic_cutoff", 100, 0.0
    #     )
    def test_van_genuchten_hydrostatic(self):
        self._assert_fully_saturated_flow(
            "Muskat", "van_genuchten_hydrostatic", 100, 11400
        )

if __name__ == "__main__":
    KratosUnittest.main()
