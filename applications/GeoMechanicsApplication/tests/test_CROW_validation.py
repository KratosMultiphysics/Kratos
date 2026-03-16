import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.GeoMechanicsApplication.geomechanics_analysis as analysis
import KratosMultiphysics.GeoMechanicsApplication.context_managers as context_managers
from KratosMultiphysics.GeoMechanicsApplication.gid_output_file_reader import (
    GiDOutputFileReader,
)
from KratosMultiphysics.GeoMechanicsApplication.unit_conversions import unit_to_k_unit
import test_helper

import os
import json
from pathlib import Path

import KratosMultiphysics.GeoMechanicsApplication.geo_plot_utilities as plot_utils


def get_sheetpile_node_ids():
    return [
        411,
        423,
        435,
        454,
        472,
        489,
        511,
        525,
        544,
        560,
        579,
        597,
        615,
        628,
        648,
        663,
        685,
        703,
        715,
        740,
        758,
        768,
        778,
        799,
        817,
        827,
        840,
    ]


class KratosGeoMechanicsCrowValidation(KratosUnittest.TestCase):
    def setUp(self):
        super().setUp()

        self.stages_info = {
            "initial_stage": {"end_time": -1.0, "base_name": "1_Initial_stage"},
            "null_step": {"end_time": 0.0, "base_name": "2_Null_step"},
            "final_equilibrium": {"end_time": 1.0, "base_name": "3_Final_equilibrium"},
        }

    def run_simulation_and_checks(self, sub_directory_name):
        project_path = test_helper.get_file_path(
            os.path.join("crow_validation", sub_directory_name)
        )

        with context_managers.set_cwd_to(project_path):
            model = Kratos.Model()
            for stage in self.stages_info.values():
                with open(f"{stage['base_name']}.json", "r") as f:
                    stage_parameters = Kratos.Parameters(f.read())
                analysis.GeoMechanicsAnalysis(model, stage_parameters).Run()

        self.create_sheetpile_plots(project_path)

    def get_plot_stages(self):
        return [self.stages_info["final_equilibrium"]]

    def get_y_coords(self, project_path, stage_base_name, node_ids):
        coordinates = test_helper.read_coordinates_from_post_msh_file(
            Path(project_path) / "gid_output" / f"{stage_base_name}.post.msh", node_ids=node_ids
        )
        return [coord[1] for coord in coordinates]

    def get_wall_variable_collections_per_stage(
            self,
            kratos_variable_label,
            variable_plot_label,
            project_path,
            plot_stages,
    ):
        node_ids = get_sheetpile_node_ids()

        y_coords = self.get_y_coords(
            project_path, plot_stages[0]["base_name"], node_ids
        )

        variable_data_collections = []
        reader = GiDOutputFileReader()

        for stage in plot_stages:
            output_data = reader.read_output_from(
                Path(project_path) / "gid_output" / f"{stage['base_name']}.post.res"
            )

            variable_kratos_data = reader.nodal_values_at_time(
                kratos_variable_label,
                stage["end_time"],
                output_data,
                node_ids=node_ids,
            )

            variable_kratos_data = [unit_to_k_unit(value) for value in variable_kratos_data]

            data_series_collection = [
                plot_utils.DataSeries(
                    zip(variable_kratos_data, y_coords),
                    f"{variable_plot_label} [Kratos]",
                    line_style="-",
                    marker=".",
                )
            ]
            variable_data_collections.append(data_series_collection)

        return variable_data_collections

    def get_wall_horizontal_displacement_collections_per_stage(
            self,
            project_path,
            plot_stages,
    ):
        node_ids = get_sheetpile_node_ids()

        y_coords = self.get_y_coords(
            project_path, plot_stages[0]["base_name"], node_ids
        )

        variable_data_collections = []
        reader = GiDOutputFileReader()

        for stage in plot_stages:
            output_data = reader.read_output_from(
                Path(project_path) / "gid_output" / f"{stage['base_name']}.post.res"
            )

            displacements = reader.nodal_values_at_time(
                "DISPLACEMENT",
                stage["end_time"],
                output_data,
                node_ids=node_ids,
            )

            horizontal_displacements = [value[0] for value in displacements]

            data_series_collection = [
                plot_utils.DataSeries(
                    zip(horizontal_displacements, y_coords),
                    "Horizontal displacement [Kratos]",
                    line_style="-",
                    marker=".",
                )
            ]
            variable_data_collections.append(data_series_collection)

        return variable_data_collections

    def create_sheetpile_plots(self, project_path):
        plot_stages = self.get_plot_stages()
        plot_titles = [stage["base_name"] for stage in plot_stages]

        bending_moment_collections = self.get_wall_variable_collections_per_stage(
            "BENDING_MOMENT",
            "Bending moment",
            project_path,
            plot_stages,
        )
        plot_utils.make_sub_plots(
            bending_moment_collections,
            Path(project_path) / "bending_moments.svg",
            titles=plot_titles,
            xlabel="Bending moment [kNm/m]",
            ylabel="y [m]",
            )

        shear_force_collections = self.get_wall_variable_collections_per_stage(
            "SHEAR_FORCE",
            "Shear force",
            project_path,
            plot_stages,
        )
        plot_utils.make_sub_plots(
            shear_force_collections,
            Path(project_path) / "shear_forces.svg",
            titles=plot_titles,
            xlabel="Shear force [kN/m]",
            ylabel="y [m]",
            )

        horizontal_displacement_collections = (
            self.get_wall_horizontal_displacement_collections_per_stage(
                project_path,
                plot_stages,
            )
        )
        plot_utils.make_sub_plots(
            horizontal_displacement_collections,
            Path(project_path) / "horizontal_displacements.svg",
            titles=plot_titles,
            xlabel="Horizontal displacement [m]",
            ylabel="y [m]",
            )

    def test_simulation_without_excavation(self):
        self.run_simulation_and_checks("without_excavation")


if __name__ == "__main__":
    KratosUnittest.main()
