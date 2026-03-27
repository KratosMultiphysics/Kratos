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
        3966,
        3981,
        3996,
        4010,
        4024,
        4036,
        4052,
        4069,
        4076,
        4093,
        4105,
        4123,
        4138,
        4150,
        4164,
        4176,
        4192,
        4207,
        4222,
        4235,
        4245,
        4260,
        4275,
        4284,
        4301,
        4312,
        4326,
        4342,
        4356,
        4365,
        4379,
        4394,
        4407,
        4419,
        4434,
        4448,
        4461,
        4474,
        4487,
        4501,
        4513,
        4526,
        4540,
        4551,
        4564,
        4575,
        4586,
        4599,
        4608,
        4618,
        4628,
        4637,
        4644,
    ]
def get_soil_side_node_ids_of_right_interfaces():
    return [
        4765,
        4766,
        4767,
        4768,
        4769,
        4770,
        4771,
        4772,
        4773,
        4774,
        4775,
        4776,
        4777,
        4778,
        4779,
        4780,
        4781,
        4782,
        4783,
        4784,
        4785,
        4786,
        4787,
        4788,
        4789,
        4790,
        4791,
        4792,
        4793,
        4794,
        4795,
        4796,
        4797,
        4798,
        4799,
        4800,
        4801,
        4802,
        4803,
        4804,
        4805,
        4806,
        4807,
        4808,
        4809,
        4810,
        4811,
        4812,
        4813,
        4814,
        4815,
        4816,
        4817,
    ]

def _extract_x_and_y_from_line(line, index_of_x=0, index_of_y=1, x_transform=None):
    line = line.strip().lstrip("\ufeff").replace("ï»¿", "")
    words = [word.strip() for word in (line.split(",") if "," in line else line.split())]

    x_ = float(words[index_of_x])
    if x_transform:
        x_ = x_transform(x_)
    y_ = float(words[index_of_y])

    return x_, y_

def extract_bending_moment_and_y_from_line(line):
    return _extract_x_and_y_from_line(line, index_of_x=1, index_of_y=0)


def extract_horizontal_displacements_from_line(line):
    return _extract_x_and_y_from_line(line, index_of_x=3, index_of_y=0, x_transform=lambda x: x / 1000.0)


def extract_shear_force_and_y_from_line(line):
    return _extract_x_and_y_from_line(
        line, index_of_x=2, index_of_y=0
    )

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
        self.create_interface_plots(project_path)

    def read_json_output(self, project_path, stage):
        with open(
                os.path.join(project_path, f"{stage['base_name']}_interface_output.json"),
                "r",
        ) as output_file:
            return json.load(output_file)

    def get_variable_collections_per_stage(
            self,
            kratos_variable_label,
            variable_plot_label,
            project_path,
            structural_stages,
            object_name
    ):
        node_ids = get_soil_side_node_ids_of_right_interfaces()


        # Since the coordinates do not change between stages, we base them on the first stage
        y_coords = self.get_y_coords(
            project_path, structural_stages[0]["base_name"], node_ids
        )

        variable_data_collections = []
        for stage in structural_stages:
            json_data = self.read_json_output(project_path, stage)
            variable_kratos_data = []
            index = 0 if variable_plot_label == "Normal traction" else 1
            for node_label in [f"NODE_{node_id}" for node_id in node_ids]:
                variable_kratos_data.append(
                    json_data[node_label][kratos_variable_label][0][index]
                )

            variable_kratos_data = [
                unit_to_k_unit(value) for value in variable_kratos_data
            ]
            sorted_y, sorted_data = zip(*sorted(zip(y_coords, variable_kratos_data)))
            data_series_collection = []
            data_series_collection.append(
                plot_utils.DataSeries(
                    zip(sorted_data, sorted_y),
                    f"{variable_plot_label} [Kratos]",
                    line_style="-",
                    marker=".",
                )
            )
            variable_data_collections.append(data_series_collection)
        return variable_data_collections

    def create_interface_plots(self, project_path):
        structural_stages = [self.stages_info["final_equilibrium"]]

        normal_traction_kratos_label = "GEO_EFFECTIVE_TRACTION_VECTOR"
        normal_traction_plot_label = "Normal traction"
        normal_traction_collections = self.get_variable_collections_per_stage(
            normal_traction_kratos_label,
            normal_traction_plot_label,
            project_path,
            structural_stages,
            "right_interface"
        )

        plot_titles = [stage["base_name"] for stage in structural_stages]
        plot_utils.make_sub_plots(
            normal_traction_collections,
            Path(project_path) / "normal_tractions_all_stages.svg",
            titles=plot_titles,
            xlabel=r"Normal Traction [$\mathrm{kN} / \mathrm{m}^2$]",
            ylabel="y [m]",
            # figsize=(4, 6),
        )

        shear_traction_plot_label = "Shear traction"
        shear_traction_kratos_label = "GEO_EFFECTIVE_TRACTION_VECTOR"
        shear_traction_collections = self.get_variable_collections_per_stage(
            shear_traction_kratos_label,
            shear_traction_plot_label,
            project_path,
            structural_stages,
            "left_interface",
        )

        plot_utils.make_sub_plots(
            shear_traction_collections,
            Path(project_path) / "shear_tractions_all_stages.svg",
            titles=plot_titles,
            xlabel=r"Shear Traction [$\mathrm{kN} / \mathrm{m}^2$]",
            ylabel="y [m]",
            # figsize=(4, 6),
            )

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
            data_extractor
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

            if (Path(project_path) / f"{stage['base_name']}_comparison.csv").exists():
                comparison_variable = test_helper.get_data_points_from_file(
                    Path(project_path) / f"{stage['base_name']}_comparison.csv", data_extractor
                )
                data_series_collection.append(
                    plot_utils.DataSeries(
                        comparison_variable,
                        f"{variable_plot_label} [Comparison]",
                        marker="1",
                    )
                )

            fem_comparison_csv = Path(project_path) / f"{stage['base_name']}_comparison_FEM.csv"
            if fem_comparison_csv.exists():
                fem_comparison_variable = test_helper.get_data_points_from_file(
                    fem_comparison_csv, data_extractor
                )
                data_series_collection.append(
                    plot_utils.DataSeries(
                        fem_comparison_variable,
                        f"{variable_plot_label} [FEM Comparison]",
                        marker="2",
                    )
                )

            fem_comparison_csv = Path(project_path) / f"{stage['base_name']}_comparison_FEM_with_excavation_stages.csv"
            if fem_comparison_csv.exists():
                fem_comparison_variable = test_helper.get_data_points_from_file(
                    fem_comparison_csv, data_extractor
                )
                data_series_collection.append(
                    plot_utils.DataSeries(
                        fem_comparison_variable,
                        f"{variable_plot_label} [FEM Comparison_with excavation stages]",
                        marker="3",
                    )
                )

            variable_data_collections.append(data_series_collection)

        return variable_data_collections

    def get_wall_horizontal_displacement_collections_per_stage(
            self,
            project_path,
            plot_stages,
            data_extractor
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

            if (Path(project_path) / f"{stage['base_name']}_comparison.csv").exists():
                comparison_variable = test_helper.get_data_points_from_file(
                    Path(project_path) / f"{stage['base_name']}_comparison.csv", data_extractor
                )
                data_series_collection.append(
                    plot_utils.DataSeries(
                        comparison_variable,
                        "Horizontal displacement [Comparison]",
                        marker="1",
                    )
                )

            if (Path(project_path) / f"{stage['base_name']}_comparison_FEM.csv").exists():
                comparison_variable = test_helper.get_data_points_from_file(
                    Path(project_path) / f"{stage['base_name']}_comparison_FEM.csv", data_extractor
                )
                data_series_collection.append(
                    plot_utils.DataSeries(
                        comparison_variable,
                        "Horizontal displacement [FEM Comparison]",
                        marker="2",
                    )
                )

            if (Path(project_path) / f"{stage['base_name']}_comparison_FEM_with_excavation_stages.csv").exists():
                comparison_variable = test_helper.get_data_points_from_file(
                    Path(project_path) / f"{stage['base_name']}_comparison_FEM_with_excavation_stages.csv", data_extractor
                )
                data_series_collection.append(
                    plot_utils.DataSeries(
                        comparison_variable,
                        "Horizontal displacement [FEM Comparison_with_excavation_stages]",
                        marker="3",
                    )
                )

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
            extract_bending_moment_and_y_from_line
        )
        plot_utils.make_sub_plots(
            bending_moment_collections,
            Path(project_path) / "bending_moments.svg",
            titles=plot_titles,
            xlabel="Bending moment [kNm/m]",
            ylabel="y [m]",
            # figsize=(4, 6),
            )

        shear_force_collections = self.get_wall_variable_collections_per_stage(
            "SHEAR_FORCE",
            "Shear force",
            project_path,
            plot_stages,
            extract_shear_force_and_y_from_line
        )
        plot_utils.make_sub_plots(
            shear_force_collections,
            Path(project_path) / "shear_forces.svg",
            titles=plot_titles,
            xlabel="Shear force [kN/m]",
            ylabel="y [m]",
            # figsize=(4, 6),
            )

        horizontal_displacement_collections = (
            self.get_wall_horizontal_displacement_collections_per_stage(
                project_path,
                plot_stages,
                extract_horizontal_displacements_from_line
            )
        )
        plot_utils.make_sub_plots(
            horizontal_displacement_collections,
            Path(project_path) / "horizontal_displacements.svg",
            titles=plot_titles,
            xlabel="Horizontal displacement [m]",
            ylabel="y [m]",
            # figsize=(4, 6),
        )

    def test_simulation_without_excavation(self):
        self.run_simulation_and_checks("without_excavation")


if __name__ == "__main__":
    KratosUnittest.main()
