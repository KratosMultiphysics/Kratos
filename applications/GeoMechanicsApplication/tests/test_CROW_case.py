import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.project import Project
import importlib
import KratosMultiphysics.GeoMechanicsApplication.context_managers as context_managers
from KratosMultiphysics.GeoMechanicsApplication.gid_output_file_reader import (
    GiDOutputFileReader,
)
from KratosMultiphysics.GeoMechanicsApplication.unit_conversions import unit_to_k_unit
import test_helper

import argparse
import csv
import os
import json
from pathlib import Path
import sys

if test_helper.want_test_plots():
    import KratosMultiphysics.GeoMechanicsApplication.geo_plot_utilities as plot_utils


csv_fieldname_node = "node"
csv_fieldname_bending_moment = "bending_moment_in_Nm_per_m"
csv_fieldname_shear_force = "shear_force_in_N_per_m"
csv_fieldname_normal_force = "normal_force_in_N_per_m"
csv_fieldname_horizontal_total_displacement = "horizontal_total_displacement_in_m"

csv_fieldnames = [
    csv_fieldname_node,  # this one must come first
    csv_fieldname_bending_moment,
    csv_fieldname_shear_force,
    csv_fieldname_normal_force,
    csv_fieldname_horizontal_total_displacement,
]


def get_sheetpile_node_ids():
    return [
        3968,
        3982,
        3996,
        4010,
        4026,
        4036,
        4051,
        4066,
        4075,
        4090,
        4102,
        4120,
        4135,
        4148,
        4161,
        4172,
        4188,
        4204,
        4218,
        4231,
        4241,
        4256,
        4271,
        4280,
        4297,
        4308,
        4322,
        4338,
        4352,
        4361,
        4375,
        4390,
        4403,
        4415,
        4430,
        4444,
        4457,
        4470,
        4483,
        4497,
        4509,
        4522,
        4536,
        4547,
        4560,
        4571,
        4582,
        4595,
        4604,
        4614,
        4624,
        4633,
        4639,
    ]


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


def extract_bending_moment_and_y_from_line(line):
    return _extract_x_and_y_from_line(
        line, index_of_x=1, index_of_y=0, x_transform=lambda x: -x
    )


def extract_horizontal_displacements_from_line(line):
    return _extract_x_and_y_from_line(
        line, index_of_x=4, index_of_y=0, x_transform=lambda x: x / 1000.0
    )


def extract_horizontal_displacements_from_dsheetpiling_line(line):
    return _extract_x_and_y_from_line(
        line, index_of_x=3, index_of_y=0, x_transform=lambda x: x / 1000.0
    )


def extract_shear_force_and_y_from_line(line):
    # The shear force sign in the comparison data is opposite to the Kratos sign
    return _extract_x_and_y_from_line(
        line, index_of_x=2, index_of_y=0, x_transform=lambda x: -x
    )


def extract_normal_force_and_y_from_line(line):
    return _extract_x_and_y_from_line(line, index_of_x=5, index_of_y=0)


def extract_normal_traction_and_y_from_line(line):
    return _extract_x_and_y_from_line(line, index_of_x=1, index_of_y=0)


def extract_shear_traction_and_y_from_line(line):
    return _extract_x_and_y_from_line(
        line,
        index_of_x=2,
        index_of_y=0,
    )


def get_expected_results_from_csv(csv_filepath):
    result_fieldnames = csv_fieldnames[1:]  # node ID is not a result
    with open(csv_filepath, newline="") as csv_file:
        reader = csv.DictReader(csv_file)
        return {
            int(row[csv_fieldname_node]): {
                fieldname: float(row[fieldname]) for fieldname in result_fieldnames
            }
            for row in reader
        }


class KratosGeoMechanicsCrowValidation(KratosUnittest.TestCase):
    def setUp(self):
        super().setUp()

        self.stages_info = None  # Will be populated by the specific simulation runs

    def run_staged_construction_analysis_and_checks(self, material_model_dir_name):
        # fmt: off
        self.stages_info = {
            "initial_stage":      {"end_time": -1.0, "base_name": "1_Initial_stage"},
            "null_step":          {"end_time":  0.0, "base_name": "2_Null_step"},
            "wall_installation":  {"end_time":  1.0, "base_name": "3_Wall_installation"},
            "first_excavation":   {"end_time":  2.0, "base_name": "4_First_excavation"},
            "strut_installation": {"end_time":  3.0, "base_name": "5_Anchor_installation"},
            "second_excavation":  {"end_time":  4.0, "base_name": "6_Second_excavation"},
            "third_excavation":   {"end_time":  5.0, "base_name": "7_Third_excavation"},
        }
        # fmt: on

        self.run_simulation_and_checks(
            Path(material_model_dir_name) / "staged_construction",
            "staged_construction.json",
        )

    def run_simulation_and_checks(self, relative_test_path, analysis_filename):
        project_path = test_helper.get_file_path(
            Path("crow_validation") / relative_test_path
        )

        with context_managers.set_cwd_to(project_path):
            with open(
                Path("..") / ".." / "common" / analysis_filename, "r"
            ) as parameter_file:
                project_parameters = Kratos.Parameters(parameter_file.read())
                project = Project(project_parameters)
                orchestrator_reg_entry = Kratos.Registry[
                    project.GetSettings()["orchestrator"]["name"].GetString()
                ]
                orchestrator_module = importlib.import_module(
                    orchestrator_reg_entry["ModuleName"]
                )
                orchestrator_class = getattr(
                    orchestrator_module, orchestrator_reg_entry["ClassName"]
                )
                orchestrator_instance = orchestrator_class(project)
                orchestrator_instance.Run()

        model = project.GetModel()
        sheet_pile_wall = model.GetModelPart("PorousDomain.Sheet_Pile_Wall")

        if test_helper.want_test_plots():
            self.create_wall_plots(sheet_pile_wall.Nodes, project_path)
            self.create_interface_plots(model.GetModelPart("PorousDomain.Right_Side_Of_Right_Interfaces").Nodes, project_path)

        # Check the obtained results
        reader = GiDOutputFileReader()
        node_ids = [node.Id for node in sheet_pile_wall.Nodes]

        for stage_name in [
            "wall_installation",
            "first_excavation",
            "second_excavation",
            "third_excavation",
        ]:
            base_name = self.stages_info[stage_name]["base_name"]
            csv_filepath = Path(project_path) / f"{base_name}__expected_results_wall.csv"
            expected_results = get_expected_results_from_csv(csv_filepath)

            output_data = reader.read_output_from(
                Path(project_path) / "gid_output" / f"{base_name}.post.res"
            )
            end_time = self.stages_info[stage_name]["end_time"]
            bending_moments = reader.nodal_values_at_time(
                "BENDING_MOMENT", end_time, output_data, node_ids=node_ids
            )
            shear_forces = reader.nodal_values_at_time(
                "SHEAR_FORCE", end_time, output_data, node_ids=node_ids
            )
            horizontal_total_displacements = [
                total_displacement_vector[0]
                for total_displacement_vector in reader.nodal_values_at_time(
                    "TOTAL_DISPLACEMENT", end_time, output_data, node_ids=node_ids
                )
            ]

            for node_id, bending_moment, shear_force, horizontal_total_displacement in zip(
                node_ids, bending_moments, shear_forces, horizontal_total_displacements
            ):
                expected_nodal_results = expected_results[node_id]
                self.assertAlmostEqual(
                    bending_moment,
                    expected_nodal_results[csv_fieldname_bending_moment],
                    msg=f"Bending moment at node {node_id} in stage '{stage_name}'",
                )
                self.assertAlmostEqual(
                    shear_force,
                    expected_nodal_results[csv_fieldname_shear_force],
                    msg=f"Shear force at node {node_id} in stage '{stage_name}'",
                )
                self.assertAlmostEqual(
                    horizontal_total_displacement,
                    expected_nodal_results[csv_fieldname_horizontal_total_displacement],
                    msg=f"Horizontal total displacement at node {node_id} in stage '{stage_name}'",
                )

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
        nodes,
        data_extractor,
    ):
        node_ids = [node.Id for node in nodes]
        y_coords = [node.Y for node in nodes]

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
                    "Kratos",
                    line_style="-",
                    marker=".",
                )
            )
            fem_comparison_csv = (
                Path(project_path) / f"{stage['base_name']}__FE_comparison_interface.csv"
            )
            if fem_comparison_csv.exists():
                fem_comparison_variable = test_helper.get_data_points_from_file(
                    fem_comparison_csv, data_extractor
                )

                data_series_collection.append(
                    plot_utils.DataSeries(
                        fem_comparison_variable,
                        "Commercial FE package",
                        marker="3",
                    )
                )
            variable_data_collections.append(data_series_collection)

        return variable_data_collections

    def create_interface_plots(self, interface_nodes, project_path):
        structural_stages = self.get_plot_stages()

        normal_traction_kratos_label = "GEO_EFFECTIVE_TRACTION_VECTOR"
        normal_traction_plot_label = "Normal traction"
        normal_traction_collections = self.get_variable_collections_per_stage(
            normal_traction_kratos_label,
            normal_traction_plot_label,
            project_path,
            structural_stages,
            interface_nodes,
            extract_normal_traction_and_y_from_line,
        )

        plot_titles = [stage["base_name"] for stage in structural_stages]
        plot_utils.make_sub_plots(
            normal_traction_collections,
            Path(project_path) / "normal_tractions_all_stages.svg",
            titles=plot_titles,
            xlabel=r"Normal Traction [$\mathrm{kN} / \mathrm{m}^2$]",
            ylabel=r"$y$ [m]",
        )

        shear_traction_plot_label = "Shear traction"
        shear_traction_kratos_label = "GEO_EFFECTIVE_TRACTION_VECTOR"
        shear_traction_collections = self.get_variable_collections_per_stage(
            shear_traction_kratos_label,
            shear_traction_plot_label,
            project_path,
            structural_stages,
            interface_nodes,
            extract_shear_traction_and_y_from_line,
        )

        plot_utils.make_sub_plots(
            shear_traction_collections,
            Path(project_path) / "shear_tractions_all_stages.svg",
            titles=plot_titles,
            xlabel=r"Shear Traction [$\mathrm{kN} / \mathrm{m}^2$]",
            ylabel=r"$y$ [m]",
        )

    def get_plot_stages(self):
        return [
            self.stages_info["wall_installation"],
            self.stages_info["first_excavation"],
            self.stages_info["strut_installation"],
            self.stages_info["second_excavation"],
            self.stages_info["third_excavation"],
        ]

    def get_wall_variable_collections_per_stage(
        self,
        kratos_variable_label,
        variable_plot_label,
        project_path,
        plot_stages,
        nodes,
        data_extractor_fem_comparison,
        data_extractor_dsheetpiling=None,
    ):
        node_ids = [node.Id for node in nodes]
        y_coords = [node.Y for node in nodes]

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

            variable_kratos_data = [
                unit_to_k_unit(value) for value in variable_kratos_data
            ]

            data_series_collection = [
                plot_utils.DataSeries(
                    zip(variable_kratos_data, y_coords),
                    "Kratos",
                    line_style="-",
                    marker=".",
                )
            ]

            fem_comparison_csv = (
                Path(project_path) / f"{stage['base_name']}__FE_comparison_wall.csv"
            )
            if fem_comparison_csv.exists():
                data_points = test_helper.get_data_points_from_file(
                    fem_comparison_csv, data_extractor_fem_comparison
                )
                data_series_collection.append(
                    plot_utils.DataSeries(
                        data_points,
                        label="Commercial FE package",
                        marker="3",
                    )
                )

            if (Path(project_path) / f"{stage['base_name']}__DSheetPiling_results.csv").exists() and data_extractor_dsheetpiling:
                comparison_variable = test_helper.get_data_points_from_file(
                    Path(project_path) / f"{stage['base_name']}__DSheetPiling_results.csv",
                    data_extractor_dsheetpiling,
                )
                data_series_collection.append(
                    plot_utils.DataSeries(
                        comparison_variable,
                        "D-Sheet Piling",
                        marker="1",
                    )
                )

            variable_data_collections.append(data_series_collection)

        return variable_data_collections


    def get_bending_moments_data_series_per_stage(self, nodes, project_path):
        return self.get_wall_variable_collections_per_stage("BENDING_MOMENT",
                                                            "Bending moment",
                                                            project_path,
                                                            self.get_plot_stages(),
                                                            nodes,
                                                            extract_bending_moment_and_y_from_line,
                                                            data_extractor_dsheetpiling=extract_bending_moment_and_y_from_line)

    def get_wall_horizontal_total_displacement_collections_per_stage(
        self, project_path, plot_stages, nodes, data_extractor_fem_comparison, data_extractor_dsheetpiling=None
    ):
        node_ids = [node.Id for node in nodes]
        y_coords = [node.Y for node in nodes]

        variable_data_collections = []
        reader = GiDOutputFileReader()

        for stage in plot_stages:
            output_data = reader.read_output_from(
                Path(project_path) / "gid_output" / f"{stage['base_name']}.post.res"
            )

            total_displacement_vectors = reader.nodal_values_at_time(
                "TOTAL_DISPLACEMENT",
                stage["end_time"],
                output_data,
                node_ids=node_ids,
            )

            horizontal_total_displacements = [value[0] for value in total_displacement_vectors]

            data_series_collection = [
                plot_utils.DataSeries(
                    zip(horizontal_total_displacements, y_coords),
                    "Kratos",
                    line_style="-",
                    marker=".",
                )
            ]

            if (
                Path(project_path) / f"{stage['base_name']}__FE_comparison_wall.csv"
            ).exists():
                comparison_variable = test_helper.get_data_points_from_file(
                    Path(project_path) / f"{stage['base_name']}__FE_comparison_wall.csv",
                    data_extractor_fem_comparison,
                )
                data_series_collection.append(
                    plot_utils.DataSeries(
                        comparison_variable,
                        "Commercial FE package",
                        marker="2",
                    )
                )

            if (Path(project_path) / f"{stage['base_name']}__DSheetPiling_results.csv").exists() and data_extractor_dsheetpiling:
                comparison_variable = test_helper.get_data_points_from_file(
                    Path(project_path) / f"{stage['base_name']}__DSheetPiling_results.csv",
                    data_extractor_dsheetpiling,
                )
                data_series_collection.append(
                    plot_utils.DataSeries(
                        comparison_variable,
                        "D-Sheet Piling",
                        marker="1",
                    )
                )

            variable_data_collections.append(data_series_collection)

        return variable_data_collections

    def create_wall_plots(self, nodes, project_path):
        plot_stages = self.get_plot_stages()
        plot_titles = [stage["base_name"] for stage in plot_stages]

        bending_moment_collections = self.get_bending_moments_data_series_per_stage(nodes, project_path)
        plot_utils.make_sub_plots(
            bending_moment_collections,
            Path(project_path) / "bending_moments.svg",
            titles=plot_titles,
            xlabel="Bending moment [kNm/m]",
            ylabel=r"$y$ [m]",
        )

        shear_force_collections = self.get_wall_variable_collections_per_stage(
            "SHEAR_FORCE",
            "Shear force",
            project_path,
            plot_stages,
            nodes,
            extract_shear_force_and_y_from_line,
            data_extractor_dsheetpiling=extract_shear_force_and_y_from_line,
        )
        plot_utils.make_sub_plots(
            shear_force_collections,
            Path(project_path) / "shear_forces.svg",
            titles=plot_titles,
            xlabel="Shear force [kN/m]",
            ylabel=r"$y$ [m]",
        )

        normal_force_collections = self.get_wall_variable_collections_per_stage(
            "AXIAL_FORCE",
            "Normal force",
            project_path,
            plot_stages,
            nodes,
            extract_normal_force_and_y_from_line)
        plot_utils.make_sub_plots(
            normal_force_collections,
            Path(project_path) / "normal_forces.svg",
            titles=plot_titles,
            xlabel="Normal force [kN/m]",
            ylabel=r"$y$ [m]")

        horizontal_total_displacement_collections = (
            self.get_wall_horizontal_total_displacement_collections_per_stage(
                project_path, plot_stages, nodes, extract_horizontal_displacements_from_line, data_extractor_dsheetpiling=extract_horizontal_displacements_from_dsheetpiling_line
            )
        )
        plot_utils.make_sub_plots(
            horizontal_total_displacement_collections,
            Path(project_path) / "horizontal_total_displacements.svg",
            titles=plot_titles,
            xlabel="Horizontal total displacement [m]",
            ylabel=r"$y$ [m]",
        )

    def update_all_expected_results(self):
        print("Updating the expected results...")

        # fmt: off
        self.stages_info = {
            "initial_stage":      {"end_time": -1.0, "base_name": "1_Initial_stage"},
            "null_step":          {"end_time":  0.0, "base_name": "2_Null_step"},
            "wall_installation":  {"end_time":  1.0, "base_name": "3_Wall_installation"},
            "first_excavation":   {"end_time":  2.0, "base_name": "4_First_excavation"},
            "strut_installation": {"end_time":  3.0, "base_name": "5_Anchor_installation"},
            "second_excavation":  {"end_time":  4.0, "base_name": "6_Second_excavation"},
            "third_excavation":   {"end_time":  5.0, "base_name": "7_Third_excavation"},
        }
        # fmt: on

        for case_name in ["linear_elastic", "mohr_coulomb_clay-sand"]:
            self.update_expected_results(case_name)


    def update_expected_results(self, case_name):
        target_dir = Path(
            test_helper.get_file_path(
                Path("crow_validation") / case_name / "staged_construction"
            )
        )

        reader = GiDOutputFileReader()
        node_ids = get_sheetpile_node_ids()

        for stage_name in [
            "wall_installation",
            "first_excavation",
            "second_excavation",
            "third_excavation",
        ]:
            base_name = self.stages_info[stage_name]["base_name"]
            output_data = reader.read_output_from(
                target_dir / "gid_output" / f"{base_name}.post.res"
            )
            end_time = self.stages_info[stage_name]["end_time"]
            bending_moments = reader.nodal_values_at_time(
                "BENDING_MOMENT", end_time, output_data, node_ids=node_ids
            )
            shear_forces = reader.nodal_values_at_time(
                "SHEAR_FORCE", end_time, output_data, node_ids=node_ids
            )
            normal_forces = reader.nodal_values_at_time(
                "AXIAL_FORCE", end_time, output_data, node_ids=node_ids
            )
            horizontal_total_displacements = [
                total_displacement_vector[0]
                for total_displacement_vector in reader.nodal_values_at_time(
                    "TOTAL_DISPLACEMENT", end_time, output_data, node_ids=node_ids
                )
            ]

            with open(
                target_dir / f"{base_name}__expected_results_wall.csv",
                "w",
                newline="",
            ) as csv_file:
                writer = csv.DictWriter(csv_file, fieldnames=csv_fieldnames)

                writer.writeheader()
                for (
                    node_id,
                    bending_moment,
                    shear_force,
                    normal_force,
                    horizontal_total_displacement,
                ) in zip(
                    node_ids, bending_moments, shear_forces, normal_forces, horizontal_total_displacements
                ):
                    writer.writerow(
                        {
                            csv_fieldname_node: node_id,
                            csv_fieldname_bending_moment: bending_moment,
                            csv_fieldname_shear_force: shear_force,
                            csv_fieldname_normal_force: normal_force,
                            csv_fieldname_horizontal_total_displacement: horizontal_total_displacement,
                        }
                    )

    def test_staged_construction_with_linear_elastic_behavior(self):
        self.run_staged_construction_analysis_and_checks("linear_elastic")

    def test_staged_construction_with_mohr_coulomb_clay_sand(self):
        self.run_staged_construction_analysis_and_checks("mohr_coulomb_clay-sand")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--update-expected-results", action="store_true")
    args = parser.parse_args()

    if args.update_expected_results:
        KratosGeoMechanicsCrowValidation().update_all_expected_results()
        sys.exit(0)

    KratosUnittest.main()
