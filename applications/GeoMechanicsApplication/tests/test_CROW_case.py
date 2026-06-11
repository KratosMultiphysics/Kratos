import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.GeoMechanicsApplication.context_managers as context_managers
from KratosMultiphysics.GeoMechanicsApplication.unit_conversions import unit_to_k_unit
import test_helper

import argparse
import csv
import os
import json
from pathlib import Path
import sys
from helper_utilities import run_orchestrator

if test_helper.want_test_plots():
    import KratosMultiphysics.GeoMechanicsApplication.geo_plot_utilities as plot_utils


wall_output_postfix = "output_wall"

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

        # The following attributes will be populated by the specific simulation runs
        self.stages_info = None
        self.analysis_type = None
        self.test_path = None

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
        self.analysis_type = "staged_construction"
        self.test_path = Path(
            test_helper.get_file_path(
                Path("crow_validation") / material_model_dir_name / self.analysis_type
            )
        )

        self.run_simulation_and_checks()

    def run_simulation_and_checks(self):
        with context_managers.set_cwd_to(self.test_path):
            with open(
                Path("..") / ".." / "common" / f"{self.analysis_type}.json", "r"
            ) as analysis_file:
                project = run_orchestrator(Kratos.Parameters(analysis_file.read()))

        model = project.GetModel()
        sheet_pile_wall = model.GetModelPart("PorousDomain.Sheet_Pile_Wall")

        if test_helper.want_test_plots():
            self.create_wall_plots(sheet_pile_wall.Nodes)
            self.create_interface_plots(
                model.GetModelPart("PorousDomain.Right_Side_Of_Right_Interfaces").Nodes
            )

        # Check the obtained results
        node_ids = [node.Id for node in sheet_pile_wall.Nodes]

        for stage_name in [
            "wall_installation",
            "first_excavation",
            "second_excavation",
            "third_excavation",
        ]:
            json_output = self.read_json_output(
                self.stages_info[stage_name], postfix=wall_output_postfix
            )

            base_name = self.stages_info[stage_name]["base_name"]
            csv_filepath = self.test_path / f"{base_name}__expected_results_wall.csv"
            expected_results = get_expected_results_from_csv(csv_filepath)

            relative_tolerance = (
                100.0 * test_helper.default_relative_tolerance_for_assertions
            )

            for (
                node_id,
                bending_moment,
                shear_force,
                horizontal_total_displacement,
            ) in zip(
                node_ids,
                test_helper.get_bending_moments_from_json_output(json_output, node_ids),
                test_helper.get_shear_forces_from_json_output(json_output, node_ids),
                test_helper.get_total_displacement_x_from_json_output(
                    json_output, node_ids
                ),
            ):
                expected_nodal_results = expected_results[node_id]

                expected_bending_moment = expected_nodal_results[
                    csv_fieldname_bending_moment
                ]
                self.assertAlmostEqual(
                    bending_moment,
                    expected_bending_moment,
                    places=None,
                    delta=test_helper.calculate_delta(
                        expected_bending_moment, relative_tolerance=relative_tolerance
                    ),
                    msg=f"Bending moment at node {node_id} in stage '{stage_name}'",
                )

                expected_shear_force = expected_nodal_results[csv_fieldname_shear_force]
                self.assertAlmostEqual(
                    shear_force,
                    expected_shear_force,
                    places=None,
                    delta=test_helper.calculate_delta(
                        expected_shear_force, relative_tolerance=relative_tolerance
                    ),
                    msg=f"Shear force at node {node_id} in stage '{stage_name}'",
                )

                expected_horizontal_total_displacement = expected_nodal_results[
                    csv_fieldname_horizontal_total_displacement
                ]
                self.assertAlmostEqual(
                    horizontal_total_displacement,
                    expected_horizontal_total_displacement,
                    places=None,
                    delta=test_helper.calculate_delta(
                        expected_horizontal_total_displacement,
                        relative_tolerance=relative_tolerance,
                    ),
                    msg=f"Horizontal total displacement at node {node_id} in stage '{stage_name}'",
                )

    def read_json_output(self, stage, postfix):
        with open(
            os.path.join(self.test_path, f"{stage['base_name']}__{postfix}.json"),
            "r",
        ) as output_file:
            return json.load(output_file)

    def get_variable_collections_per_stage(
        self,
        kratos_variable_label,
        plot_stages,
        nodes,
        postfix_json_output,
        transform_output,
        postfix_fem_comparison_csv,
        data_extractor_fem_comparison,
        data_extractor_dsheetpiling=None,
    ):
        node_ids = [node.Id for node in nodes]
        y_coords = [node.Y for node in nodes]

        variable_data_collections = []
        for stage in plot_stages:
            json_data = self.read_json_output(stage, postfix_json_output)
            variable_kratos_data = test_helper.get_nodal_values_from_json_output(
                json_data, kratos_variable_label, node_ids
            )

            variable_kratos_data = [
                transform_output(value) for value in variable_kratos_data
            ]
            sorted_y, sorted_data = zip(*sorted(zip(y_coords, variable_kratos_data)))
            data_series_collection = [
                plot_utils.DataSeries(
                    zip(sorted_data, sorted_y),
                    "Kratos",
                    line_style="-",
                    marker=".",
                )
            ]

            fem_comparison_csv = (
                self.test_path
                / f"{stage['base_name']}__{postfix_fem_comparison_csv}.csv"
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

            if (
                self.test_path / f"{stage['base_name']}__DSheetPiling_results.csv"
            ).exists() and data_extractor_dsheetpiling:
                data_points = test_helper.get_data_points_from_file(
                    self.test_path / f"{stage['base_name']}__DSheetPiling_results.csv",
                    data_extractor_dsheetpiling,
                )
                data_series_collection.append(
                    plot_utils.DataSeries(
                        data_points,
                        "D-Sheet Piling",
                        marker="1",
                    )
                )

            variable_data_collections.append(data_series_collection)

        return variable_data_collections

    def create_interface_plots(self, interface_nodes):
        structural_stages = self.get_plot_stages()

        effective_traction_vector_label = "GEO_EFFECTIVE_TRACTION_VECTOR"
        postfix_fem_comparison_csv = "FE_comparison_interface"
        postfix_json_output = "output_interface"
        to_normal_traction_in_kPa = lambda traction_vector: unit_to_k_unit(
            traction_vector[0]
        )
        normal_traction_collections = self.get_variable_collections_per_stage(
            effective_traction_vector_label,
            structural_stages,
            interface_nodes,
            postfix_json_output,
            to_normal_traction_in_kPa,
            postfix_fem_comparison_csv,
            extract_normal_traction_and_y_from_line,
        )

        plot_titles = [stage["base_name"] for stage in structural_stages]
        plot_utils.make_sub_plots(
            normal_traction_collections,
            self.test_path / "normal_tractions_all_stages.svg",
            titles=plot_titles,
            xlabel=r"Normal Traction [$\mathrm{kN} / \mathrm{m}^2$]",
            ylabel=r"$y$ [m]",
        )

        to_shear_traction_in_kPa = lambda traction_vector: unit_to_k_unit(
            traction_vector[1]
        )
        shear_traction_collections = self.get_variable_collections_per_stage(
            effective_traction_vector_label,
            structural_stages,
            interface_nodes,
            postfix_json_output,
            to_shear_traction_in_kPa,
            postfix_fem_comparison_csv,
            extract_shear_traction_and_y_from_line,
        )

        plot_utils.make_sub_plots(
            shear_traction_collections,
            self.test_path / "shear_tractions_all_stages.svg",
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

    def get_bending_moments_data_series_per_stage(self, nodes):
        return self.get_variable_collections_per_stage(
            "BENDING_MOMENT",
            self.get_plot_stages(),
            nodes,
            wall_output_postfix,
            unit_to_k_unit,
            "FE_comparison_wall",
            extract_bending_moment_and_y_from_line,
            data_extractor_dsheetpiling=extract_bending_moment_and_y_from_line,
        )

    def create_wall_plots(self, nodes):
        plot_stages = self.get_plot_stages()
        plot_titles = [stage["base_name"] for stage in plot_stages]

        bending_moment_collections = self.get_bending_moments_data_series_per_stage(
            nodes
        )
        plot_utils.make_sub_plots(
            bending_moment_collections,
            self.test_path / "bending_moments.svg",
            titles=plot_titles,
            xlabel="Bending moment [kNm/m]",
            ylabel=r"$y$ [m]",
        )

        postfix_fem_comparison_csv = "FE_comparison_wall"
        postfix_json_output = wall_output_postfix
        shear_force_collections = self.get_variable_collections_per_stage(
            "SHEAR_FORCE",
            plot_stages,
            nodes,
            postfix_json_output,
            unit_to_k_unit,
            postfix_fem_comparison_csv,
            extract_shear_force_and_y_from_line,
            data_extractor_dsheetpiling=extract_shear_force_and_y_from_line,
        )
        plot_utils.make_sub_plots(
            shear_force_collections,
            self.test_path / "shear_forces.svg",
            titles=plot_titles,
            xlabel="Shear force [kN/m]",
            ylabel=r"$y$ [m]",
        )

        normal_force_collections = self.get_variable_collections_per_stage(
            "AXIAL_FORCE",
            plot_stages,
            nodes,
            postfix_json_output,
            unit_to_k_unit,
            postfix_fem_comparison_csv,
            extract_normal_force_and_y_from_line,
        )
        plot_utils.make_sub_plots(
            normal_force_collections,
            self.test_path / "normal_forces.svg",
            titles=plot_titles,
            xlabel="Normal force [kN/m]",
            ylabel=r"$y$ [m]",
        )

        no_transformation = lambda value: value
        horizontal_total_displacement_collections = self.get_variable_collections_per_stage(
            "TOTAL_DISPLACEMENT_X",
            plot_stages,
            nodes,
            postfix_json_output,
            no_transformation,
            postfix_fem_comparison_csv,
            extract_horizontal_displacements_from_line,
            data_extractor_dsheetpiling=extract_horizontal_displacements_from_dsheetpiling_line,
        )
        plot_utils.make_sub_plots(
            horizontal_total_displacement_collections,
            self.test_path / "horizontal_total_displacements.svg",
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
        self.analysis_type = "staged_construction"
        self.test_path = Path(
            test_helper.get_file_path(
                Path("crow_validation") / case_name / self.analysis_type
            )
        )

        mdpa_file_path_without_file_extension = test_helper.get_file_path(
            Path("crow_validation") / "common" / "model"
        )
        model = Kratos.Model()
        main_model_part = model.CreateModelPart("PorousDomain")
        Kratos.ModelPartIO(mdpa_file_path_without_file_extension).ReadModelPart(
            main_model_part
        )
        node_ids = [
            node.Id for node in model.GetModelPart("PorousDomain.Sheet_Pile_Wall").Nodes
        ]

        for stage_name in [
            "wall_installation",
            "first_excavation",
            "second_excavation",
            "third_excavation",
        ]:
            json_output = self.read_json_output(
                self.stages_info[stage_name], postfix=wall_output_postfix
            )

            base_name = self.stages_info[stage_name]["base_name"]
            with open(
                self.test_path / f"{base_name}__expected_results_wall.csv",
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
                    node_ids,
                    test_helper.get_bending_moments_from_json_output(
                        json_output, node_ids
                    ),
                    test_helper.get_shear_forces_from_json_output(
                        json_output, node_ids
                    ),
                    test_helper.get_normal_forces_from_json_output(
                        json_output, node_ids
                    ),
                    test_helper.get_total_displacement_x_from_json_output(
                        json_output, node_ids
                    ),
                ):
                    writer.writerow(
                        {
                            csv_fieldname_node: node_id,
                            csv_fieldname_bending_moment: f"{bending_moment:.6}",
                            csv_fieldname_shear_force: f"{shear_force:.6}",
                            csv_fieldname_normal_force: f"{normal_force:.6}",
                            csv_fieldname_horizontal_total_displacement: f"{horizontal_total_displacement:.6}",
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
