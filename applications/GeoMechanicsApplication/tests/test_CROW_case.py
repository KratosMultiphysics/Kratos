import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.json_utilities as json_utilities
import KratosMultiphysics.GeoMechanicsApplication.context_managers as context_managers
from KratosMultiphysics.GeoMechanicsApplication.unit_conversions import unit_to_k_unit
import test_helper
import helper_utilities

import argparse
import csv
from pathlib import Path
import sys

if test_helper.want_test_plots():
    import KratosMultiphysics.GeoMechanicsApplication.geo_plot_utilities as plot_utils


linear_elastic_dir_name = "linear_elastic"
mohr_coulomb_clay_sand_dir_name = "mohr_coulomb_clay-sand"
staged_construction = "staged_construction"

wall_output_postfix = "output_wall"
interface_output_postfix = "output_interface"

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

stages_to_be_checked = [
    "3_Wall_installation",
    "4_First_excavation",
    "6_Second_excavation",
    "7_Third_excavation",
]
stages_to_be_plotted = [
    "3_Wall_installation",
    "4_First_excavation",
    "5_Anchor_installation",
    "6_Second_excavation",
    "7_Third_excavation",
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


def max_abs_nodal_value(values_by_node_id_and_result_item, node_ids, result_item_name):
    return max(
        [
            abs(values_by_node_id_and_result_item[node_id][result_item_name])
            for node_id in node_ids
        ]
    )


def set_quasi_newton_method(project_parameters, quasi_newton_method):
    # Add/set quasi-Newton method settings that apply to all stages
    for stage in project_parameters["stages"].values():
        solver_settings = stage["stage_settings"]["solver_settings"]
        solver_settings["strategy_type"].SetString("quasi_newton_raphson")
        solver_settings.AddString("quasi_newton_type", quasi_newton_method)
        solver_settings.AddInt("quasi_newton_raphson_restart_interval", 100)

    # Add stage-dependent quasi-Newton method settings
    stage_dependent_settings = {
        "1_Initial_stage": {
            "quasi_newton_raphson_max_rank": 10,
        },
        "2_Null_step": {
            "quasi_newton_raphson_max_rank": 10,
        },
        "3_Wall_installation": {
            "quasi_newton_raphson_max_rank": 10,
        },
        "4_First_excavation": {
            "quasi_newton_raphson_max_rank": 10,
        },
        "5_Anchor_installation": {
            "quasi_newton_raphson_max_rank": 20,
        },
        "6_Second_excavation": {
            "quasi_newton_raphson_max_rank": 10,
        },
        "7_Third_excavation": {
            "quasi_newton_raphson_max_rank": 20,
        },
    }
    for stage_name, settings in stage_dependent_settings.items():
        solver_settings = project_parameters["stages"][stage_name]["stage_settings"][
            "solver_settings"
        ]
        for key, value in settings.items():
            solver_settings.AddInt(key, value)


class KratosGeoMechanicsCrowValidation(KratosUnittest.TestCase):
    def setUp(self):
        super().setUp()

        self.common_test_files_dir = (
            Path(test_helper.get_file_path("crow_validation")) / "common"
        )

        # The following attributes will be populated by the specific simulation runs
        self.analysis_type = None
        self.test_path = None
        self.csv_files_dir = None
        self.modify_project_parameters = None

    def prepare_test_run(
        self,
        material_model_dir_name,
        analysis_type,
        variant,
        modify_project_parameters=None,
        run_analysis=helper_utilities.run_orchestrator,
    ):
        self.analysis_type = analysis_type

        base_test_path = Path(
            test_helper.get_file_path(
                Path("crow_validation") / material_model_dir_name / analysis_type
            )
        )
        self.test_path = base_test_path / variant
        self.test_path.mkdir(exist_ok=True)
        self.csv_files_dir = base_test_path
        self.modify_project_parameters = modify_project_parameters
        self.run_analysis = run_analysis

    def run_simulation_and_checks(self):
        with context_managers.set_cwd_to(self.test_path):
            with open(
                self.common_test_files_dir / f"{self.analysis_type}.json", "r"
            ) as analysis_file:
                project_parameters = Kratos.Parameters(analysis_file.read())

            if self.modify_project_parameters is not None:
                self.modify_project_parameters(project_parameters)

            project = self.run_analysis(project_parameters)

        model = project.GetModel()
        sheet_pile_wall = model.GetModelPart("PorousDomain.Sheet_Pile_Wall")

        if test_helper.want_test_plots():
            self.create_wall_plots(sheet_pile_wall.Nodes)
            self.create_interface_plots(
                model.GetModelPart("PorousDomain.Right_Side_Of_Right_Interfaces").Nodes
            )

        # Check the obtained results
        node_ids = [node.Id for node in sheet_pile_wall.Nodes]

        for stage_name in stages_to_be_checked:
            json_output = json_utilities.read_external_json(
                self.file_path_to_json_output(stage_name, wall_output_postfix)
            )

            csv_filepath = (
                self.csv_files_dir / f"{stage_name}__expected_results_wall.csv"
            )
            expected_results = get_expected_results_from_csv(csv_filepath)

            expected_max_abs_bending_moment = max_abs_nodal_value(
                expected_results, node_ids, csv_fieldname_bending_moment
            )
            expected_max_abs_shear_force = max_abs_nodal_value(
                expected_results, node_ids, csv_fieldname_shear_force
            )
            expected_max_abs_horizontal_total_displacement = max_abs_nodal_value(
                expected_results, node_ids, csv_fieldname_horizontal_total_displacement
            )

            relative_tolerance = 0.02

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
                        expected_bending_moment,
                        absolute_tolerance=relative_tolerance
                        * expected_max_abs_bending_moment,
                        relative_tolerance=relative_tolerance,
                    ),
                    msg=f"Bending moment at node {node_id} in stage '{stage_name}'",
                )

                expected_shear_force = expected_nodal_results[csv_fieldname_shear_force]
                self.assertAlmostEqual(
                    shear_force,
                    expected_shear_force,
                    places=None,
                    delta=test_helper.calculate_delta(
                        expected_shear_force,
                        absolute_tolerance=relative_tolerance
                        * expected_max_abs_shear_force,
                        relative_tolerance=relative_tolerance,
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
                        absolute_tolerance=relative_tolerance
                        * expected_max_abs_horizontal_total_displacement,
                        relative_tolerance=relative_tolerance,
                    ),
                    msg=f"Horizontal total displacement at node {node_id} in stage '{stage_name}'",
                )

    def file_path_to_json_output(self, stage_name, postfix):
        return self.test_path / f"{stage_name}__{postfix}.json"

    def get_variable_collections_per_stage(
        self,
        kratos_variable_label,
        stage_names,
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
        for stage_name in stage_names:
            json_data = json_utilities.read_external_json(
                self.file_path_to_json_output(stage_name, postfix_json_output)
            )
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

            # Get results of comparison FE analysis (if they exist)
            csv_file_path = (
                self.csv_files_dir / f"{stage_name}__{postfix_fem_comparison_csv}.csv"
            )
            if csv_file_path.exists():
                data_points = test_helper.get_data_points_from_file(
                    csv_file_path, data_extractor_fem_comparison
                )
                data_series_collection.append(
                    plot_utils.DataSeries(
                        data_points,
                        label="Commercial FE package",
                        marker="3",
                    )
                )

            # Get results of D-Sheet Piling analysis (if they exist)
            csv_file_path = (
                self.csv_files_dir / f"{stage_name}__DSheetPiling_results.csv"
            )
            if csv_file_path.exists() and data_extractor_dsheetpiling:
                data_points = test_helper.get_data_points_from_file(
                    csv_file_path,
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
        effective_traction_vector_label = "GEO_EFFECTIVE_TRACTION_VECTOR"
        postfix_fem_comparison_csv = "FE_comparison_interface"
        to_normal_traction_in_kPa = lambda traction_vector: unit_to_k_unit(
            traction_vector[0]
        )
        normal_traction_collections = self.get_variable_collections_per_stage(
            effective_traction_vector_label,
            stages_to_be_plotted,
            interface_nodes,
            interface_output_postfix,
            to_normal_traction_in_kPa,
            postfix_fem_comparison_csv,
            extract_normal_traction_and_y_from_line,
        )

        plot_utils.make_sub_plots(
            normal_traction_collections,
            self.test_path / "normal_tractions_all_stages.svg",
            titles=stages_to_be_plotted,
            xlabel=r"Normal Traction [$\mathrm{kN} / \mathrm{m}^2$]",
            ylabel=r"$y$ [m]",
        )

        to_shear_traction_in_kPa = lambda traction_vector: unit_to_k_unit(
            traction_vector[1]
        )
        shear_traction_collections = self.get_variable_collections_per_stage(
            effective_traction_vector_label,
            stages_to_be_plotted,
            interface_nodes,
            interface_output_postfix,
            to_shear_traction_in_kPa,
            postfix_fem_comparison_csv,
            extract_shear_traction_and_y_from_line,
        )

        plot_utils.make_sub_plots(
            shear_traction_collections,
            self.test_path / "shear_tractions_all_stages.svg",
            titles=stages_to_be_plotted,
            xlabel=r"Shear Traction [$\mathrm{kN} / \mathrm{m}^2$]",
            ylabel=r"$y$ [m]",
        )

    def get_bending_moment_data_series_per_stage(
        self, nodes, postfix_fem_comparison_csv
    ):
        return self.get_variable_collections_per_stage(
            "BENDING_MOMENT",
            stages_to_be_plotted,
            nodes,
            wall_output_postfix,
            unit_to_k_unit,
            postfix_fem_comparison_csv,
            extract_bending_moment_and_y_from_line,
            data_extractor_dsheetpiling=extract_bending_moment_and_y_from_line,
        )

    def get_shear_force_data_series_per_stage(self, nodes, postfix_fem_comparison_csv):
        return self.get_variable_collections_per_stage(
            "SHEAR_FORCE",
            stages_to_be_plotted,
            nodes,
            wall_output_postfix,
            unit_to_k_unit,
            postfix_fem_comparison_csv,
            extract_shear_force_and_y_from_line,
            data_extractor_dsheetpiling=extract_shear_force_and_y_from_line,
        )

    def get_normal_force_data_series_per_stage(self, nodes, postfix_fem_comparison_csv):
        return self.get_variable_collections_per_stage(
            "AXIAL_FORCE",
            stages_to_be_plotted,
            nodes,
            wall_output_postfix,
            unit_to_k_unit,
            postfix_fem_comparison_csv,
            extract_normal_force_and_y_from_line,
        )

    def get_horizontal_total_displacement_data_series_per_stage(
        self, nodes, postfix_fem_comparison_csv
    ):
        no_transformation = lambda value: value
        return self.get_variable_collections_per_stage(
            "TOTAL_DISPLACEMENT_X",
            stages_to_be_plotted,
            nodes,
            wall_output_postfix,
            no_transformation,
            postfix_fem_comparison_csv,
            extract_horizontal_displacements_from_line,
            data_extractor_dsheetpiling=extract_horizontal_displacements_from_dsheetpiling_line,
        )

    def create_wall_plots(self, nodes):
        postfix_fem_comparison_csv = "FE_comparison_wall"

        bending_moment_collections = self.get_bending_moment_data_series_per_stage(
            nodes, postfix_fem_comparison_csv
        )
        plot_utils.make_sub_plots(
            bending_moment_collections,
            self.test_path / "bending_moments.svg",
            titles=stages_to_be_plotted,
            xlabel="Bending moment [kNm/m]",
            ylabel=r"$y$ [m]",
        )

        shear_force_collections = self.get_shear_force_data_series_per_stage(
            nodes, postfix_fem_comparison_csv
        )
        plot_utils.make_sub_plots(
            shear_force_collections,
            self.test_path / "shear_forces.svg",
            titles=stages_to_be_plotted,
            xlabel="Shear force [kN/m]",
            ylabel=r"$y$ [m]",
        )

        normal_force_collections = self.get_normal_force_data_series_per_stage(
            nodes, postfix_fem_comparison_csv
        )
        plot_utils.make_sub_plots(
            normal_force_collections,
            self.test_path / "normal_forces.svg",
            titles=stages_to_be_plotted,
            xlabel="Normal force [kN/m]",
            ylabel=r"$y$ [m]",
        )

        horizontal_total_displacement_collections = (
            self.get_horizontal_total_displacement_data_series_per_stage(
                nodes, postfix_fem_comparison_csv
            )
        )
        plot_utils.make_sub_plots(
            horizontal_total_displacement_collections,
            self.test_path / "horizontal_total_displacements.svg",
            titles=stages_to_be_plotted,
            xlabel="Horizontal total displacement [m]",
            ylabel=r"$y$ [m]",
        )

    def update_all_expected_results(self):
        print("Updating the expected results...")

        self.update_expected_results(linear_elastic_dir_name)
        self.update_expected_results(
            mohr_coulomb_clay_sand_dir_name, analysis_type="staged_construction_broyden"
        )
        self.update_expected_results(
            mohr_coulomb_clay_sand_dir_name, analysis_type="staged_construction_lbfgs"
        )

    def update_expected_results(self, case_name, analysis_type=staged_construction):
        self.analysis_type = analysis_type
        self.test_path = Path(
            test_helper.get_file_path(
                Path("crow_validation") / case_name / self.analysis_type
            )
        )

        mdpa_file_path_without_file_extension = test_helper.get_file_path(
            self.common_test_files_dir / "model"
        )
        model = Kratos.Model()
        main_model_part = model.CreateModelPart("PorousDomain")
        Kratos.ModelPartIO(mdpa_file_path_without_file_extension).ReadModelPart(
            main_model_part
        )
        node_ids = [
            node.Id for node in model.GetModelPart("PorousDomain.Sheet_Pile_Wall").Nodes
        ]

        for stage_name in stages_to_be_checked:
            json_output = json_utilities.read_external_json(
                self.file_path_to_json_output(stage_name, wall_output_postfix)
            )

            with open(
                self.test_path / f"{stage_name}__expected_results_wall.csv",
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
        self.prepare_test_run(
            material_model_dir_name=linear_elastic_dir_name,
            analysis_type=staged_construction,
            variant="as-is",
        )
        self.run_simulation_and_checks()

    def test_staged_construction_with_mohr_coulomb_clay_sand_broyden(self):
        self.prepare_test_run(
            material_model_dir_name=mohr_coulomb_clay_sand_dir_name,
            analysis_type=staged_construction,
            variant="broyden",
            modify_project_parameters=(
                lambda project_parameters: set_quasi_newton_method(
                    project_parameters, "broyden"
                )
            ),
        )
        self.run_simulation_and_checks()

    def test_staged_construction_with_mohr_coulomb_clay_sand_lbfgs(self):
        self.prepare_test_run(
            material_model_dir_name=mohr_coulomb_clay_sand_dir_name,
            analysis_type=staged_construction,
            variant="lbfgs",
            modify_project_parameters=(
                lambda project_parameters: set_quasi_newton_method(
                    project_parameters, "lbfgs"
                )
            ),
        )
        self.run_simulation_and_checks()

    def test_staged_construction_with_mohr_coulomb_clay_sand_linear_iteration(self):
        self.prepare_test_run(
            material_model_dir_name=mohr_coulomb_clay_sand_dir_name,
            analysis_type=staged_construction,
            variant="linear_iteration",
        )
        self.run_simulation_and_checks()

    def test_staged_construction_with_mohr_coulomb_clay_sand_using_save_and_load(
        self,
    ):
        self.prepare_test_run(
            material_model_dir_name=mohr_coulomb_clay_sand_dir_name,
            analysis_type=staged_construction,
            variant="with_save_and_load",
            run_analysis=(
                helper_utilities.run_multistage_analysis_with_intermediate_save_and_load
            ),
        )
        self.run_simulation_and_checks()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--update-expected-results", action="store_true")
    args = parser.parse_args()

    if args.update_expected_results:
        KratosGeoMechanicsCrowValidation().update_all_expected_results()
        sys.exit(0)

    KratosUnittest.main()
