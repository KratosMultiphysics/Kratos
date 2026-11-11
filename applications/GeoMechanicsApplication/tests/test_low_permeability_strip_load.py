import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.json_utilities as json_utilities
import KratosMultiphysics.GeoMechanicsApplication.context_managers as context_managers
import test_helper
import helper_utilities

import csv
from pathlib import Path


stage_names = ["1_Initial", "2_Load"]

benchmark_dir = Path(test_helper.get_file_path("low_permeability_strip_load"))
common_test_files_dir = benchmark_dir / "common"

# The hydrostatic water pressure process uses this specific weight and a reference coordinate of
# y = 0.0 (the top surface), such that the hydrostatic water pressure equals specific_weight * y.
specific_weight_of_water = 10000.0

csv_fieldname_node = "node"
csv_fieldname_total_displacement_x = "total_displacement_x_in_m"
csv_fieldname_total_displacement_y = "total_displacement_y_in_m"

absolute_tolerance_water_pressure = 1.0  # [Pa], on water pressures up to 400 kPa
absolute_tolerance_displacement = 1.0e-6  # [m], on settlements up to 0.28 m
relative_tolerance_displacement = 1.0e-4


def solver_settings_of_all_stages(project_parameters):
    return [
        stage["stage_settings"]["solver_settings"]
        for stage in project_parameters["stages"].values()
    ]


def set_linear_solver_scaling(project_parameters, use_scaling):
    for solver_settings in solver_settings_of_all_stages(project_parameters):
        solver_settings["linear_solver_settings"]["scaling"].SetBool(use_scaling)


def set_quasi_newton_method_with_scaling(project_parameters, quasi_newton_method):
    for solver_settings in solver_settings_of_all_stages(project_parameters):
        solver_settings["strategy_type"].SetString("quasi_newton")
        solver_settings.AddString("quasi_newton_type", quasi_newton_method)
        solver_settings.AddInt("quasi_newton_restart_interval", 100)
        solver_settings.AddInt("quasi_newton_max_rank", 10)
        solver_settings.AddBool("quasi_newton_scaling", True)


def fix_hydrostatic_water_pressure(project_parameters):
    # Prescribe the (known) hydrostatic water pressure at all nodes, such that only the force
    # equations remain to be solved. This gives the drained reference solution.
    for stage in project_parameters["stages"].values():
        constraints = stage["stage_settings"]["processes"]["constraints_process_list"]
        for i in range(constraints.size()):
            parameters = constraints[i]["Parameters"]
            if (
                parameters["variable_name"].GetString() == "WATER_PRESSURE"
                and parameters.Has("fluid_pressure_type")
                and parameters["fluid_pressure_type"].GetString() == "Hydrostatic"
            ):
                parameters["is_fixed"].SetBool(True)


def get_expected_displacements_from_csv(csv_filepath):
    with open(csv_filepath, newline="") as csv_file:
        return {
            int(row[csv_fieldname_node]): (
                float(row[csv_fieldname_total_displacement_x]),
                float(row[csv_fieldname_total_displacement_y]),
            )
            for row in csv.DictReader(csv_file)
        }


class KratosGeoMechanicsLowPermeabilityStripLoadTests(KratosUnittest.TestCase):
    """
    Strip load on a layered soil, modelled with fully coupled U-Pw elements with a very low
    permeability (1e-16 m^2) and solved with the static scheme. The entries of the flow equations
    are many orders of magnitude smaller than those of the force equations, such that the linear
    system can only be solved accurately when it is scaled.

    Since the top surface is drained and all other boundaries are impermeable, the static flow
    problem has a hydrostatic solution. Therefore, the water pressures must be hydrostatic, and the
    displacements must equal those of the drained reference, in which the water pressures are
    prescribed. The expected displacements have been obtained from that drained reference.
    """

    def run_simulation_and_checks(self, variant, modify_project_parameters):
        test_path = benchmark_dir / variant
        test_path.mkdir(exist_ok=True)

        with context_managers.set_cwd_to(test_path):
            with open(common_test_files_dir / "strip_load.json", "r") as project_parameters_file:
                project_parameters = Kratos.Parameters(project_parameters_file.read())

            modify_project_parameters(project_parameters)
            project = helper_utilities.run_orchestrator(project_parameters)

        nodes = project.GetModel().GetModelPart("PorousDomain").Nodes
        node_ids = [node.Id for node in nodes]
        y_coordinates = [node.Y0 for node in nodes]

        for stage_name in stage_names:
            json_output = json_utilities.read_external_json(
                str(test_path / f"{stage_name}__output.json")
            )
            self.check_water_pressures_are_hydrostatic(
                json_output, node_ids, y_coordinates, stage_name
            )
            self.check_displacements(json_output, node_ids, stage_name)

    def check_water_pressures_are_hydrostatic(
        self, json_output, node_ids, y_coordinates, stage_name
    ):
        water_pressures = test_helper.get_nodal_values_from_json_output(
            json_output, "WATER_PRESSURE", node_ids, index=-1
        )
        for node_id, y, water_pressure in zip(node_ids, y_coordinates, water_pressures):
            self.assertAlmostEqual(
                water_pressure,
                specific_weight_of_water * y,
                places=None,
                delta=absolute_tolerance_water_pressure,
                msg=f"Water pressure at node {node_id} in stage '{stage_name}'",
            )

    def check_displacements(self, json_output, node_ids, stage_name):
        expected_displacements = get_expected_displacements_from_csv(
            benchmark_dir / f"{stage_name}__expected_results.csv"
        )
        displacements_x = test_helper.get_nodal_values_from_json_output(
            json_output, "TOTAL_DISPLACEMENT_X", node_ids, index=-1
        )
        displacements_y = test_helper.get_nodal_values_from_json_output(
            json_output, "TOTAL_DISPLACEMENT_Y", node_ids, index=-1
        )
        for node_id, displacement_x, displacement_y in zip(
            node_ids, displacements_x, displacements_y
        ):
            expected_x, expected_y = expected_displacements[node_id]
            for direction, actual, expected in (
                ("x", displacement_x, expected_x),
                ("y", displacement_y, expected_y),
            ):
                self.assertAlmostEqual(
                    actual,
                    expected,
                    places=None,
                    delta=test_helper.calculate_delta(
                        expected,
                        absolute_tolerance=absolute_tolerance_displacement,
                        relative_tolerance=relative_tolerance_displacement,
                    ),
                    msg=f"Total displacement {direction} at node {node_id} in stage '{stage_name}'",
                )

    def test_drained_reference(self):
        # With prescribed water pressures, the result does not depend on scaling
        def modify(project_parameters):
            fix_hydrostatic_water_pressure(project_parameters)
            set_linear_solver_scaling(project_parameters, False)

        self.run_simulation_and_checks("drained_reference", modify)

    def test_newton_raphson_with_scaling(self):
        self.run_simulation_and_checks(
            "newton_raphson_with_scaling",
            lambda project_parameters: set_linear_solver_scaling(project_parameters, True),
        )

    def test_quasi_newton_lbfgs_with_scaling(self):
        self.run_simulation_and_checks(
            "quasi_newton_lbfgs_with_scaling",
            lambda project_parameters: set_quasi_newton_method_with_scaling(
                project_parameters, "lbfgs"
            ),
        )

    def test_quasi_newton_broyden_with_scaling(self):
        self.run_simulation_and_checks(
            "quasi_newton_broyden_with_scaling",
            lambda project_parameters: set_quasi_newton_method_with_scaling(
                project_parameters, "broyden"
            ),
        )


if __name__ == "__main__":
    KratosUnittest.main()
