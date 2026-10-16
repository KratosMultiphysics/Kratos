import json
import os

import KratosMultiphysics as KM
import KratosMultiphysics.DamApplication as KratosDam
import KratosMultiphysics.PoromechanicsApplication as KratosPoro
import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.DamApplication.dam_analysis import DamAnalysis


def _load_construction_parameters():
    """Loads the old-style construction ProjectParameters file (an existing/older
    user ProjectParameters file."""
    parameters_file = os.path.join(
        os.path.dirname(__file__), "construction", "construction_parameters.json")
    with open(parameters_file) as f:
        return json.load(f)


def _run_construction():
    """Runs the construction analysis and returns the resulting Model."""
    parameters = _load_construction_parameters()
    model = KM.Model()
    simulation = DamAnalysis(model, KM.Parameters(json.dumps(parameters)))
    simulation.Run()
    return model


def _mechanical_solver_settings():
    """An old-style mechanical solver settings dict."""
    return KM.Parameters(json.dumps({
        "solver_type": "dam_mechanical_solver",
        "model_import_settings": {"input_type": "mdpa", "input_filename": "x", "input_file_label": 0},
        "echo_level": 0,
        "buffer_size": 2,
        "processes_sub_model_part_list": [""],
        "mechanical_solver_settings": {
            "echo_level": 0,
            "reform_dofs_at_each_step": False,
            "clear_storage": False,
            "compute_reactions": False,
            "move_mesh_flag": True,
            "solution_type": "Quasi-Static",
            "scheme_type": "Newmark",
            "rayleigh_m": 0.0,
            "rayleigh_k": 0.0,
            "strategy_type": "Newton-Raphson",
            "convergence_criterion": "Displacement_criterion",
            "displacement_relative_tolerance": 1.0e-4,
            "displacement_absolute_tolerance": 1.0e-9,
            "residual_relative_tolerance": 1.0e-4,
            "residual_absolute_tolerance": 1.0e-9,
            "max_iteration": 15,
            "desired_iterations": 4,
            "max_radius_factor": 20.0,
            "min_radius_factor": 0.5,
            "block_builder": True,
            "nonlocal_damage": False,
            "characteristic_length": 0.05,
            "search_neighbours_step": False,
            "linear_solver_settings": {
                "solver_type": "bicgstab",
                "max_iteration": 200,
                "tolerance": 1.0e-7,
                "preconditioner_type": "ilu0"
            },
            "problem_domain_sub_model_part_list": [""],
            "body_domain_sub_model_part_list": [],
            "mechanical_loads_sub_model_part_list": [],
            "loads_sub_model_part_list": [],
            "loads_variable_list": []
        }
    }))


def _get_computing_domain(model):
    root = model.GetModelPart("MainModelPart")
    return root, root.GetSubModelPart("mechanical_computing_domain")


def _extract_results(model):
    """Extracts the observable nodal quantities after the analysis: per-node
    normalized NODAL_CAUCHY_STRESS_TENSOR / NODAL_AREA and the expected
    (single-accumulation) NODAL_AREA."""
    root, computing_domain = _get_computing_domain(model)

    expected_area = {}
    for element in computing_domain.Elements:
        volume = element.GetGeometry().Volume()
        for node in element.GetNodes():
            expected_area[node.Id] = expected_area.get(node.Id, 0.0) + volume

    node_data = {}
    for node in computing_domain.Nodes:
        stress = node.GetSolutionStepValue(KratosPoro.NODAL_CAUCHY_STRESS_TENSOR)
        node_data[node.Id] = {
            "stress": [stress[i, j] for i in range(stress.Size1()) for j in range(stress.Size2())],
            "area": node.GetSolutionStepValue(KM.NODAL_AREA),
            "expected_area": expected_area[node.Id],
        }
    return node_data


class DamProcessBasedNodalSmoothingTest(KratosUnittest.TestCase):

    def test_old_project_parameters_without_setting(self):
        """An old ProjectParameters file must run normally and produce the
        smoothed nodal Cauchy stress through the single scheme/process
        implementation, with a single accumulation."""
        with KratosUnittest.WorkFolderScope(".", __file__):
            model = _run_construction()
            node_data = _extract_results(model)

            # Single accumulation of the element measures (no element-level
            # accumulation and no double accumulation).
            max_area_error = 0.0
            for node in node_data.values():
                max_area_error = max(max_area_error, abs(node["area"] - node["expected_area"]))
            self.assertLessEqual(max_area_error, 1.0e-9)

            # The final normalized nodal Cauchy stress is available.
            max_stress = 0.0
            for node in node_data.values():
                max_stress = max(max_stress, max(abs(s) for s in node["stress"]))
            self.assertGreater(max_stress, 1.0e-6)

    def test_mechanical_solver_construction(self):
        """The Dam mechanical solver constructs from an old-style
        ProjectParameters mechanical settings block."""
        from KratosMultiphysics.DamApplication import dam_mechanical_solver

        model = KM.Model()
        model_part = model.CreateModelPart("MainModelPart")
        solver = dam_mechanical_solver.DamMechanicalSolver(
            model_part, _mechanical_solver_settings())
        solver._ConstructScheme("Newmark", "Quasi-Static")

    def test_smoothing_scheme_construction_all_schemes(self):
        """All three Dam smoothing schemes construct through the solver without
        without any Dam-specific options (single process implementation)."""
        from KratosMultiphysics.DamApplication import dam_mechanical_solver

        def make_params(rayleigh):
            return KM.Parameters(json.dumps({
                "solver_type": "dam_mechanical_solver",
                "model_import_settings": {"input_type": "mdpa", "input_filename": "x", "input_file_label": 0},
                "echo_level": 0,
                "buffer_size": 2,
                "processes_sub_model_part_list": [""],
                "mechanical_solver_settings": {
                    "echo_level": 0,
                    "reform_dofs_at_each_step": False,
                    "clear_storage": False,
                    "compute_reactions": False,
                    "move_mesh_flag": True,
                    "solution_type": "Quasi-Static",
                    "scheme_type": "Newmark",
                    "rayleigh_m": rayleigh,
                    "rayleigh_k": rayleigh,
                    "strategy_type": "Newton-Raphson",
                    "convergence_criterion": "Displacement_criterion",
                    "displacement_relative_tolerance": 1.0e-4,
                    "displacement_absolute_tolerance": 1.0e-9,
                    "residual_relative_tolerance": 1.0e-4,
                    "residual_absolute_tolerance": 1.0e-9,
                    "max_iteration": 15,
                    "desired_iterations": 4,
                    "max_radius_factor": 20.0,
                    "min_radius_factor": 0.5,
                    "block_builder": True,
                    "nonlocal_damage": False,
                    "characteristic_length": 0.05,
                    "search_neighbours_step": False,
                    "linear_solver_settings": {
                        "solver_type": "bicgstab",
                        "max_iteration": 200,
                        "tolerance": 1.0e-7,
                        "preconditioner_type": "ilu0"
                    },
                    "problem_domain_sub_model_part_list": [""],
                    "body_domain_sub_model_part_list": [],
                    "mechanical_loads_sub_model_part_list": [],
                    "loads_sub_model_part_list": [],
                    "loads_variable_list": []
                }
            }))

        for rayleigh, expected_scheme in [
                (0.0, "IncrementalUpdateStaticSmoothingScheme"),
                (1.0e-3, "IncrementalUpdateStaticDampedSmoothingScheme")]:
            model = KM.Model()
            model_part = model.CreateModelPart("S" + str(rayleigh))
            solver = dam_mechanical_solver.DamMechanicalSolver(model_part, make_params(rayleigh))
            scheme = solver._ConstructScheme("Newmark", "Quasi-Static")
            self.assertEqual(type(scheme).__name__, expected_scheme)

        model = KM.Model()
        model_part = model.CreateModelPart("Bossak")
        solver = dam_mechanical_solver.DamMechanicalSolver(model_part, make_params(0.0))
        scheme = solver._ConstructScheme("Newmark", "Dynamic")
        self.assertEqual(type(scheme).__name__, "BossakDisplacementSmoothingScheme")


if __name__ == '__main__':
    KratosUnittest.runTests(KratosUnittest.KratosSuites)
