import json

import KratosMultiphysics as KM
import KratosMultiphysics.DamApplication as KratosDam
import KratosMultiphysics.KratosUnittest as KratosUnittest


def _make_solver_settings(nonlocal_damage):
    """Builds an OLD-style mechanical solver settings dict: it contains only the
    existing 'nonlocal_damage' setting and NO internal ownership option."""
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
            "nonlocal_damage": nonlocal_damage,
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


class DamNonlocalOwnershipTest(KratosUnittest.TestCase):

    def test_old_project_parameters_nonlocal_damage(self):
        """An old ProjectParameters with 'nonlocal_damage': true (and no new
        ownership setting) must resolve the internal process-based LOCAL
        ownership flag to true through the existing solver."""
        from KratosMultiphysics.DamApplication import dam_mechanical_solver

        model = KM.Model()
        model_part = model.CreateModelPart("OldNonlocal")
        solver = dam_mechanical_solver.DamMechanicalSolver(model_part, _make_solver_settings(True))
        solver._ConstructScheme("Newmark", "Quasi-Static")
        self.assertTrue(
            model_part.ProcessInfo[KratosDam.USE_PROCESS_BASED_LOCAL_EQUIVALENT_STRAIN])

    def test_nonlocal_damage_false_no_local_ownership(self):
        """'nonlocal_damage': false must NOT activate the process-based LOCAL
        ownership."""
        from KratosMultiphysics.DamApplication import dam_mechanical_solver

        model = KM.Model()
        model_part = model.CreateModelPart("NoNonlocal")
        solver = dam_mechanical_solver.DamMechanicalSolver(model_part, _make_solver_settings(False))
        solver._ConstructScheme("Newmark", "Quasi-Static")
        self.assertFalse(
            model_part.ProcessInfo[KratosDam.USE_PROCESS_BASED_LOCAL_EQUIVALENT_STRAIN])


if __name__ == '__main__':
    KratosUnittest.runTests(KratosUnittest.KratosSuites)
