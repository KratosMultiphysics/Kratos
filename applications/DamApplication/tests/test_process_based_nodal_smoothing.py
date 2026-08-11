import json
import os

import KratosMultiphysics as KM
import KratosMultiphysics.DamApplication as KratosDam
import KratosMultiphysics.PoromechanicsApplication as KratosPoro
import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.DamApplication.dam_analysis import DamAnalysis


def _load_construction_parameters():
    """Loads the old-style construction ProjectParameters file.

    This file does NOT contain 'use_process_based_nodal_stress_extrapolation',
    representing an existing/older user ProjectParameters file.
    """
    parameters_file = os.path.join(
        os.path.dirname(__file__), "construction", "construction_parameters.json")
    with open(parameters_file) as f:
        return json.load(f)


def _run_variant(mode):
    """Runs the construction analysis.

    mode = None  : the setting is NOT present (old ProjectParameters file)
    mode = True  : explicit use_process_based_nodal_stress_extrapolation = true
    mode = False : explicit use_process_based_nodal_stress_extrapolation = false
    Returns the resulting Model.
    """
    parameters = _load_construction_parameters()
    mechanical_settings = parameters["solver_settings"]["mechanical_solver_settings"]
    if mode is None:
        mechanical_settings.pop("use_process_based_nodal_stress_extrapolation", None)
    else:
        mechanical_settings["use_process_based_nodal_stress_extrapolation"] = bool(mode)

    model = KM.Model()
    simulation = DamAnalysis(model, KM.Parameters(json.dumps(parameters)))
    simulation.Run()
    return model


def _get_computing_domain(model):
    root = model.GetModelPart("MainModelPart")
    return root, root.GetSubModelPart("mechanical_computing_domain")


def _extract_results(model):
    """Extracts the observable nodal and integration-point quantities after the
    analysis: ProcessInfo resolution, per-node DISPLACEMENT / normalized
    NODAL_CAUCHY_STRESS_TENSOR / NODAL_AREA, the expected (single-accumulation)
    NODAL_AREA, and the per-element integration-point CAUCHY_STRESS_TENSOR.
    """
    root, computing_domain = _get_computing_domain(model)
    resolved_flag = root.ProcessInfo[KratosDam.USE_PROCESS_BASED_NODAL_CAUCHY_STRESS_EXTRAPOLATION]

    expected_area = {}
    for element in computing_domain.Elements:
        volume = element.GetGeometry().Volume()
        for node in element.GetNodes():
            expected_area[node.Id] = expected_area.get(node.Id, 0.0) + volume

    node_data = {}
    for node in computing_domain.Nodes:
        displacement = node.GetSolutionStepValue(KM.DISPLACEMENT)
        stress = node.GetSolutionStepValue(KratosPoro.NODAL_CAUCHY_STRESS_TENSOR)
        node_data[node.Id] = {
            "displacement": [float(displacement[i]) for i in range(len(displacement))],
            "stress": [stress[i, j] for i in range(stress.Size1()) for j in range(stress.Size2())],
            "area": node.GetSolutionStepValue(KM.NODAL_AREA),
            "expected_area": expected_area[node.Id],
        }

    ip_stress = {}
    for element in computing_domain.Elements:
        out = element.CalculateOnIntegrationPoints(KM.CAUCHY_STRESS_TENSOR, computing_domain.ProcessInfo)
        ip_stress[element.Id] = out

    return {
        "resolved_flag": resolved_flag,
        "node_data": node_data,
        "ip_stress": ip_stress,
    }


def _tolerance(reference_value):
    return max(1.0e-12, 1.0e-10 * abs(reference_value))


def _assert_close(a, b, label):
    for key in ("displacement", "stress", "area"):
        av = a[key]
        bv = b[key]
        if key == "area":
            assert abs(av - bv) <= _tolerance(bv), f"{label}: {key} {av} vs {bv}"
        else:
            assert len(av) == len(bv), f"{label}: {key} size"
            for i in range(len(av)):
                assert abs(av[i] - bv[i]) <= _tolerance(bv[i]), f"{label}: {key}[{i}] {av[i]} vs {bv[i]}"


class DamProcessBasedNodalSmoothingTest(KratosUnittest.TestCase):

    def test_old_project_parameters_without_setting(self):
        """Old ProjectParameters without the setting must resolve to the
        process-based default (true), run without error, and produce a single
        nodal accumulation even with legacy thermo-mechanical elements."""
        with KratosUnittest.WorkFolderScope(".", __file__):
            model = _run_variant(None)
            results = _extract_results(model)

            # The absent setting resolves to the process-based default.
            self.assertTrue(results["resolved_flag"])

            # Legacy Dam thermo-mechanical elements are still present in the
            # model. Their Info() string is "Large Displacement Element",
            # which is distinct from the StructuralMechanics element Info()
            # string ("Small Displacement Solid Element").
            root, computing_domain = _get_computing_domain(model)
            element_names = {e.Info() for e in computing_domain.Elements}
            self.assertTrue(any("Large Displacement Element" in name for name in element_names))

            # ...and the nodal smoothing shows a SINGLE accumulation of the
            # element measures (no legacy+process double accumulation).
            max_area_error = 0.0
            for node_data in results["node_data"].values():
                max_area_error = max(max_area_error, abs(node_data["area"] - node_data["expected_area"]))
            self.assertLessEqual(max_area_error, 1.0e-9)

            # The final normalized nodal Cauchy stress is available for the
            # downstream consumers (selfweight/initial-stress transfer).
            max_stress = 0.0
            for node_data in results["node_data"].values():
                max_stress = max(max_stress, max(abs(s) for s in node_data["stress"]))
            self.assertGreater(max_stress, 1.0e-6)

    def test_default_true_false_three_way_equivalence(self):
        """A (absent -> true), B (explicit true) and C (explicit false) must be
        numerically equivalent; A and B follow the same process-based path."""
        with KratosUnittest.WorkFolderScope(".", __file__):
            model_a = _run_variant(None)
            model_b = _run_variant(True)
            model_c = _run_variant(False)
            results_a = _extract_results(model_a)
            results_b = _extract_results(model_b)
            results_c = _extract_results(model_c)

            # A and B resolve to process-based; C to legacy.
            self.assertTrue(results_a["resolved_flag"])
            self.assertTrue(results_b["resolved_flag"])
            self.assertFalse(results_c["resolved_flag"])

            # A and B are the same process-based path: identical results.
            for nid in results_a["node_data"]:
                _assert_close(results_a["node_data"][nid], results_b["node_data"][nid], f"A/B node {nid}")
            self.assertEqual(sorted(results_a["ip_stress"]), sorted(results_b["ip_stress"]))
            for eid in results_a["ip_stress"]:
                for gp, (sa, sb) in enumerate(zip(results_a["ip_stress"][eid], results_b["ip_stress"][eid])):
                    for i in range(sa.Size1()):
                        for j in range(sa.Size2()):
                            self.assertLessEqual(abs(sa[i, j] - sb[i, j]), _tolerance(sb[i, j]))

            # A, B and C numerically equivalent within the established
            # tolerances (displacement, nodal stress, NODAL_AREA, IP stress).
            for nid in results_a["node_data"]:
                _assert_close(results_a["node_data"][nid], results_c["node_data"][nid], f"A/C node {nid}")
                _assert_close(results_b["node_data"][nid], results_c["node_data"][nid], f"B/C node {nid}")
            for eid in results_a["ip_stress"]:
                for gp, (sa, sc) in enumerate(zip(results_a["ip_stress"][eid], results_c["ip_stress"][eid])):
                    for i in range(sa.Size1()):
                        for j in range(sa.Size2()):
                            self.assertLessEqual(abs(sa[i, j] - sc[i, j]), _tolerance(sc[i, j]))

    def test_default_resolution_all_smoothing_schemes(self):
        """An omitted setting selects process-based ownership for every Dam
        smoothing scheme (static, damped static and Bossak)."""
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

        # Static (Quasi-Static, no Rayleigh damping).
        model = KM.Model()
        model_part = model.CreateModelPart("Static")
        solver = dam_mechanical_solver.DamMechanicalSolver(model_part, make_params(0.0))
        scheme = solver._ConstructScheme("Newmark", "Quasi-Static")
        self.assertTrue(model_part.ProcessInfo[KratosDam.USE_PROCESS_BASED_NODAL_CAUCHY_STRESS_EXTRAPOLATION])
        self.assertEqual(type(scheme).__name__, "IncrementalUpdateStaticSmoothingScheme")

        # Damped static (Quasi-Static with Rayleigh damping).
        model = KM.Model()
        model_part = model.CreateModelPart("Damped")
        solver = dam_mechanical_solver.DamMechanicalSolver(model_part, make_params(1.0e-3))
        scheme = solver._ConstructScheme("Newmark", "Quasi-Static")
        self.assertTrue(model_part.ProcessInfo[KratosDam.USE_PROCESS_BASED_NODAL_CAUCHY_STRESS_EXTRAPOLATION])
        self.assertEqual(type(scheme).__name__, "IncrementalUpdateStaticDampedSmoothingScheme")

        # Bossak (dynamic solution type).
        model = KM.Model()
        model_part = model.CreateModelPart("Bossak")
        solver = dam_mechanical_solver.DamMechanicalSolver(model_part, make_params(0.0))
        scheme = solver._ConstructScheme("Newmark", "Dynamic")
        self.assertTrue(model_part.ProcessInfo[KratosDam.USE_PROCESS_BASED_NODAL_CAUCHY_STRESS_EXTRAPOLATION])
        self.assertEqual(type(scheme).__name__, "BossakDisplacementSmoothingScheme")

    def test_structural_mechanics_only_default_smoothing(self):
        """A model containing only StructuralMechanicsApplication
        SmallDisplacement elements (no Dam thermo-mechanical element) runs the
        real Dam smoothing scheme with the default settings (no explicit
        setting) and produces the nodal smoothing."""
        from KratosMultiphysics.DamApplication import dam_mechanical_solver

        model = KM.Model()
        model_part = model.CreateModelPart("SMAOnly", 2)
        model_part.ProcessInfo[KM.DOMAIN_SIZE] = 3
        model_part.ProcessInfo[KM.SPACE_DIMENSION] = 3
        model_part.ProcessInfo[KM.IS_RESTARTED] = False

        for variable in [KM.DISPLACEMENT, KM.VELOCITY, KM.ACCELERATION,
                         KM.VOLUME_ACCELERATION, KM.TEMPERATURE,
                         KratosDam.NODAL_REFERENCE_TEMPERATURE,
                         KratosPoro.NODAL_CAUCHY_STRESS_TENSOR, KM.NODAL_AREA,
                         KratosPoro.NODAL_JOINT_AREA, KratosPoro.NODAL_JOINT_WIDTH,
                         KratosDam.INITIAL_NODAL_CAUCHY_STRESS_TENSOR]:
            model_part.AddNodalSolutionStepVariable(variable)

        coordinates = [[0, 0, 0], [2, 0, 0], [2, 1, 0], [0, 1, 0],
                       [0, 0, 1], [2, 0, 1], [2, 1, 1], [0, 1, 1],
                       [3, 0, 0], [3, 1, 0], [3, 0, 1], [3, 1, 1]]
        for i, coordinate in enumerate(coordinates):
            node = model_part.CreateNewNode(i + 1, *coordinate)
            node.AddDof(KM.DISPLACEMENT_X)
            node.AddDof(KM.DISPLACEMENT_Y)
            node.AddDof(KM.DISPLACEMENT_Z)
            node.SetSolutionStepValue(KratosDam.NODAL_REFERENCE_TEMPERATURE, 0, 20.0)

        properties = model_part.CreateNewProperties(1)
        properties[KM.YOUNG_MODULUS] = 2.0e7
        properties[KM.POISSON_RATIO] = 0.2
        properties[KratosDam.THERMAL_EXPANSION] = 1.0e-5
        properties[KM.DENSITY] = 2400.0
        properties.SetValue(KM.CONSTITUTIVE_LAW, KratosDam.ThermalLinearElastic3DLaw())

        model_part.CreateNewElement("SmallDisplacementElement3D8N", 1, [1, 2, 3, 4, 5, 6, 7, 8], properties)
        model_part.CreateNewElement("SmallDisplacementElement3D8N", 2, [2, 9, 10, 3, 6, 11, 12, 7], properties)

        for i, node in enumerate(model_part.Nodes):
            x0 = node.GetInitialPosition()
            displacement = node.GetSolutionStepValue(KM.DISPLACEMENT)
            displacement[0] = 1.0e-3 * x0[0] + 3.0e-4 * x0[1] + 1.0e-4 * x0[0] * x0[1] + 2.0e-4 * x0[0] * x0[2]
            displacement[1] = 3.0e-4 * x0[0] - 1.0e-3 * x0[1] + 2.0e-4 * x0[0] ** 2 + 1.0e-4 * x0[1] * x0[2]
            displacement[2] = 5.0e-4 * x0[2] + 2.0e-4 * x0[0] + 1.0e-4 * x0[1] + 1.0e-4 * x0[0] * x0[1]
            node.SetSolutionStepValue(KM.DISPLACEMENT, 0, displacement)
            node.X = x0[0] + displacement[0]
            node.Y = x0[1] + displacement[1]
            node.Z = x0[2] + displacement[2]
            node.SetSolutionStepValue(KM.TEMPERATURE, 0, 20.0 + 5.0 * i)

        for element in model_part.Elements:
            element.Initialize(model_part.ProcessInfo)

        # Default resolution (no explicit setting) via the real Dam solver.
        solver = dam_mechanical_solver.DamMechanicalSolver(model_part, self._default_solver_settings())
        scheme = solver._ConstructScheme("Newmark", "Quasi-Static")
        self.assertTrue(model_part.ProcessInfo[KratosDam.USE_PROCESS_BASED_NODAL_CAUCHY_STRESS_EXTRAPOLATION])

        # Run the real Dam smoothing scheme finalization.
        system_matrix = KM.CompressedMatrix()
        system_vector = KM.Vector()
        scheme.FinalizeSolutionStep(model_part, system_matrix, system_vector, system_vector)

        # Single accumulation of the element measures.
        expected_area = {}
        for element in model_part.Elements:
            volume = element.GetGeometry().Volume()
            for node in element.GetNodes():
                expected_area[node.Id] = expected_area.get(node.Id, 0.0) + volume
        for node in model_part.Nodes:
            self.assertAlmostEqual(
                node.GetSolutionStepValue(KM.NODAL_AREA),
                expected_area[node.Id],
                delta=1.0e-9)

        # Non-trivial smoothed nodal Cauchy stress.
        for node in model_part.Nodes:
            stress = node.GetSolutionStepValue(KratosPoro.NODAL_CAUCHY_STRESS_TENSOR)
            self.assertGreater(max(abs(stress[i, j]) for i in range(3) for j in range(3)), 1.0e-6)

    @staticmethod
    def _default_solver_settings():
        import json as _json
        return KM.Parameters(_json.dumps({
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


if __name__ == '__main__':
    KratosUnittest.runTests(KratosUnittest.KratosSuites)
