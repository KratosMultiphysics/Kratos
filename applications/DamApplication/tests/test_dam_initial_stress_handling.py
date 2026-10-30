import os
import sys
import json
import subprocess
import tempfile

import KratosMultiphysics as KM
import KratosMultiphysics.DamApplication as KratosDam
import KratosMultiphysics.PoromechanicsApplication as KratosPoro
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.testing.utilities import GetPython3Command

from KratosMultiphysics.DamApplication.dam_mechanical_solver import DamMechanicalSolver
from KratosMultiphysics.DamApplication.dam_UP_mechanical_solver import DamUPMechanicalSolver

TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
if TESTS_DIR not in sys.path:
    sys.path.insert(0, TESTS_DIR)

import test_interface_tangent_consistency as interface_test

_MDPA = """Begin ModelPartData
End ModelPartData

Begin Properties 0
End Properties

Begin Properties 1
    DENSITY 2400.0
    YOUNG_MODULUS 26000000000.0
    POISSON_RATIO 0.2
    THERMAL_EXPANSION 1.0e-5
    CONSTITUTIVE_LAW_NAME ThermalLinearElastic2DPlaneStrain
End Properties

Begin Nodes
1 0.0 0.0 0.0
2 1.0 0.0 0.0
3 1.0 1.0 0.0
4 0.0 1.0 0.0
End Nodes

Begin Elements SmallDisplacementThermoMechanicElement2D4N
1 1 1 2 3 4
End Elements

Begin SubModelPart Parts_Concrete
    Begin SubModelPartNodes
    1
    2
    3
    4
    End SubModelPartNodes
    Begin SubModelPartElements
    1
    End SubModelPartElements
End SubModelPart

Begin SubModelPart DISPLACEMENT_Displacement_Condition
    Begin SubModelPartNodes
    1
    2
    3
    4
    End SubModelPartNodes
End SubModelPart
"""

_PREVIOUS_MDPA = _MDPA + """
Begin NodalData DISPLACEMENT
1 0 [3](0.0,0.0,0.0)
2 0 [3](0.01,0.0,0.0)
3 0 [3](0.01,0.0,0.0)
4 0 [3](0.0,0.0,0.0)
End NodalData

Begin NodalData NODAL_CAUCHY_STRESS_TENSOR
1 0 [2,2]((2.888888888888889e8,0.0),(0.0,7.222222222222222e7))
2 0 [2,2]((2.888888888888889e8,0.0),(0.0,7.222222222222222e7))
3 0 [2,2]((2.888888888888889e8,0.0),(0.0,7.222222222222222e7))
4 0 [2,2]((2.888888888888889e8,0.0),(0.0,7.222222222222222e7))
End NodalData
"""

_CHILD = r"""import json
import os

import KratosMultiphysics as KM
from KratosMultiphysics.DamApplication.dam_analysis import DamAnalysis

with open(os.path.join(os.getcwd(), "parameters.json"), "r") as parameter_file:
    parameters = KM.Parameters(parameter_file.read())

DamAnalysis(KM.Model(), parameters).Run()
"""

_EXPECTED_SIGMA_XX = 2.888888888888889e8
_EXPECTED_SIGMA_YY = 7.222222222222222e7


class TestDamInitialStressHandling(KratosUnittest.TestCase):
    """Regression tests for INITIAL_STRESS_TENSOR storage in the Dam schemes."""

    @staticmethod
    def _make_model_part():
        model = KM.Model()
        model_part = model.CreateModelPart("initial_stress_test")
        model_part.AddNodalSolutionStepVariable(KratosPoro.INITIAL_STRESS_TENSOR)
        model_part.ProcessInfo.SetValue(KM.DELTA_TIME, 0.1)
        model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
        return model_part

    @staticmethod
    def _schemes():
        return [
            KratosDam.IncrementalUpdateStaticSmoothingScheme(),
            KratosDam.IncrementalUpdateStaticDampedSmoothingScheme(0.0, 0.0),
            KratosDam.BossakDisplacementSmoothingScheme(0.0, 0.0, 0.0),
            KratosDam.DamUPScheme(0.25, 0.5, 0.0, 0.0),
            KratosDam.DamPScheme(0.25, 0.5),
        ]

    def _assert_tensor(self, tensor, expected):
        self.assertEqual(tensor.Size1(), 3)
        self.assertEqual(tensor.Size2(), 3)
        for i in range(3):
            for j in range(3):
                self.assertAlmostEqual(tensor[i, j], expected[i][j], places=14)

    def test_schemes_preserve_preloaded_3x3_initial_stress(self):
        expected = [[1.0e6, 2.0e6, 0.0], [2.0e6, 3.0e6, 0.0], [0.0, 0.0, 4.0e6]]
        for scheme in self._schemes():
            with self.subTest(scheme=type(scheme).__name__):
                model_part = self._make_model_part()
                for node in model_part.Nodes:
                    node.SetSolutionStepValue(KratosPoro.INITIAL_STRESS_TENSOR,
                                              KM.Matrix(expected))
                scheme.Initialize(model_part)
                for node in model_part.Nodes:
                    self._assert_tensor(
                        node.GetSolutionStepValue(KratosPoro.INITIAL_STRESS_TENSOR),
                        expected)

    def test_schemes_initialize_empty_tensor_to_zero_3x3(self):
        for scheme in self._schemes():
            with self.subTest(scheme=type(scheme).__name__):
                model_part = self._make_model_part()
                scheme.Initialize(model_part)
                for node in model_part.Nodes:
                    self._assert_tensor(
                        node.GetSolutionStepValue(KratosPoro.INITIAL_STRESS_TENSOR),
                        [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]])

    def test_schemes_convert_2x2_initial_stress_to_3x3(self):
        expected = [[1.0e6, 2.0e6, 0.0], [2.0e6, 3.0e6, 0.0], [0.0, 0.0, 0.0]]
        for scheme in self._schemes():
            with self.subTest(scheme=type(scheme).__name__):
                model_part = self._make_model_part()
                initial_stress = KM.Matrix(2, 2)
                initial_stress[0, 0] = 1.0e6
                initial_stress[0, 1] = 2.0e6
                initial_stress[1, 0] = 2.0e6
                initial_stress[1, 1] = 3.0e6
                for node in model_part.Nodes:
                    node.SetSolutionStepValue(KratosPoro.INITIAL_STRESS_TENSOR,
                                              initial_stress)
                scheme.Initialize(model_part)
                scheme.Initialize(model_part)  # idempotence
                for node in model_part.Nodes:
                    self._assert_tensor(
                        node.GetSolutionStepValue(KratosPoro.INITIAL_STRESS_TENSOR),
                        expected)

    def test_schemes_reject_unsupported_tensor_sizes(self):
        for unsupported_size in ((1, 1), (4, 4)):
            for scheme in self._schemes():
                with self.subTest(scheme=type(scheme).__name__,
                                  size=unsupported_size):
                    model_part = self._make_model_part()
                    matrix = KM.Matrix(unsupported_size[0], unsupported_size[1])
                    for node in model_part.Nodes:
                        node.SetSolutionStepValue(KratosPoro.INITIAL_STRESS_TENSOR,
                                                  matrix)
                    with self.assertRaises(RuntimeError):
                        scheme.Initialize(model_part)

    def test_interface_initial_stress_survives_scheme_initialization(self):
        test_case = interface_test.TestSmallDisplacementInterfaceTangent(
            "test_initial_gap_uses_reference_geometry")

        initial_stress = KM.Matrix(3, 3)
        for i in range(3):
            for j in range(3):
                initial_stress[i, j] = 1.0e6 if i < 2 and j < 2 else 0.0

        model_part, element = test_case._make_element()
        for node in model_part.Nodes:
            node.SetSolutionStepValue(KratosPoro.INITIAL_STRESS_TENSOR,
                                      initial_stress)
        _, rhs_before = test_case._local_system(element, model_part.ProcessInfo)

        KratosDam.IncrementalUpdateStaticSmoothingScheme().Initialize(model_part)

        for node in model_part.Nodes:
            self._assert_tensor(
                node.GetSolutionStepValue(KratosPoro.INITIAL_STRESS_TENSOR),
                [[1.0e6, 1.0e6, 0.0], [1.0e6, 1.0e6, 0.0], [0.0, 0.0, 0.0]])

        _, rhs_after = test_case._local_system(element, model_part.ProcessInfo)
        test_case._assert_vector_near(rhs_before, rhs_after)

        zero_model_part, zero_element = test_case._make_element()
        _, rhs_zero = test_case._local_system(
            zero_element, zero_model_part.ProcessInfo)
        self.assertGreater(
            max(abs(a - b) for a, b in zip(rhs_before, rhs_zero)), 1.0)

    @staticmethod
    def _base_parameters(problem_name, add_displacement, end_time):
        displacement_value = ["0.01*x", "0.0", "0.0"] if add_displacement else ["0.0", "0.0", "0.0"]
        return {
            "problem_data": {
                "problem_name": problem_name,
                "model_part_name": "MainModelPart",
                "domain_size": 2,
                "parallel_type": "OpenMP",
                "number_of_threads": 1,
                "start_time": 0.0,
                "end_time": end_time,
                "time_step": 1.0,
                "time_scale": "Seconds",
                "consider_construction": False,
                "streamlines_utility": False,
            },
            "solver_settings": {
                "solver_type": "dam_mechanical_solver",
                "model_import_settings": {
                    "input_type": "mdpa",
                    "input_filename": "mini",
                    "input_file_label": 0,
                },
                "echo_level": 0,
                "buffer_size": 2,
                "processes_sub_model_part_list": ["DISPLACEMENT_Displacement_Condition"],
                "mechanical_solver_settings": {
                    "echo_level": 0,
                    "reform_dofs_at_each_step": False,
                    "clear_storage": False,
                    "compute_reactions": True,
                    "move_mesh_flag": False,
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
                        "preconditioner_type": "ilu0",
                        "scaling": False,
                    },
                    "problem_domain_sub_model_part_list": ["Parts_Concrete"],
                    "body_domain_sub_model_part_list": ["Parts_Concrete"],
                    "mechanical_loads_sub_model_part_list": [],
                    "loads_sub_model_part_list": [],
                    "loads_variable_list": [],
                },
            },
            "output_configuration": {
                "result_file_configuration": {
                    "gidpost_flags": {
                        "GiDPostMode": "GiD_PostAscii",
                        "WriteDeformedMeshFlag": "WriteDeformed",
                        "WriteConditionsFlag": "WriteElementsOnly",
                        "MultiFileFlag": "SingleFile",
                    },
                    "output_control_type": "time_s",
                    "output_interval": 1.0,
                    "start_output_results": 0,
                    "body_output": True,
                    "node_output": False,
                    "skin_output": False,
                    "plane_output": [],
                    "nodal_results": ["DISPLACEMENT", "NODAL_CAUCHY_STRESS_TENSOR"],
                    "gauss_point_results": ["CAUCHY_STRESS_TENSOR"],
                },
                "point_data_configuration": [],
            },
            "constraints_process_list": [{
                "python_module": "assign_vector_variable_process",
                "kratos_module": "KratosMultiphysics",
                "process_name": "AssignVectorVariableProcess",
                "Parameters": {
                    "model_part_name": "MainModelPart.DISPLACEMENT_Displacement_Condition",
                    "variable_name": "DISPLACEMENT",
                    "constrained": True,
                    "interval": [0.0, 1.0e30],
                    "value": displacement_value,
                },
            }],
            "loads_process_list": [],
            "construction_process": [],
            "transfer_results_process": {
                "save_intermediate_mechanical_variables": False,
                "save_intermediate_thermal_variables": False,
                "save_intermediate_variables_step": 0,
                "save_final_mechanical_variables": False,
                "save_final_thermal_variables": False,
                "add_previous_results": True,
                "type_of_results": "Mechanical",
                "add_displacement": add_displacement,
                "add_stress": True,
                "add_temperature": True,
                "add_reference_temperature": True,
                "file_name_mechanical": "mini_mechanical.mdpa",
                "file_name_thermal": "mini_thermal.mdpa",
            },
            "temperature_by_device_list": [],
            "output_device_list": [],
        }

    def _run_transfer_case(self, problem_name, add_displacement, end_time):
        with tempfile.TemporaryDirectory() as working_directory:
            with open(os.path.join(working_directory, "mini.mdpa"), "w") as mdpa_file:
                mdpa_file.write(_MDPA)
            with open(os.path.join(working_directory, "mini_mechanical.mdpa"), "w") as mdpa_file:
                mdpa_file.write(_PREVIOUS_MDPA)
            with open(os.path.join(working_directory, "parameters.json"), "w") as parameter_file:
                json.dump(self._base_parameters(problem_name, add_displacement, end_time),
                          parameter_file, indent=4)
            with open(os.path.join(working_directory, "run_case.py"), "w") as run_file:
                run_file.write(_CHILD)

            process = subprocess.Popen(
                [GetPython3Command(), "run_case.py"],
                cwd=working_directory,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT)
            output = process.communicate(timeout=300)[0].decode("ascii", "replace")
            self.assertEqual(
                process.returncode, 0,
                msg="Dam initial-stress transfer case failed (exit {}):\n{}".format(
                    process.returncode, output))

            result_file = os.path.join(working_directory, problem_name + ".post.res")
            self.assertTrue(os.path.isfile(result_file), "GiD result file was not written")
            return {
                step: {
                    "nodal": self._read_nodal_result(result_file, step),
                    "gp": self._read_first_gauss_point(result_file, step),
                }
                for step in range(end_time + 1)
            }

    @staticmethod
    def _read_nodal_result(result_file, label):
        header = 'Result "NODAL_CAUCHY_STRESS_TENSOR" "Kratos" {} Matrix OnNodes'.format(label)
        with open(result_file, "r") as file:
            lines = [line.strip() for line in file]
        header_index = lines.index(header)
        values = {}
        for line in lines[header_index + 2:]:
            if line == "End Values":
                break
            parts = line.split()
            if parts:
                values[int(parts[0])] = [float(value) for value in parts[1:]]
        return values

    @staticmethod
    def _read_first_gauss_point(result_file, label):
        header = 'Result "CAUCHY_STRESS_TENSOR" "Kratos" {} Matrix OnGaussPoints'.format(label)
        with open(result_file, "r") as file:
            lines = file.read().splitlines()
        for index, line in enumerate(lines):
            if line.startswith(header):
                break
        for line in lines[index + 2:]:
            if line.strip() == "End Values":
                break
            if line.strip():
                parts = line.split()
                if line[0] != " ":
                    return [float(value) for value in parts[1:]]
                return [float(value) for value in parts]
        return []

    def _assert_stress_xx_yy(self, values, expected_xx, expected_yy):
        self.assertAlmostEqual(values[0], expected_xx, delta=abs(expected_xx) * 1.0e-5)
        self.assertAlmostEqual(values[1], expected_yy, delta=abs(expected_yy) * 1.0e-5)

    def test_nodal_stress_not_double_counted_with_mapped_displacement(self):
        results = self._run_transfer_case("dc_displacement", True, 2)

        for step in (0, 1, 2):
            self._assert_stress_xx_yy(
                results[step]["nodal"][1], _EXPECTED_SIGMA_XX, _EXPECTED_SIGMA_YY)

        for step in (1, 2):
            self._assert_stress_xx_yy(
                results[step]["gp"], _EXPECTED_SIGMA_XX, _EXPECTED_SIGMA_YY)
            self._assert_stress_xx_yy(
                results[step]["nodal"][1],
                results[step]["gp"][0],
                results[step]["gp"][1])

    def test_stress_only_transfer_preserves_initial_nodal_output(self):
        results = self._run_transfer_case("dc_stress_only", False, 2)

        for step in (0, 1, 2):
            self._assert_stress_xx_yy(
                results[step]["nodal"][1], _EXPECTED_SIGMA_XX, _EXPECTED_SIGMA_YY)

        for step in (1, 2):
            self.assertLess(abs(results[step]["gp"][0]), 1.0e-6)
            self.assertLess(abs(results[step]["gp"][1]), 1.0e-6)

    def test_nodal_stress_recomputation_flag_is_scheme_aware(self):
        settings = self._base_parameters("flag", True, 1)["solver_settings"]

        mechanical_settings = KM.Parameters(json.dumps(settings))
        mechanical_solver = DamMechanicalSolver(
            KM.Model().CreateModelPart("mechanical"), mechanical_settings)
        mechanical_solver._ConstructScheme("Newmark", "Quasi-Static")
        self.assertTrue(mechanical_solver._recomputes_nodal_cauchy_stress)

        up_settings_dict = self._base_parameters("flag", True, 1)["solver_settings"]
        up_settings_dict["solver_type"] = "dam_UP_mechanical_solver"
        up_solver = DamUPMechanicalSolver(
            KM.Model().CreateModelPart("up"), KM.Parameters(json.dumps(up_settings_dict)))
        up_solver._ConstructScheme("Newmark", "Quasi-Static")
        self.assertFalse(getattr(up_solver, "_recomputes_nodal_cauchy_stress", False))


if __name__ == "__main__":
    KratosUnittest.main()
