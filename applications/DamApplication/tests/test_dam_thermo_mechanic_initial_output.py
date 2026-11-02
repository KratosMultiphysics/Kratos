import json
import os
import subprocess
import tempfile

import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.testing.utilities import GetPython3Command

# Minimal two-dimensional thermo-mechanical case. The historical Dam element
# names resolve to the StructuralMechanicsApplication SmallDisplacement runtime
# element, and the thermal mesh uses the ConvectionDiffusion EulerianConvDiff2D.
# A single triangle keeps the case fast while still exercising the real
# constitutive-law initialization and the real GiD Gauss-point writer.
_MDPA = """Begin ModelPartData
End ModelPartData

Begin Properties 0
End Properties

Begin Properties 1
    DENSITY 2400.0
    YOUNG_MODULUS 26000000000.0
    POISSON_RATIO 0.2
    THERMAL_EXPANSION 1.0e-5
    CONDUCTIVITY 1.0
    SPECIFIC_HEAT 1000.0
    CONSTITUTIVE_LAW_NAME ThermalLinearElastic2DPlaneStrain
End Properties

Begin Nodes
1 0.0 0.0 0.0
2 1.0 0.0 0.0
3 1.0 1.0 0.0
End Nodes

Begin Elements SmallDisplacementThermoMechanicElement2D3N
1 1 1 2 3
End Elements

Begin Elements EulerianConvDiff2D
2 1 1 2 3
End Elements

Begin SubModelPart Parts_Concrete
    Begin SubModelPartNodes
    1
    2
    3
    End SubModelPartNodes
    Begin SubModelPartElements
    1
    End SubModelPartElements
End SubModelPart

Begin SubModelPart Thermal_Concrete
    Begin SubModelPartNodes
    1
    2
    3
    End SubModelPartNodes
    Begin SubModelPartElements
    2
    End SubModelPartElements
End SubModelPart
"""

# Previous-results mesh used by the initial-state transfer test. The imposed
# displacement is a uniaxial stretch u_x = 0.01 * x, so the plane-strain stress
# produced by ThermalLinearElastic2DPlaneStrain is known analytically. The
# nodal stress tensor is deliberately different from that stress so the test can
# distinguish the constitutive response from the post-processing field.
_PREVIOUS_MDPA = """Begin ModelPartData
End ModelPartData

Begin Properties 0
End Properties

Begin Properties 1
    DENSITY 2400.0
    YOUNG_MODULUS 26000000000.0
    POISSON_RATIO 0.2
    THERMAL_EXPANSION 1.0e-5
    CONDUCTIVITY 1.0
    SPECIFIC_HEAT 1000.0
    CONSTITUTIVE_LAW_NAME ThermalLinearElastic2DPlaneStrain
End Properties

Begin Nodes
1 0.0 0.0 0.0
2 1.0 0.0 0.0
3 1.0 1.0 0.0
End Nodes

Begin Elements SmallDisplacementThermoMechanicElement2D3N
1 1 1 2 3
End Elements

Begin NodalData DISPLACEMENT
1 0 [3](0.0,0.0,0.0)
2 0 [3](0.01,0.0,0.0)
3 0 [3](0.01,0.0,0.0)
End NodalData

Begin NodalData NODAL_CAUCHY_STRESS_TENSOR
1 0 [2,2]((1.0e8,0.0),(0.0,0.0))
2 0 [2,2]((1.0e8,0.0),(0.0,0.0))
3 0 [2,2]((1.0e8,0.0),(0.0,0.0))
End NodalData
"""

# The analysis is run in a child interpreter so that a regression to the
# initialization bug (a SIGSEGV) is reported as a test failure instead of
# killing the whole test suite.
_CHILD = r"""import json
import os

import KratosMultiphysics as KM
from KratosMultiphysics.DamApplication.dam_analysis import DamAnalysis

with open(os.path.join(os.getcwd(), "parameters.json"), "r") as parameter_file:
    parameters = KM.Parameters(parameter_file.read())

DamAnalysis(KM.Model(), parameters).Run()
"""

_EXPECTED_UNIAXIAL_STRESS = (
    288888888.8888889,  # sigma_xx
    72222222.22222222,  # sigma_yy
    0.0,  # sigma_zz
    0.0,  # sigma_xy
    0.0,  # sigma_yz
    0.0,  # sigma_xz
)
_IMPORTED_NODAL_STRESS_XX = 1.0e8


class TestDamThermoMechanicInitialOutput(KratosUnittest.TestCase):
    """Regression tests for the initial Gauss-point output of the Dam solver.

    The dual thermal/mechanical solver must initialize its strategies before
    the initial result frame is written. Otherwise the elements still have an
    empty constitutive-law vector when the GiD writer requests
    CAUCHY_STRESS_TENSOR at t=0.
    """

    def _base_parameters(self, problem_name, add_previous_results=False):
        return {
            "problem_data": {
                "problem_name": problem_name,
                "model_part_name": "MainModelPart",
                "domain_size": 2,
                "parallel_type": "OpenMP",
                "number_of_threads": 1,
                "start_time": 0.0,
                # No time loop: the test only needs the initial frame written
                # before the first solution step.
                "end_time": 0.0,
                "time_step": 1.0,
                "time_scale": "Seconds",
                "consider_construction": False,
                "streamlines_utility": False
            },
            "solver_settings": {
                "solver_type": "dam_thermo_mechanic_solver",
                "model_import_settings": {
                    "input_type": "mdpa",
                    "input_filename": "mini",
                    "input_file_label": 0
                },
                "echo_level": 0,
                "buffer_size": 2,
                "processes_sub_model_part_list": [],
                "thermal_solver_settings": {
                    "echo_level": 0,
                    "reform_dofs_at_each_step": True,
                    "clear_storage": False,
                    "compute_reactions": False,
                    "move_mesh_flag": True,
                    "compute_norm_dx_flag": False,
                    "theta_scheme": 1.0,
                    "block_builder": True,
                    "linear_solver_settings": {
                        "solver_type": "bicgstab",
                        "max_iteration": 200,
                        "tolerance": 1.0e-7,
                        "preconditioner_type": "ilu0",
                        "scaling": False
                    },
                    "problem_domain_sub_model_part_list": ["Thermal_Concrete"],
                    "thermal_loads_sub_model_part_list": []
                },
                "mechanical_solver_settings": {
                    "solution_type": "Quasi-Static",
                    "strategy_type": "Newton-Raphson",
                    "scheme_type": "Newmark",
                    "convergence_criterion": "And_criterion",
                    "displacement_relative_tolerance": 1.0e-4,
                    "displacement_absolute_tolerance": 1.0e-9,
                    "residual_relative_tolerance": 1.0e-4,
                    "residual_absolute_tolerance": 1.0e-9,
                    "max_iteration": 10,
                    "echo_level": 0,
                    "buffer_size": 2,
                    "compute_reactions": True,
                    "reform_dofs_at_each_step": False,
                    "move_mesh_flag": True,
                    "block_builder": True,
                    "clear_storage": False,
                    "rayleigh_m": 0.0,
                    "rayleigh_k": 0.0,
                    "nonlocal_damage": False,
                    "linear_solver_settings": {
                        "solver_type": "bicgstab",
                        "max_iteration": 200,
                        "tolerance": 1.0e-7,
                        "preconditioner_type": "ilu0",
                        "scaling": False
                    },
                    "problem_domain_sub_model_part_list": ["Parts_Concrete"],
                    "body_domain_sub_model_part_list": ["Parts_Concrete"],
                    "mechanical_loads_sub_model_part_list": [],
                    "loads_sub_model_part_list": [],
                    "loads_variable_list": []
                }
            },
            "output_configuration": {
                "result_file_configuration": {
                    "gidpost_flags": {
                        "GiDPostMode": "GiD_PostAscii",
                        "WriteDeformedMeshFlag": "WriteDeformed",
                        "WriteConditionsFlag": "WriteElementsOnly",
                        "MultiFileFlag": "SingleFile"
                    },
                    "output_control_type": "time_s",
                    "output_interval": 1.0,
                    "start_output_results": 0,
                    "body_output": True,
                    "node_output": False,
                    "skin_output": False,
                    "plane_output": [],
                    "nodal_results": ["DISPLACEMENT", "TEMPERATURE"],
                    "gauss_point_results": ["CAUCHY_STRESS_TENSOR"]
                },
                "point_data_configuration": []
            },
            "constraints_process_list": [],
            "loads_process_list": [],
            "construction_process": [],
            "transfer_results_process": {
                "save_intermediate_mechanical_variables": False,
                "save_intermediate_thermal_variables": False,
                "save_intermediate_variables_step": 0,
                "save_final_mechanical_variables": False,
                "save_final_thermal_variables": False,
                "add_previous_results": add_previous_results,
                "type_of_results": "Mechanical",
                "add_displacement": True,
                "add_stress": True,
                "add_temperature": True,
                "add_reference_temperature": True,
                "file_name_mechanical": "mini_mechanical.mdpa",
                "file_name_thermal": "mini_thermal.mdpa"
            },
            "temperature_by_device_list": [],
            "output_device_list": []
        }

    def _run_analysis(self, parameters, with_previous_results=False):
        with tempfile.TemporaryDirectory() as working_directory:
            with open(os.path.join(working_directory, "mini.mdpa"), "w") as mdpa_file:
                mdpa_file.write(_MDPA)
            if with_previous_results:
                with open(os.path.join(working_directory, "mini_mechanical.mdpa"), "w") as mdpa_file:
                    mdpa_file.write(_PREVIOUS_MDPA)
            with open(os.path.join(working_directory, "parameters.json"), "w") as parameter_file:
                json.dump(parameters, parameter_file, indent=4)
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
                msg="Dam thermo-mechanical analysis failed (exit {}):\n{}".format(
                    process.returncode, output))

            result_file = os.path.join(
                working_directory, parameters["problem_data"]["problem_name"] + ".post.res")
            self.assertTrue(os.path.isfile(result_file), "GiD result file was not written")
            gauss_point_values = self._read_gauss_point_frame(result_file)
            nodal_results = parameters["output_configuration"][
                "result_file_configuration"]["nodal_results"]
            nodal_stress = None
            if "NODAL_CAUCHY_STRESS_TENSOR" in nodal_results:
                nodal_stress = self._read_nodal_values(
                    result_file, "NODAL_CAUCHY_STRESS_TENSOR")
            return gauss_point_values, nodal_stress

    def _read_gauss_point_frame(self, result_file):
        with open(result_file, "r") as file:
            lines = [line.strip() for line in file]

        header = 'Result "CAUCHY_STRESS_TENSOR" "Kratos" 0 Matrix OnGaussPoints "tri1_element_gp"'
        self.assertIn(header, lines)
        header_index = lines.index(header)
        self.assertEqual(lines[header_index + 1], "Values")

        values = []
        for line in lines[header_index + 2:]:
            if line == "End Values":
                break
            parts = line.split()
            if parts:
                values.append([float(value) for value in parts[1:]])
        return values

    def _read_nodal_values(self, result_file, variable_name):
        with open(result_file, "r") as file:
            lines = [line.strip() for line in file]

        header = 'Result "{}" "Kratos" 0 Matrix OnNodes'.format(variable_name)
        self.assertIn(header, lines)
        header_index = lines.index(header)
        self.assertEqual(lines[header_index + 1], "Values")

        values = {}
        for line in lines[header_index + 2:]:
            if line == "End Values":
                break
            parts = line.split()
            if parts:
                values[int(parts[0])] = [float(value) for value in parts[1:]]
        return values

    def test_gauss_point_output_at_t0_without_previous_results(self):
        parameters = self._base_parameters("initial_gp_output")
        values, _ = self._run_analysis(parameters)

        # A single triangle has one Gauss point. The constitutive law must have
        # been initialized, and the undeformed state must produce zero stress.
        # GiD stores the tensor as six components.
        self.assertEqual(len(values), 1)
        self.assertEqual(len(values[0]), 6)
        for value in values[0]:
            self.assertLess(abs(value), 1.0e-9)

    def test_initial_stress_transfer_drives_constitutive_response(self):
        parameters = self._base_parameters(
            "initial_stress_transfer", add_previous_results=True)
        parameters["output_configuration"]["result_file_configuration"][
            "nodal_results"].append("NODAL_CAUCHY_STRESS_TENSOR")

        values, nodal_stress = self._run_analysis(parameters, with_previous_results=True)
        self.assertIsNotNone(nodal_stress)

        # The transferred displacement must reach the constitutive law: the
        # Gauss-point stress must equal the analytically expected plane-strain
        # stress of u_x = 0.01 * x, and not the imported nodal tensor. GiD
        # stores the tensor as six components (xx, yy, zz, xy, yz, xz).
        self.assertEqual(len(values), 1)
        self.assertEqual(len(values[0]), 6)
        for computed, expected in zip(values[0], _EXPECTED_UNIAXIAL_STRESS):
            # The ASCII GiD output keeps six significant digits.
            self.assertAlmostEqual(computed, expected, delta=abs(expected) * 1.0e-5)

        for node_id in (1, 2, 3):
            self.assertAlmostEqual(
                nodal_stress[node_id][0], _IMPORTED_NODAL_STRESS_XX,
                delta=abs(_IMPORTED_NODAL_STRESS_XX) * 1.0e-12)
            self.assertAlmostEqual(
                nodal_stress[node_id][1], 0.0, delta=1.0e-12)


if __name__ == "__main__":
    KratosUnittest.main()
