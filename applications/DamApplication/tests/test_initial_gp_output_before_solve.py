import json
import os
import subprocess
import tempfile

import KratosMultiphysics as KM
import KratosMultiphysics.DamApplication  # noqa: F401
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.testing.utilities import GetPython3Command

# A minimal single-quadrilateral continuum mdpa. The historical Dam name
# SmallDisplacementThermoMechanicElement2D4N resolves to the
# StructuralMechanicsApplication SmallDisplacement element on this branch.
_MDPA = """Begin ModelPartData
End ModelPartData
Begin Properties 1
DENSITY 1000.0
YOUNG_MODULUS 1000000000.0
POISSON_RATIO 0.2
THICKNESS 1.0
CONSTITUTIVE_LAW_NAME LinearElasticPlaneStrain2DLaw
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

Begin SubModelPart Parts_Parts_Auto1
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
"""

# The child reproducer drives the real DamMechanicalSolver lifecycle and, after
# DamMechanicalSolver.Initialize() and BEFORE any solve, queries a constitutive
# Gauss-point quantity. It prints the query result for the parent to validate.
_CHILD = r"""
import json
import os

import KratosMultiphysics as KM
from KratosMultiphysics.DamApplication import dam_mechanical_solver

workdir = os.getcwd()

settings = KM.Parameters(json.dumps({
    "solver_type": "dam_mechanical_solver",
    "model_import_settings": {"input_type": "mdpa", "input_filename": os.path.join(workdir, "q4"), "input_file_label": 0},
    "echo_level": 0,
    "buffer_size": 2,
    "processes_sub_model_part_list": [""],
    "mechanical_solver_settings": {
        "echo_level": 0,
        "reform_dofs_at_each_step": False,
        "clear_storage": False,
        "compute_reactions": True,
        "move_mesh_flag": True,
        "solution_type": "Quasi-Static",
        "scheme_type": "Newmark",
        "rayleigh_m": 0.0,
        "rayleigh_k": 0.0,
        "nonlocal_damage": False,
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
        "problem_domain_sub_model_part_list": ["Parts_Parts_Auto1"],
        "body_domain_sub_model_part_list": ["Parts_Parts_Auto1"],
        "mechanical_loads_sub_model_part_list": [],
        "loads_sub_model_part_list": [],
        "loads_variable_list": []
    }
}))

model = KM.Model()
main_model_part = model.CreateModelPart("Main")
main_model_part.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 2)

solver = dam_mechanical_solver.DamMechanicalSolver(main_model_part, settings)
solver.AddVariables()
solver.ImportModelPart()
solver.AddDofs()
solver.Initialize()

# Initial Gauss-point constitutive output, intentionally before the first solve.
element = main_model_part.GetElement(1)
stress = element.CalculateOnIntegrationPoints(
    KM.CAUCHY_STRESS_VECTOR, main_model_part.ProcessInfo)

print("NGP", len(stress))
for entry in stress:
    print("STRS", " ".join("{:.10e}".format(v) for v in entry))
print("DONE")
"""


class TestInitialGaussPointOutputBeforeSolve(KratosUnittest.TestCase):

    def test_gp_cauchy_output_before_first_solve(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            mdpa_path = os.path.join(temp_dir, "q4.mdpa")
            child_path = os.path.join(temp_dir, "child_reproducer.py")
            with open(mdpa_path, "w") as f:
                f.write(_MDPA)
            with open(child_path, "w") as f:
                f.write(_CHILD)

            proc = subprocess.Popen(
                [GetPython3Command(), child_path],
                cwd=temp_dir,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
            )
            p_out = proc.communicate(timeout=120)
            self.assertEqual(
                proc.returncode, 0,
                msg="Child solver-lifecycle reproducer failed (exit {}):\n{}".format(
                    proc.returncode, p_out[0].decode('ascii', 'replace')))

        gathered = self.__parse_gp_values(p_out[0].decode('ascii'))
        self.assertEqual(len(gathered), 4)  # 2x2 Gauss points on the Q4
        for values in gathered:
            self.assertEqual(len(values), 3)  # plane-strain stress vector
            for value in values:
                self.assertTrue(value == value, "non-finite stress")  # NaN check
                self.assertLess(abs(value), 1.0e-9)  # zero initial stress

    @staticmethod
    def __parse_gp_values(output):
        gathered = []
        for line in output.splitlines():
            parts = line.split()
            if len(parts) != 4 or parts[0] != "STRS":
                continue
            gathered.append([float(parts[1]), float(parts[2]), float(parts[3])])
        return gathered


if __name__ == "__main__":
    KratosUnittest.main()
