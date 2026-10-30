import json
import os
import subprocess
import tempfile

import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.testing.utilities import GetPython3Command
import test_dam_initial_stress_handling as initial_stress_tests


_MDPA = """Begin ModelPartData
End ModelPartData
Begin Properties 0
End Properties
Begin Properties 1
DENSITY 2.0
YOUNG_MODULUS 10000.0
POISSON_RATIO 0.2
THICKNESS 1.0
CONSTITUTIVE_LAW_NAME LinearElasticPlaneStrain2DLaw
End Properties
Begin Nodes
1 0 0 0
2 1 0 0
3 0 1 0
4 1 1 0
5 0 2 0
6 1 2 0
End Nodes
Begin Elements SmallDisplacementElement2D4N
1 1 1 2 4 3
2 1 3 4 6 5
End Elements
Begin NodalData DISPLACEMENT_X
1 1 0.0
2 1 0.0
End NodalData
Begin NodalData DISPLACEMENT_Y
1 1 0.0
2 1 0.0
End NodalData
Begin SubModelPart Parts_Concrete
Begin SubModelPartNodes
1
2
3
4
5
6
End SubModelPartNodes
Begin SubModelPartElements
1
2
End SubModelPartElements
End SubModelPart
"""

_CHILD = r"""
import json
import os
import sys
import KratosMultiphysics as KM
import KratosMultiphysics.DamApplication as Dam

# Prefer source scripts in a source checkout; installed test packages use the
# installed application normally.
if os.path.isdir(sys.argv[1]):
    Dam.__path__.insert(0, sys.argv[1])
from KratosMultiphysics.DamApplication.dam_analysis import DamAnalysis

with open('parameters.json') as f:
    parameters = KM.Parameters(f.read())
simulation = DamAnalysis(KM.Model(), parameters)
simulation.Initialize()
part = simulation.main_model_part

def snapshot():
    residual = [0.] * (2 * part.NumberOfNodes())
    stresses = []
    for element in part.Elements:
        stresses.extend([list(v) for v in element.CalculateOnIntegrationPoints(
            KM.CAUCHY_STRESS_VECTOR, part.ProcessInfo)])
        lhs, rhs = KM.Matrix(), KM.Vector()
        element.CalculateLocalSystem(lhs, rhs, part.ProcessInfo)
        for i, node in enumerate(element.GetGeometry()):
            for j in range(2):
                residual[2 * (node.Id - 1) + j] += rhs[2 * i + j]
    return {
        'time': part.ProcessInfo[KM.TIME], 'step': part.ProcessInfo[KM.STEP],
        'displacement': [list(n.GetSolutionStepValue(KM.DISPLACEMENT)) for n in part.Nodes],
        'previous': [list(n.GetSolutionStepValue(KM.DISPLACEMENT, 1)) for n in part.Nodes],
        'reaction': [list(n.GetSolutionStepValue(KM.REACTION)) for n in part.Nodes],
        'stress': stresses, 'residual': residual,
        'fixity': [[n.IsFixed(KM.DISPLACEMENT_X), n.IsFixed(KM.DISPLACEMENT_Y)] for n in part.Nodes],
    }

initial = snapshot()
steps = []
solve = simulation.solver.Solve
def recording_solve():
    solve()
    steps.append(snapshot())
simulation.solver.Solve = recording_solve
simulation.RunMainTemporalLoop()
simulation.Finalize()
with open('snapshots.json', 'w') as f:
    json.dump({'initial': initial, 'steps': steps}, f)
"""


class TestDamStandardRestart(KratosUnittest.TestCase):
    """Real DamAnalysis checkpoint and continuation for a gravity-loaded Q4 body."""

    @staticmethod
    def _load(value, interval):
        return {
            "python_module": "assign_vector_variable_process",
            "kratos_module": "KratosMultiphysics",
            "Parameters": {
                "model_part_name": "MainModelPart.Parts_Concrete",
                "variable_name": "VOLUME_ACCELERATION", "constrained": False,
                "interval": interval, "value": value,
            },
        }

    def _parameters(self, name, end_time, model_input, checkpoint):
        settings = initial_stress_tests.TestDamInitialStressHandling._base_parameters(name, False, end_time)
        settings["solver_settings"]["model_import_settings"] = model_input
        settings["solver_settings"]["processes_sub_model_part_list"] = []
        settings["transfer_results_process"]["add_previous_results"] = False
        settings["constraints_process_list"] = []
        settings["loads_process_list"] = []
        if model_input["input_type"] == "mdpa":
            settings["loads_process_list"].append(self._load([0., -10., 0.], [0., "End"]))
        settings["loads_process_list"].append(self._load([2., None, None], [3., "End"]))
        if checkpoint:
            settings["output_processes"] = {"restart": [{
                "python_module": "save_restart_process", "kratos_module": "KratosMultiphysics",
                "Parameters": {"model_part_name": "MainModelPart", "restart_save_frequency": 0.0,
                               "output_path": "checkpoints", "save_restart_files_in_folder": True},
            }]}
        return settings

    def _run(self, directory, parameters, expect_success=True):
        os.makedirs(directory, exist_ok=True)
        with open(os.path.join(directory, "mini.mdpa"), "w") as f:
            f.write(_MDPA)
        with open(os.path.join(directory, "parameters.json"), "w") as f:
            json.dump(parameters, f)
        with open(os.path.join(directory, "run.py"), "w") as f:
            f.write(_CHILD)
        source_scripts = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "python_scripts"))
        result = subprocess.run([GetPython3Command(), "run.py", source_scripts], cwd=directory,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=120)
        if not expect_success:
            self.assertNotEqual(result.returncode, 0, result.stdout)
            return result.stdout
        self.assertEqual(result.returncode, 0, result.stdout)
        with open(os.path.join(directory, "snapshots.json")) as f:
            return json.load(f)

    def _compare(self, actual, expected):
        if isinstance(expected, list):
            self.assertEqual(len(actual), len(expected))
            for a, e in zip(actual, expected):
                self._compare(a, e)
        else:
            self.assertAlmostEqual(actual, expected, delta=1.e-10 + 1.e-10 * abs(expected))

    @staticmethod
    def _max_difference(actual, expected):
        if isinstance(expected, list):
            return max((TestDamStandardRestart._max_difference(a, e)
                        for a, e in zip(actual, expected)), default=0.)
        return abs(actual - expected)

    def test_mechanical_restart_continuation(self):
        with tempfile.TemporaryDirectory() as root:
            source = {"input_type": "mdpa", "input_filename": "mini", "input_file_label": 0}
            reference = self._run(os.path.join(root, "reference"),
                                  self._parameters("reference", 3, source, True))
            first_dir = os.path.join(root, "first")
            first = self._run(first_dir, self._parameters("original", 1, source, True))
            checkpoint = os.path.join(first_dir, "checkpoints", "MainModelPart_1.0.rest")
            self.assertTrue(os.path.isfile(checkpoint))
            original_results = os.path.join(first_dir, "original.post.res")
            with open(original_results, "rb") as f:
                original_bytes = f.read()
            restart_input = {"input_type": "rest", "input_filename": os.path.join(first_dir, "MainModelPart"),
                             "input_output_path": "checkpoints", "restart_load_file_label": "1.0"}
            resumed = self._run(os.path.join(root, "restart"),
                                self._parameters("continued", 3, restart_input, False))
            with open(os.path.join(root, "restart", "continued.post.res")) as f:
                self.assertIn('Result "CAUCHY_STRESS_TENSOR" "Kratos" 1 Matrix OnGaussPoints', f.read())
            with open(original_results, "rb") as f:
                self.assertEqual(f.read(), original_bytes)
            self.assertEqual([s["time"] for s in reference["steps"]], [1., 2., 3.])
            self.assertEqual([s["time"] for s in resumed["steps"]], [2., 3.])
            self.assertEqual([s["step"] for s in resumed["steps"]], [2, 3])
            self.assertEqual(resumed["initial"]["time"], 1.)
            self.assertEqual(resumed["initial"]["step"], 1)
            for key in ("displacement", "previous", "stress", "reaction", "residual", "fixity"):
                self._compare(resumed["initial"][key], first["steps"][-1][key])
            for actual, expected in zip(resumed["steps"], reference["steps"][1:]):
                for key in ("displacement", "stress", "reaction", "residual"):
                    self._compare(actual[key], expected[key])
                for residual in actual["residual"][4:]:
                    self.assertAlmostEqual(residual, 0., delta=1.e-10)
            for key in ("displacement", "stress", "reaction", "residual"):
                error = max(self._max_difference(a[key], e[key])
                            for a, e in zip(resumed["steps"], reference["steps"][1:]))
                print("restart maximum {} error: {:.6e}".format(key, error))
            self._compare(resumed["steps"][0]["displacement"], resumed["initial"]["displacement"])
            unsupported = self._parameters("unsupported", 3, restart_input, False)
            unsupported["constraints_process_list"] = [{"python_module": "unknown_process", "Parameters": {}}]
            output = self._run(os.path.join(root, "unsupported"), unsupported, False)
            self.assertIn("does not support processes in constraints_process_list", output)
            unsupported = self._parameters("nonlocal", 3, restart_input, False)
            unsupported["solver_settings"]["mechanical_solver_settings"]["nonlocal_damage"] = True
            output = self._run(os.path.join(root, "nonlocal"), unsupported, False)
            self.assertIn("supports only 2D quasistatic local elasticity", output)
            output = self._run(first_dir, self._parameters("original", 3, restart_input, False), False)
            self.assertIn("separate directory without existing GiD files", output)
            with open(original_results, "rb") as f:
                self.assertEqual(f.read(), original_bytes)


if __name__ == "__main__":
    KratosUnittest.main()
