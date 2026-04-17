import importlib
import copy
import glob
import json
import os
from pathlib import Path
from typing import Dict, List, Union
import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.project import Project


NodalValue = Union[float, List[float]]
NodalSeries = Dict[float, Dict[int, NodalValue]]

class KratosGeoMechanicsUDSMSerializationTest(KratosUnittest.TestCase):

    def _cleanup_case_run_artifacts(self, case_dir: str) -> None:
        checkpoint_dir = os.path.join(case_dir, "checkpoints")
        if os.path.exists(checkpoint_dir):
            for path in glob.glob(os.path.join(checkpoint_dir, "*")):
                if os.path.isfile(path):
                    os.remove(path)

        result_patterns = [
            "stage1.post.*",
            "stage2.post.*",
            "stage3.post.*",
            "stage3_from_checkpoint.post.*",
            "mesh.post.*",
        ]
        for pattern in result_patterns:
            for path in glob.glob(os.path.join(case_dir, pattern)):
                if os.path.isfile(path):
                    os.remove(path)

    def _read_project_settings(self, filepath: Path) -> dict:
        with open(filepath, "r", encoding="utf-8") as f:
            return json.load(f)

    def _run_orchestrator_with_settings(self, project_settings: dict) -> None:
        project_parameters = Kratos.Parameters(json.dumps(project_settings))
        project = Project(project_parameters)

        orchestrator_reg_entry = Kratos.Registry[
            project.GetSettings()["orchestrator"]["name"].GetString()
        ]
        orchestrator_module = importlib.import_module(
            orchestrator_reg_entry["ModuleName"]
        )
        orchestrator_class = getattr(
            orchestrator_module, orchestrator_reg_entry["ClassName"]
        )
        orchestrator_instance = orchestrator_class(project)
        orchestrator_instance.Run()

    def _make_save_checkpoint_settings(self, base_settings: dict) -> dict:
        settings = copy.deepcopy(base_settings)
        settings["orchestrator"]["settings"]["execution_list"] = ["stage_1", "stage_2", "stage_3"]
        settings["orchestrator"]["settings"]["load_from_checkpoint"] = None
        settings["orchestrator"]["settings"]["stage_checkpoints"] = True
        settings["orchestrator"]["settings"]["stage_checkpoints_folder"] = "checkpoints"
        return settings

    def _make_restart_from_checkpoint_settings(self, base_settings: dict) -> dict:
        settings = copy.deepcopy(base_settings)
        settings["orchestrator"]["settings"]["execution_list"] = ["stage_3"]
        settings["orchestrator"]["settings"][
            "load_from_checkpoint"
        ] = "checkpoints/stage_2"
        settings["orchestrator"]["settings"]["stage_checkpoints"] = False

        # Avoid overwriting the stage-3 result from the full run.
        settings["stages"]["stage_3"]["stage_settings"]["output_processes"][
            "gid_output"
        ][0]["Parameters"]["output_name"] = "stage3_from_checkpoint"
        return settings

    def _parse_gid_ascii_nodal_variable(
        self, filepath: Path, variable_name: str
    ) -> NodalSeries:
        with open(filepath, "r", encoding="utf-8") as f:
            lines = f.readlines()

        result: NodalSeries = {}
        active = False
        current_time = None

        for line in lines:
            stripped = line.strip()

            if active and stripped == "End Values":
                active = False
                current_time = None
                continue

            if active:
                if stripped == "Values" or not stripped:
                    continue
                parts = stripped.split()
                node_id = int(parts[0])
                values = [float(v) for v in parts[1:]]
                if len(values) == 1:
                    result[current_time][node_id] = values[0]
                else:
                    result[current_time][node_id] = values
                continue

            if f'"{variable_name}"' in stripped:
                parts = stripped.split()
                if len(parts) < 4:
                    raise RuntimeError(f"Could not parse result header line: {line}")
                current_time = float(parts[3])
                result[current_time] = {}
                active = True

        if not result:
            raise RuntimeError(f"Variable '{variable_name}' not found in {filepath}")

        return result

    def _max_common_time(self, lhs: NodalSeries, rhs: NodalSeries) -> float:
        common = set(lhs.keys()).intersection(rhs.keys())
        if not common:
            raise RuntimeError("No common time steps found between result files")
        return max(common)

    def _compare_nodal_series_at_time(
        self,
        lhs: NodalSeries,
        rhs: NodalSeries,
        time: float,
        abs_tol: float,
        rel_tol: float,
    ) -> List[str]:
        errors: List[str] = []

        lhs_nodes = lhs[time]
        rhs_nodes = rhs[time]

        if set(lhs_nodes.keys()) != set(rhs_nodes.keys()):
            only_lhs = sorted(set(lhs_nodes.keys()) - set(rhs_nodes.keys()))
            only_rhs = sorted(set(rhs_nodes.keys()) - set(lhs_nodes.keys()))
            errors.append(
                f"Node id mismatch at time {time}: only in lhs={only_lhs[:10]}, only in rhs={only_rhs[:10]}"
            )
            return errors

        for node_id in sorted(lhs_nodes.keys()):
            lval = lhs_nodes[node_id]
            rval = rhs_nodes[node_id]

            if isinstance(lval, list) != isinstance(rval, list):
                errors.append(
                    f"Type mismatch at node {node_id}: lhs={lval}, rhs={rval}"
                )
                continue

            if isinstance(lval, list):
                if len(lval) != len(rval):
                    errors.append(
                        f"Component count mismatch at node {node_id}: lhs={lval}, rhs={rval}"
                    )
                    continue
                for comp, (lv, rv) in enumerate(zip(lval, rval)):
                    diff = abs(lv - rv)
                    tol = max(abs_tol, rel_tol * max(abs(lv), abs(rv), 1.0))
                    if diff > tol:
                        errors.append(
                            f"Node {node_id}, comp {comp}: lhs={lv}, rhs={rv}, diff={diff}, tol={tol}"
                        )
            else:
                diff = abs(lval - rval)
                tol = max(abs_tol, rel_tol * max(abs(lval), abs(rval), 1.0))
                if diff > tol:
                    errors.append(
                        f"Node {node_id}: lhs={lval}, rhs={rval}, diff={diff}, tol={tol}"
                    )

        return errors

    def _compare_case_outputs(
        self,
        full_stage2_res: Path,
        checkpoint_stage2_res: Path,
        abs_tol: float,
        rel_tol: float,
    ) -> None:
        variables = ["TOTAL_DISPLACEMENT", "WATER_PRESSURE"]

        failures: List[str] = []

        for variable in variables:
            legacy_series = self._parse_gid_ascii_nodal_variable(
                full_stage2_res, variable
            )
            orchestrator_series = self._parse_gid_ascii_nodal_variable(
                checkpoint_stage2_res, variable
            )

            compare_time = self._max_common_time(legacy_series, orchestrator_series)
            errors = self._compare_nodal_series_at_time(
                legacy_series, orchestrator_series, compare_time, abs_tol, rel_tol
            )

            if errors:
                failures.append(
                    f"Variable {variable} mismatch at time {compare_time}. First 20 differences:\n"
                    + "\n".join(errors[:20])
                )

        if failures:
            raise AssertionError("\n\n".join(failures))

    def _resolve_case_dir(self) -> str:
        script_dir = str(Path(__file__).resolve().parent)
        TEST_NAME = "UDSM-serialization"
        case_dir = os.path.join(script_dir, TEST_NAME)
        if not os.path.exists(case_dir):
            raise FileNotFoundError(
                f"Could not find case directory for test '{TEST_NAME}': {case_dir}"
            )
        return case_dir

    def _checkpoint_path_exists(self, case_dir: str) -> bool:
        checkpoint_dir = os.path.join(case_dir, "checkpoints")
        if not os.path.exists(checkpoint_dir):
            return False

        return any(glob.glob(os.path.join(checkpoint_dir, "stage_2*")))

    def test_udsm_serialization(self):
        abs_tol = 1e-10
        rel_tol = 1e-8
        case_dir = self._resolve_case_dir()
        base_project_settings = self._read_project_settings(
            os.path.join(case_dir, "ProjectParameters.json")
        )

        self._cleanup_case_run_artifacts(case_dir)

        print(f"Case directory: {case_dir}")
        print(f"Run dir: {case_dir}")

        cwd = os.getcwd()

        try:
            os.chdir(case_dir)

            print("Running orchestrator with stage checkpoints enabled...")
            save_checkpoint_settings = self._make_save_checkpoint_settings(
                base_project_settings
            )
            self._run_orchestrator_with_settings(save_checkpoint_settings)

            if not self._checkpoint_path_exists(case_dir):
                raise FileNotFoundError(
                    "Stage-2 checkpoint dump was not created in 'checkpoints'."
                )

            print("Running stage 3 from stage-2 checkpoint dump...")
            restart_settings = self._make_restart_from_checkpoint_settings(
                base_project_settings
            )
            self._run_orchestrator_with_settings(restart_settings)

            full_stage3 = os.path.join(case_dir, "stage3.post.res")
            checkpoint_stage3 = os.path.join(case_dir, "stage3_from_checkpoint.post.res")

            if not os.path.exists(full_stage3):
                raise FileNotFoundError(
                    f"Full-run stage-3 result not found: {full_stage3}"
                )
            if not os.path.exists(checkpoint_stage3):
                raise FileNotFoundError(
                    f"Checkpoint-run stage-3 result not found: {checkpoint_stage3}"
                )

            print("Comparing stage-3 results (full run vs checkpoint restart)...")
            self._compare_case_outputs(full_stage3, checkpoint_stage3, abs_tol, rel_tol)

            print(
                "SUCCESS: Stage-3 results match between full run and checkpoint restart."
            )
        finally:
            os.chdir(cwd)


if __name__ == "__main__":
    KratosUnittest.main()
