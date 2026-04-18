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
import test_helper


NodalValue = Union[float, List[float]]
NodalSeries = Dict[float, Dict[int, NodalValue]]

class KratosGeoMechanicsUDSMSerializationTest(KratosUnittest.TestCase):

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
        result = test_helper.get_nodal_variable_from_ascii(str(filepath), variable_name)

        if not result:
            raise RuntimeError(f"Variable '{variable_name}' not found in {filepath}")

        return result

    def _max_common_time(self, lhs: NodalSeries, rhs: NodalSeries) -> float:
        common = set(lhs.keys()).intersection(rhs.keys())
        if not common:
            raise RuntimeError("No common time steps found between result files")
        return max(common)

    @staticmethod
    def _get_nodes_mismatch_error(lhs_nodes: Dict[int, NodalValue],
                                  rhs_nodes: Dict[int, NodalValue],
                                  time: float) -> Union[str, None]:
        if set(lhs_nodes.keys()) == set(rhs_nodes.keys()):
            return None

        only_lhs = sorted(set(lhs_nodes.keys()) - set(rhs_nodes.keys()))
        only_rhs = sorted(set(rhs_nodes.keys()) - set(lhs_nodes.keys()))
        return (
            f"Node id mismatch at time {time}: only in lhs={only_lhs[:10]}, "
            f"only in rhs={only_rhs[:10]}"
        )

    @staticmethod
    def _get_tolerance(lhs_value: float, rhs_value: float, abs_tol: float, rel_tol: float) -> float:
        return max(abs_tol, rel_tol * max(abs(lhs_value), abs(rhs_value), 1.0))

    def _compare_scalar_value(self,
                              node_id: int,
                              lhs_value: float,
                              rhs_value: float,
                              abs_tol: float,
                              rel_tol: float) -> List[str]:
        diff = abs(lhs_value - rhs_value)
        tol = self._get_tolerance(lhs_value, rhs_value, abs_tol, rel_tol)
        if diff > tol:
            return [
                f"Node {node_id}: lhs={lhs_value}, rhs={rhs_value}, diff={diff}, tol={tol}"
            ]
        return []

    def _compare_vector_value(self,
                              node_id: int,
                              lhs_value: List[float],
                              rhs_value: List[float],
                              abs_tol: float,
                              rel_tol: float) -> List[str]:
        if len(lhs_value) != len(rhs_value):
            return [
                f"Component count mismatch at node {node_id}: lhs={lhs_value}, rhs={rhs_value}"
            ]

        errors: List[str] = []
        for comp, (lv, rv) in enumerate(zip(lhs_value, rhs_value)):
            diff = abs(lv - rv)
            tol = self._get_tolerance(lv, rv, abs_tol, rel_tol)
            if diff > tol:
                errors.append(
                    f"Node {node_id}, comp {comp}: lhs={lv}, rhs={rv}, diff={diff}, tol={tol}"
                )
        return errors

    def _compare_node_value(self,
                            node_id: int,
                            lhs_value: NodalValue,
                            rhs_value: NodalValue,
                            abs_tol: float,
                            rel_tol: float) -> List[str]:
        if isinstance(lhs_value, list) != isinstance(rhs_value, list):
            return [f"Type mismatch at node {node_id}: lhs={lhs_value}, rhs={rhs_value}"]

        if isinstance(lhs_value, list):
            return self._compare_vector_value(node_id, lhs_value, rhs_value, abs_tol, rel_tol)

        return self._compare_scalar_value(node_id, lhs_value, rhs_value, abs_tol, rel_tol)

    def _compare_nodal_series_at_time(
        self,
        lhs: NodalSeries,
        rhs: NodalSeries,
        time: float,
        abs_tol: float,
        rel_tol: float,
    ) -> List[str]:
        lhs_nodes = lhs[time]
        rhs_nodes = rhs[time]

        mismatch_error = self._get_nodes_mismatch_error(lhs_nodes, rhs_nodes, time)
        if mismatch_error is not None:
            return [mismatch_error]

        errors: List[str] = []
        for node_id in sorted(lhs_nodes.keys()):
            errors.extend(
                self._compare_node_value(
                    node_id, lhs_nodes[node_id], rhs_nodes[node_id], abs_tol, rel_tol
                )
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
