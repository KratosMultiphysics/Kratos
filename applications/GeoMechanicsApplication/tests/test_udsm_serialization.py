import importlib
import copy
import glob
import json
import os
from pathlib import Path
import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.project import Project
from helper_utilities import _compare_case_outputs


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

    def test_udsm_serialization(self):
        script_dir = str(Path(__file__).resolve().parent)
        case_dir = os.path.join(script_dir, "UDSM_serialization")
        base_project_settings = self._read_project_settings(
            os.path.join(case_dir, "orchestrator_stages.json")
        )

        cwd = os.getcwd()

        try:
            os.chdir(case_dir)

            save_checkpoint_settings = self._make_save_checkpoint_settings(
                base_project_settings
            )
            self._run_orchestrator_with_settings(save_checkpoint_settings)

            checkpoint_dir = os.path.join(case_dir, "checkpoints")
            if not os.path.exists(checkpoint_dir) or not any(
                glob.glob(os.path.join(checkpoint_dir, "stage_2*"))
            ):
                raise FileNotFoundError(
                    "Stage-2 checkpoint dump was not created in 'checkpoints'."
                )

            restart_settings = self._make_restart_from_checkpoint_settings(
                base_project_settings
            )
            self._run_orchestrator_with_settings(restart_settings)

            full_run_stage3 = os.path.join(case_dir, "stage3.post.res")
            checkpoint_stage3 = os.path.join(case_dir, "stage3_from_checkpoint.post.res")

            if not os.path.exists(full_run_stage3):
                raise FileNotFoundError(
                    f"Full-run stage-3 result not found: {full_run_stage3}"
                )
            if not os.path.exists(checkpoint_stage3):
                raise FileNotFoundError(
                    f"Checkpoint-run stage-3 result not found: {checkpoint_stage3}"
                )
            abs_tol = 1e-10
            rel_tol = 1e-8
            variables = ["TOTAL_DISPLACEMENT", "WATER_PRESSURE"]
            _compare_case_outputs(
                full_run_stage3,
                checkpoint_stage3,
                abs_tol,
                rel_tol,
                variables,
            )

        finally:
            os.chdir(cwd)


if __name__ == "__main__":
    KratosUnittest.main()
