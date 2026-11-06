"""Subprocess execution backend — the primary labeling mode.

Each sample gets an isolated copy of a case template directory; its
ProjectParameters.json is patched with the sample's parameter overrides and
the case is solved by an external command (a plain python run by default; an
"srun"/"sbatch --wait" prefix turns it into an HPC job submission with no
extra code). Results travel back through the .npz files written by this
application's DatasetExportProcess, which the case template is expected to
list among its processes.

Running externally keeps Kratos's MPI ranks and any torch.distributed ranks
of the surrounding training loop in separate OS processes — this is the
reason this backend is the recommended one, not just an HPC convenience.

Pure Python: no torch/physicsnemo imports.
"""

import concurrent.futures
import shutil
import subprocess
import sys
from pathlib import Path

import KratosMultiphysics as Kratos
from KratosMultiphysics.PhysicsNeMoApplication.active_learning.execution_backends.base_backend import KratosExecutionBackend
from KratosMultiphysics.PhysicsNeMoApplication.active_learning.sample_io import (
    KratosALSample, ApplyParameterOverrides, LoadFieldsFromNpzDirectory)


class SubprocessBackend(KratosExecutionBackend):
    """Runs each sample's case as an external subprocess."""

    def __init__(self, settings: Kratos.Parameters) -> None:
        default_settings = Kratos.Parameters("""{
            "template_directory"      : "PLEASE_SPECIFY_TEMPLATE_DIRECTORY",
            "run_command"             : [],
            "project_parameters_file" : "ProjectParameters.json",
            "working_directory"       : "physics_nemo_al_cases",
            "results_directory"       : "physics_nemo_dataset",
            "timeout_seconds"         : 3600,
            "max_retries"             : 1,
            "max_parallel_jobs"       : 1
        }""")
        settings.ValidateAndAssignDefaults(default_settings)

        self.template_directory = Path(settings["template_directory"].GetString())
        self.run_command = settings["run_command"].GetStringArray()
        if not self.run_command:
            self.run_command = [sys.executable, "MainKratos.py"]
        self.project_parameters_file = settings["project_parameters_file"].GetString()
        self.working_directory = Path(settings["working_directory"].GetString())
        self.results_directory = settings["results_directory"].GetString()
        self.timeout_seconds = settings["timeout_seconds"].GetInt()
        self.max_retries = settings["max_retries"].GetInt()
        self.max_parallel_jobs = settings["max_parallel_jobs"].GetInt()
        if self.max_parallel_jobs < 1:
            raise ValueError(
                f"\"max_parallel_jobs\" must be >= 1 [ max_parallel_jobs = {self.max_parallel_jobs} ].")

        if not self.template_directory.is_dir():
            raise FileNotFoundError(f"Template directory \"{self.template_directory}\" does not exist.")

    @property
    def is_external(self) -> bool:
        return True

    def RunCases(self, samples) -> list:
        """Labels a batch of samples, fanning out up to max_parallel_jobs
        concurrent subprocesses (each sample already gets its own isolated
        case directory, so parallel runs cannot collide). Threads suffice:
        subprocess.run releases the GIL while the external solve runs.
        """
        samples = list(samples)
        if self.max_parallel_jobs == 1 or len(samples) <= 1:
            return super().RunCases(samples)

        def run_one(sample):
            try:
                return self.RunCase(sample)
            except Exception as error:  # noqa: BLE001 - reported to the caller
                return error

        with concurrent.futures.ThreadPoolExecutor(max_workers=self.max_parallel_jobs) as executor:
            return list(executor.map(run_one, samples))

    def RunCase(self, sample: KratosALSample) -> KratosALSample:
        case_directory = self.working_directory / sample.sample_id
        if case_directory.exists():
            shutil.rmtree(case_directory)
        shutil.copytree(self.template_directory, case_directory)

        parameters_path = case_directory / self.project_parameters_file
        with open(parameters_path, "r") as f:
            case_parameters = Kratos.Parameters(f.read())
        ApplyParameterOverrides(case_parameters, sample.parameters)
        with open(parameters_path, "w") as f:
            f.write(case_parameters.PrettyPrintJsonString())

        attempts = 0
        while True:
            attempts += 1
            result = subprocess.run(
                self.run_command,
                cwd=case_directory,
                capture_output=True,
                text=True,
                timeout=self.timeout_seconds)
            if result.returncode == 0:
                break
            if attempts > self.max_retries:
                raise RuntimeError(
                    f"Case \"{sample.sample_id}\" failed after {attempts} attempt(s) "
                    f"(exit code {result.returncode}). Last stderr:\n{result.stderr[-2000:]}")
            Kratos.Logger.PrintWarning(
                "SubprocessBackend",
                f"Case \"{sample.sample_id}\" attempt {attempts} failed (exit code "
                f"{result.returncode}); retrying.")

        fields, metadata = LoadFieldsFromNpzDirectory(case_directory / self.results_directory)
        sample.fields.update(fields)
        sample.metadata.update({
            "case_directory": str(case_directory),
            "attempts": attempts,
            "TIME": metadata["TIME"],
            "STEP": metadata["STEP"],
        })
        return sample
