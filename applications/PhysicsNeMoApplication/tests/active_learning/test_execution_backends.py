import os
import sys
from pathlib import Path

import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as KratosUtilities

from KratosMultiphysics.PhysicsNeMoApplication.active_learning.sample_io import KratosALSample
from KratosMultiphysics.PhysicsNeMoApplication.active_learning.execution_backends.in_process_backend import InProcessBackend
from KratosMultiphysics.PhysicsNeMoApplication.active_learning.execution_backends.subprocess_backend import SubprocessBackend

_AUX_FILES = Path(__file__).parent / "aux_files"
sys.path.insert(0, str(_AUX_FILES))


class TestInProcessBackend(KratosUnittest.TestCase):
    def setUp(self):
        self.working_directory = Path("test_in_process_cases")
        self.parameters_file = Path("test_in_process_parameters.json")
        with open(_AUX_FILES / "template" / "ProjectParameters.json") as f:
            self.parameters_file.write_text(f.read())

    def tearDown(self):
        KratosUtilities.DeleteDirectoryIfExisting(str(self.working_directory))
        KratosUtilities.DeleteFileIfExisting(str(self.parameters_file))

    def _CreateBackend(self):
        return InProcessBackend(Kratos.Parameters("""{
            "project_parameters_file" : "test_in_process_parameters.json",
            "analysis_stage_module"   : "dummy_analysis",
            "working_directory"       : "test_in_process_cases",
            "model_part_name"         : "Main",
            "output_field_specs"      : [ { "variable_name" : "PRESSURE", "data_location" : "node_historical" } ]
        }"""))

    def test_RunCase(self):
        backend = self._CreateBackend()
        self.assertFalse(backend.is_external)

        sample = KratosALSample("case_a", parameters={"dummy_settings/alpha": 2.0})
        labeled = backend.RunCase(sample)

        self.assertTrue(labeled.is_labeled)
        pressure = labeled.fields["PRESSURE__node_historical"]
        # PRESSURE = alpha * x with x = 0..4 and alpha = 2
        self.assertTrue(numpy.allclose(pressure, 2.0 * numpy.arange(5.0)))

    def test_RunCaseRestoresCwd(self):
        backend = self._CreateBackend()
        cwd = os.getcwd()
        backend.RunCase(KratosALSample("case_b", parameters={}))
        self.assertEqual(os.getcwd(), cwd)


class TestSubprocessBackend(KratosUnittest.TestCase):
    def setUp(self):
        self.working_directory = Path("test_subprocess_cases")
        # Make the dummy_analysis fixture importable inside the subprocess.
        self._previous_pythonpath = os.environ.get("PYTHONPATH", "")
        os.environ["PYTHONPATH"] = str(_AUX_FILES) + os.pathsep + self._previous_pythonpath

    def tearDown(self):
        os.environ["PYTHONPATH"] = self._previous_pythonpath
        KratosUtilities.DeleteDirectoryIfExisting(str(self.working_directory))

    def _CreateBackend(self, run_command=None):
        settings = Kratos.Parameters("""{
            "template_directory" : "",
            "working_directory"  : "test_subprocess_cases",
            "timeout_seconds"    : 120,
            "max_retries"        : 1
        }""")
        settings["template_directory"].SetString(str(_AUX_FILES / "template"))
        if run_command is not None:
            settings.AddEmptyArray("run_command")
            for entry in run_command:
                settings["run_command"].Append(entry)
        return SubprocessBackend(settings)

    def test_RunCase(self):
        backend = self._CreateBackend([sys.executable, "MainKratos.py"])
        self.assertTrue(backend.is_external)

        labeled = backend.RunCase(KratosALSample("case_a", parameters={"dummy_settings/alpha": 3.0}))

        self.assertTrue(numpy.allclose(
            labeled.fields["PRESSURE__node_historical"], 3.0 * numpy.arange(5.0)))
        self.assertEqual(labeled.metadata["attempts"], 1)

    def test_CaseIsolation(self):
        backend = self._CreateBackend([sys.executable, "MainKratos.py"])
        backend.RunCase(KratosALSample("case_one", parameters={"dummy_settings/alpha": 1.0}))
        backend.RunCase(KratosALSample("case_two", parameters={"dummy_settings/alpha": 2.0}))
        self.assertTrue((self.working_directory / "case_one" / "ProjectParameters.json").is_file())
        self.assertTrue((self.working_directory / "case_two" / "ProjectParameters.json").is_file())

    def test_FailingCommandRaisesAfterRetries(self):
        backend = self._CreateBackend([sys.executable, "-c", "raise SystemExit(1)"])
        with self.assertRaisesRegex(RuntimeError, "failed after 2 attempt"):
            backend.RunCase(KratosALSample("case_fail", parameters={}))

    def test_MissingTemplateRaises(self):
        settings = Kratos.Parameters("""{ "template_directory" : "does_not_exist_dir" }""")
        with self.assertRaises(FileNotFoundError):
            SubprocessBackend(settings)


if __name__ == '__main__':
    KratosUnittest.main()
