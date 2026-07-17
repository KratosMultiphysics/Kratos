# Import necessary modules from KratosMultiphysics
import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.meshio_output_process import Factory as MeshioOutputProcessFactory

# Import the pathlib and tempfile modules
from pathlib import Path
import tempfile


def _PopulateModelPart(model_part):
    """Creates 4 nodes and 1 tetrahedron element."""
    props = model_part.CreateNewProperties(1)
    model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
    model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
    model_part.CreateNewNode(3, 0.0, 1.0, 0.0)
    model_part.CreateNewNode(4, 0.0, 0.0, 1.0)
    model_part.CreateNewElement("Element3D4N", 1, [1, 2, 3, 4], props)


class TestMeshioOutputProcess(KratosUnittest.TestCase):
    """Tests for the meshio++ output process."""

    def setUp(self):
        self.model = KratosMultiphysics.Model()
        self.model_part = self.model.CreateModelPart("main")
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)
        _PopulateModelPart(self.model_part)
        # The OutputController requires the control variable in the ProcessInfo
        # already at Check (the AnalysisStage initializes them the same way)
        self.model_part.ProcessInfo[KratosMultiphysics.STEP] = 0
        self.model_part.ProcessInfo[KratosMultiphysics.TIME] = 0.0

    def _RunSolutionLoop(self, process, number_of_steps):
        process.Check()
        process.ExecuteInitialize()
        process.ExecuteBeforeSolutionLoop()
        for step in range(1, number_of_steps + 1):
            self.model_part.ProcessInfo[KratosMultiphysics.STEP] = step
            self.model_part.ProcessInfo[KratosMultiphysics.TIME] = 0.1 * step
            for node in self.model_part.Nodes:
                node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, step * node.Id)
            process.ExecuteInitializeSolutionStep()
            process.ExecuteFinalizeSolutionStep()
            if process.IsOutputStep():
                process.PrintOutput()
        process.ExecuteFinalize()

    def testVtuFileSeries(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            settings = KratosMultiphysics.Parameters("""{
                "Parameters" : {
                    "model_part_name" : "main",
                    "output_name"     : "results.vtu",
                    "output_path"     : "%s"
                }
            }""" % str(Path(temp_dir) / "vtu_output").replace("\\\\", "/"))
            process = MeshioOutputProcessFactory(settings, self.model)
            self._RunSolutionLoop(process, 3)

            output_files = sorted(f.name for f in (Path(temp_dir) / "vtu_output").iterdir())
            self.assertEqual(output_files, ["results_1.vtu", "results_2.vtu", "results_3.vtu"])

    def testXdmfTimeSeries(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            settings = KratosMultiphysics.Parameters("""{
                "Parameters" : {
                    "model_part_name"                    : "main",
                    "output_name"                        : "results.xdmf",
                    "output_path"                        : "%s",
                    "output_control_type"                : "time",
                    "output_interval"                    : 0.1,
                    "xdmf_data_format"                   : "XML",
                    "nodal_solution_step_data_variables" : ["TEMPERATURE"]
                }
            }""" % str(Path(temp_dir) / "xdmf_output").replace("\\\\", "/"))
            process = MeshioOutputProcessFactory(settings, self.model)
            self._RunSolutionLoop(process, 3)

            xdmf_file = Path(temp_dir) / "xdmf_output" / "results.xdmf"
            self.assertTrue(xdmf_file.is_file())
            content = xdmf_file.read_text()
            self.assertEqual(content.count("<Time Value="), 3)
            self.assertEqual(content.count('Name="mesh"'), 1)
            self.assertEqual(content.count('Name="TEMPERATURE"'), 3)

    def testExtendedFeaturesThroughProcess(self):
        """Cell data, flags, ids, prefix/postfix and sub model part output driven by the process."""
        volume = self.model_part.CreateSubModelPart("Volume")
        volume.AddNodes([1, 2, 3, 4])
        volume.AddElements([1])
        for element in self.model_part.Elements:
            element.SetValue(KratosMultiphysics.PRESSURE, 2.0 * element.Id)

        with tempfile.TemporaryDirectory() as temp_dir:
            settings = KratosMultiphysics.Parameters("""{
                "Parameters" : {
                    "model_part_name"              : "main",
                    "output_name"                  : "results.vtu",
                    "output_path"                  : "%s",
                    "file_format"                  : "ascii",
                    "custom_name_prefix"           : "pre_",
                    "output_sub_model_parts"       : true,
                    "write_ids"                    : true,
                    "nodal_flags"                  : ["TO_ERASE"],
                    "element_data_value_variables" : ["PRESSURE"]
                }
            }""" % str(Path(temp_dir) / "extended").replace("\\\\", "/"))
            process = MeshioOutputProcessFactory(settings, self.model)
            self._RunSolutionLoop(process, 1)

            produced = sorted(f.name for f in (Path(temp_dir) / "extended").iterdir())
            self.assertEqual(produced, ["pre_results_1.vtu", "pre_results_main_Volume_1.vtu"])
            content = (Path(temp_dir) / "extended" / "pre_results_1.vtu").read_text()
            for token in ("KRATOS_NODE_ID", "KRATOS_ELEMENT_ID", "PROPERTIES_ID", "TO_ERASE", "PRESSURE"):
                self.assertIn(token, content)
            self.assertIn('format="ascii"', content)

    def testMissingExtensionRaises(self):
        settings = KratosMultiphysics.Parameters("""{
            "Parameters" : {
                "model_part_name" : "main",
                "output_name"     : "results_without_extension"
            }
        }""")
        with self.assertRaisesRegex(Exception, "file extension"):
            MeshioOutputProcessFactory(settings, self.model)

    def testProcessIsRegistered(self):
        self.assertTrue(KratosMultiphysics.Registry.HasItem("Processes.KratosMultiphysics.MeshioOutputProcess"))


if __name__ == '__main__':
    KratosUnittest.main()
