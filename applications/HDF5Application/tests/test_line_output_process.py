# Core imports
import KratosMultiphysics
import KratosMultiphysics.kratos_utilities as KratosUtils
from KratosMultiphysics import KratosUnittest as UnitTest
from KratosMultiphysics.testing.utilities import ReadModelPart

# HDF5 imports
import KratosMultiphysics.HDF5Application as HDF5
from KratosMultiphysics.HDF5Application.line_output_process import Factory as LineOutputProcessFactory
from KratosMultiphysics.HDF5Application.core.file_io import OpenHDF5File

# STD imports
import pathlib

class TestLineOutputProcess(UnitTest.TestCase):

    communicator = KratosMultiphysics.Testing.GetDefaultDataCommunicator()
    file_name = "test_line_output.h5"

    def setUp(self):
        KratosUtils.DeleteFileIfExisting(self.file_name)
        self.communicator.Barrier()

    def tearDown(self):
        # The output file is not actually checked yet in the script,
        # so if you need to validate the results, comment the line
        # below.
        self.communicator.Barrier()
        KratosUtils.DeleteFileIfExisting(self.file_name)

    def test_LineOutputProcess(self):
        model, model_part = self.MakeModel()
        parameters = self.parameters
        number_of_steps = 10

        process_parameters = KratosMultiphysics.Parameters()
        process_parameters.AddValue("Parameters", parameters)
        process = LineOutputProcessFactory(process_parameters, model)
        process.ExecuteInitialize()

        for i_step in range(number_of_steps):
            # Create new step data
            model_part.CloneTimeStep(2.0 * i_step)
            model_part.ProcessInfo[KratosMultiphysics.STEP] = i_step

            # Modify variables
            for node in model_part.Nodes:
                for variable, increment in zip((KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.VELOCITY), (1.0, [0.0,0.0,1.0])):
                    node.SetSolutionStepValue(variable,node.GetSolutionStepValue(variable) + increment)
                    node.SetValue(variable,node.GetValue(variable) + increment)

            # Print output if requested
            if process.IsOutputStep():
                process.PrintOutput()

        self.communicator.Barrier()

        # Open output file
        file_parameters = parameters["file_settings"].Clone()
        file_parameters["file_access_mode"].SetString("read_only")
        with OpenHDF5File(file_parameters, model_part) as file:
            # Check output file structure
            root = "/test_line_output_{}".format(parameters["model_part_name"].GetString())
            self.assertTrue(file.IsGroup(root))
            self.assertTrue(file.IsDataSet(root + "/Ids"))
            self.assertTrue(file.IsDataSet(root + "/Coordinates"))

            self.__CheckStep(file, root, 3)
            self.__CheckStep(file, root, 6)
            self.__CheckStep(file, root, 9)

    def __CheckStep(self, file: HDF5.HDF5File, root: str, step: int) -> None:
        self.assertTrue(file.IsDataSet(f"{root}/test_step_hist_{step}/DISPLACEMENT_X"))
        self.assertTrue(file.IsDataSet(f"{root}/test_step_hist_{step}/VELOCITY"))
        self.assertTrue(file.GetDataDimensions(f"{root}/test_step_hist_{step}/DISPLACEMENT_X"), [51])
        self.assertTrue(file.GetDataDimensions(f"{root}/test_step_hist_{step}/VELOCITY"), [51, 3])
        self.assertTrue(file.IsDataSet(f"{root}/test_step_nonhist_{step}/DISPLACEMENT_X"))
        self.assertTrue(file.IsDataSet(f"{root}/test_step_nonhist_{step}/VELOCITY"))
        self.assertTrue(file.GetDataDimensions(f"{root}/test_step_nonhist_{step}/DISPLACEMENT_X"), [51])
        self.assertTrue(file.GetDataDimensions(f"{root}/test_step_nonhist_{step}/VELOCITY"), [51, 3])

    @property
    def parameters(self) -> KratosMultiphysics.Parameters:
        parameters = KratosMultiphysics.Parameters("""{
            "model_part_name"      : "main",
            "file_settings": {
                "file_name"        : "",
                "time_format"      : "0.4f",
                "file_access_mode" : "read_write",
                "max_files_to_keep": "unlimited",
                "echo_level"       :  0
            },
            "output_time_settings": {
                "output_interval": 3.0
            },
            "line_output_settings": {
                "prefix"              : "/test_line_output_<model_part_name>",
                "time_format"         : "0.4f",
                "start_point"         : [0.0, 0.0, 0.0],
                "end_point"           : [1.0, 0.0, 0.0],
                "number_of_points"    : 51,
                "search_configuration": "initial",
                "search_tolerance"    : 1e-6
            },
            "nodal_solution_step_data_settings": {
                "prefix"           : "/test_line_output_<model_part_name>/test_step_hist_<step>/",
                "list_of_variables": ["DISPLACEMENT_X", "VELOCITY"],
                "time_format"      : "0.4f"
            },
            "nodal_data_value_settings": {
                "prefix"           : "/test_line_output_<model_part_name>/test_step_nonhist_<step>/",
                "list_of_variables": ["DISPLACEMENT_X", "VELOCITY"],
                "time_format"      : "0.4f"
            }
        }""")
        parameters["file_settings"]["file_name"].SetString(self.file_name)
        return parameters

    @staticmethod
    def MakeModel():
        model = KratosMultiphysics.Model()
        model_part = model.CreateModelPart("main")
        model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = 3

        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT_X)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)

        ReadModelPart(str(TestLineOutputProcess.GetInputMDPAPath()), model_part)

        for node in model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X, node.X)
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, [node.X, node.Y, node.Z])
            node.SetValue(KratosMultiphysics.DISPLACEMENT_X, node.X * 2)
            node.SetValue(KratosMultiphysics.VELOCITY, [node.X * 2, node.Y * 2, node.Z * 2])

        return model, model_part

    @staticmethod
    def GetInputMDPAPath() -> pathlib.Path:
        script_directory      = pathlib.Path(__file__).absolute().parent
        kratos_root_directory = script_directory.parent.parent.parent
        test_input_directory  = kratos_root_directory / "kratos" / "tests" / "auxiliar_files_for_python_unittest"
        test_file_stem        = test_input_directory / "mdpa_files" / "test_processes"
        test_file_path        = pathlib.Path(str(test_file_stem) + ".mdpa")

        if not test_file_path.is_file():
            raise FileNotFoundError("Test file not found: {}".format(test_file_path))

        return test_file_stem

if __name__ == "__main__":
    UnitTest.main()