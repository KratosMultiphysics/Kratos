# Core imports
import KratosMultiphysics
import KratosMultiphysics.kratos_utilities as KratosUtils
from KratosMultiphysics import KratosUnittest as UnitTest
from KratosMultiphysics.testing.utilities import ReadModelPart

# HDF5 imports
from KratosMultiphysics import HDF5Application as HDF5
from KratosMultiphysics.HDF5Application.point_set_output_process import Factory as PointSetOutputProcessFactory
from KratosMultiphysics.HDF5Application.core.file_io import OpenHDF5File

# STL imports
import pathlib


class TestPointSetOutputProcess(UnitTest.TestCase):

    communicator = KratosMultiphysics.Testing.GetDefaultDataCommunicator()

    def setUp(self):
        KratosUtils.DeleteFileIfExisting("test_point_set_output.h5")
        self.communicator.Barrier()


    def tearDown(self):
        # The output file is not actually checked yet in the script,
        # so if you need to validate the results, comment the line
        # below.
        self.communicator.Barrier()
        KratosUtils.DeleteFileIfExisting("test_point_set_output.h5")


    def test_PointSetOutputProcessWrite(self):
        model, model_part = self.MakeModel()
        parameters = self.parameters
        number_of_steps = 10

        # Write coordinates and variables
        process_parameters = KratosMultiphysics.Parameters()
        process_parameters.AddValue("Parameters", parameters)
        point_set_output_process = PointSetOutputProcessFactory(process_parameters, model)
        point_set_output_process.ExecuteInitialize()

        for i_step in range(number_of_steps):
            # Create new step data
            model_part.CloneTimeStep(2.0 * i_step)
            model_part.ProcessInfo[KratosMultiphysics.STEP] = i_step

            # Modify variables
            for node in model_part.Nodes:
                for variable, increment in zip((KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.VELOCITY), (1.0, [0.0,0.0,1.0])):
                    node.SetSolutionStepValue(variable,node.GetSolutionStepValue(variable) + increment)
                    node.SetValue(variable,node.GetValue(variable) + increment)

            # Print output if necessary
            if point_set_output_process.IsOutputStep():
                point_set_output_process.PrintOutput()

        self.communicator.Barrier()

        # Open output file
        file_parameters = parameters["file_settings"].Clone()
        file_parameters["file_access_mode"].SetString("read_only")
        with OpenHDF5File(file_parameters, model_part) as file:
            # Check output file structure
            root = "/test_point_set_output_{}".format(parameters["model_part_name"].GetString())
            self.assertTrue(file.IsGroup(root))
            self.assertTrue(file.IsDataSet(root + "/Ids"))
            self.assertTrue(file.IsDataSet(root + "/Coordinates"))

            self.__CheckStep(file, root, 3)
            self.__CheckStep(file, root, 6)
            self.__CheckStep(file, root, 9)

    @property
    def parameters(self) -> KratosMultiphysics.Parameters:
        parameters = KratosMultiphysics.Parameters("""{
            "model_part_name"      : "main",
            "file_settings": {
                "file_name"        : "test_point_set_output.h5",
                "time_format"      : "0.4f",
                "file_access_mode" : "read_write",
                "max_files_to_keep": "unlimited",
                "echo_level"       :  0
            },
            "output_time_settings": {
                "output_interval": 3.0
            },
            "point_output_settings": {
                "prefix"              : "/test_point_set_output_<model_part_name>",
                "time_format"         : "0.4f",
                "positions"           : [[0.0, 0.0, 0.0]],
                "search_configuration": "initial",
                "search_tolerance"    : 1e-6
            },
            "nodal_solution_step_data_settings": {
                "prefix"           : "/test_point_set_output_<model_part_name>/test_step_hist_<step>/",
                "list_of_variables": ["DISPLACEMENT_X", "VELOCITY"],
                "time_format"      : "0.4f"
            },
            "nodal_data_value_settings": {
                "prefix"           : "/test_point_set_output_<model_part_name>/test_step_nonhist_<step>/",
                "list_of_variables": ["DISPLACEMENT_X", "VELOCITY"],
                "time_format"      : "0.4f"
            }
        }""")

        return parameters

    def __CheckStep(self, file: HDF5.HDF5File, root: str, step: int) -> None:
        self.assertTrue(file.IsDataSet(f"{root}/test_step_hist_{step}/DISPLACEMENT_X"))
        self.assertTrue(file.IsDataSet(f"{root}/test_step_hist_{step}/VELOCITY"))
        self.assertTrue(file.GetDataDimensions(f"{root}/test_step_hist_{step}/DISPLACEMENT_X"), [1])
        self.assertTrue(file.GetDataDimensions(f"{root}/test_step_hist_{step}/VELOCITY"), [1, 3])
        self.assertTrue(file.IsDataSet(f"{root}/test_step_nonhist_{step}/DISPLACEMENT_X"))
        self.assertTrue(file.IsDataSet(f"{root}/test_step_nonhist_{step}/VELOCITY"))
        self.assertTrue(file.GetDataDimensions(f"{root}/test_step_nonhist_{step}/DISPLACEMENT_X"), [1])
        self.assertTrue(file.GetDataDimensions(f"{root}/test_step_nonhist_{step}/VELOCITY"), [1, 3])


    @staticmethod
    def MakeModel():
        model = KratosMultiphysics.Model()
        model_part = model.CreateModelPart("main")
        model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = 3

        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT_X)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)

        ReadModelPart(str(TestPointSetOutputProcess.GetInputMDPAPath()), model_part)

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