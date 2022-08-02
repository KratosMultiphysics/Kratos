# Core imports
import KratosMultiphysics
import KratosMultiphysics.kratos_utilities as KratosUtils
from KratosMultiphysics import KratosUnittest as UnitTest
from KratosMultiphysics.testing.utilities import ReadModelPart

# HDF5 imports
from KratosMultiphysics import HDF5Application as HDF5
from KratosMultiphysics.HDF5Application.point_set_output_process import Factory as PointSetOutputProcessFactory

# STL imports
import math
import pathlib


class TestPointSetOutputProcess(UnitTest.TestCase):

    communicator = KratosMultiphysics.Testing.GetDefaultDataCommunicator()
    file_name = "test_point_set_output.h5"


    def setUp(self):
        KratosUtils.DeleteFileIfExisting(self.file_name)
        self.communicator.Barrier()


    def tearDown(self):
        # The output file is not actually checked yet in the script,
        # so if you need to validate the results, comment the line
        # below.
        self.communicator.Barrier()
        KratosUtils.DeleteFileIfExisting(self.file_name)


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

            # Print output if necessary
            if point_set_output_process.IsOutputStep():
                point_set_output_process.PrintOutput()

        self.communicator.Barrier()

        # Open output file
        file_parameters = parameters["file_parameters"].Clone()
        file_parameters.AddString("file_access_mode","read_only")
        if self.communicator.IsDistributed():
            File = KratosMultiphysics.HDF5Application.HDF5FileParallel
            file_parameters.AddString("file_driver", "mpio")
        else:
            File = KratosMultiphysics.HDF5Application.HDF5FileSerial

        file = File(file_parameters)

        # Check output file structure
        root = "/test_point_set_output_{}".format(parameters["model_part_name"].GetString())
        self.assertTrue(file.IsGroup(root))
        self.assertTrue(file.IsDataSet(root + "/POSITION"))


    @property
    def parameters(self) -> KratosMultiphysics.Parameters:
        parameters = KratosMultiphysics.Parameters("""{
            "model_part_name"      : "main",
            "positions"            : [],
            "output_variables"     : ["DISPLACEMENT_X", "VELOCITY"],
            "output_frequency"     : 3,
            "coordinates_prefix"   : "/test_point_set_output_<model_part_name>",
            "variables_prefix"     : "/test_point_set_output_<model_part_name>/test_step_<step>",
            "file_parameters"      : {
                "file_name"        : ""
            }
        }""")
        parameters["file_parameters"]["file_name"].SetString(self.file_name)

        # Create sample points across the the element row
        positions = parameters["positions"]
        number_of_points = 50
        for i_point in range(number_of_points):
            s = float(i_point) / (number_of_points - 1)
            position = KratosMultiphysics.Vector([s, s, 0.0])
            positions.Append(position)

        return parameters


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