import KratosMultiphysics
import KratosMultiphysics.kratos_utilities as KratosUtils
from KratosMultiphysics import KratosUnittest as UnitTest
from KratosMultiphysics.HDF5Application.line_output_process import Factory as LineOutputProcessFactory
from KratosMultiphysics.testing.utilities import ReadModelPart

import math
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
        KratosUtils.DeleteFileIfExisting(self.file_name)


    def test_LineOutputProcess(self):
        number_of_elements = 100
        edge_length = 1.0

        model, model_part = self.MakeModel()

        parameters = KratosMultiphysics.Parameters("""{
            "model_part_name"   : "main",
            "start_point"       : [0.0, 0.5, 0.0],
            "end_point"         : [100.0, 0.5, 0.0],
            "number_of_points"  : 201,
            "output_variables"  : ["DISPLACEMENT", "REACTION"],
            "coordinates_prefix" : "/test_line_output_<model_part_name>",
            "variables_prefix"   : "/test_line_output_<model_part_name>/test_step_<step>",
            "file_parameters"   : {
                "file_name"     : ""
            }
        }""")
        parameters["file_parameters"]["file_name"].SetString(self.file_name)

        # Write coordinates and variables
        point_set_output_process = LineOutputProcessFactory(parameters, model)
        point_set_output_process.ExecuteInitialize()
        point_set_output_process.PrintOutput()


    @staticmethod
    def MakeModel():
        model = KratosMultiphysics.Model()
        model_part = model.CreateModelPart("main")
        model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = 3

        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT_X)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)

        ReadModelPart("/home/mate/build/kratos/Debug/install/kratos/tests/auxiliar_files_for_python_unittest/mdpa_files/test_processes", model_part)

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