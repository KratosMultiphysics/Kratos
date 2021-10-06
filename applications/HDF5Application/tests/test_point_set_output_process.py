import KratosMultiphysics
import KratosMultiphysics.kratos_utilities as KratosUtils
from KratosMultiphysics import KratosUnittest as UnitTest
from KratosMultiphysics.HDF5Application.point_set_output_process import Factory as PointSetOutputProcessFactory

import math
import pathlib


class TestPointSetOutputProcess(UnitTest.TestCase):

    communicator = KratosMultiphysics.ParallelEnvironment.GetDefaultDataCommunicator()
    file_name = "test_point_set_output.h5"


    def setUp(self):
        if self.communicator.Rank() == 1:
            KratosUtils.DeleteFileIfExisting(self.file_name)
        self.communicator.Barrier()


    def tearDown(self):
        if self.communicator.Rank() == 1:
            #KratosUtils.DeleteFileIfExisting(self.file_name)
            pass


    def test_PointSetOutputProcess(self):
        number_of_elements = 100
        edge_length = 1.0

        model, model_part = self.MakeModel(
            number_of_elements = number_of_elements,
            edge_length = edge_length)

        parameters = KratosMultiphysics.Parameters("""{
            "model_part_name"    : "main",
            "positions"          : [],
            "output_variables"   : ["DISPLACEMENT", "REACTION"],
            "coordinates_prefix" : "/test_point_set_output_<model_part_name>",
            "variables_prefix"   : "/test_point_set_output_<model_part_name>/test_step_<step>",
            "file_parameters"    : {
                "file_name"      : ""
            }
        }""")
        parameters["file_parameters"]["file_name"].SetString(self.file_name)

        # Create sample points across the the element row
        positions = parameters["positions"]
        for i_point in range(2 * (number_of_elements+1) - 1):
            position = KratosMultiphysics.Vector([i_point*edge_length/2.0, edge_length/2.0, 0.0])
            positions.Append(position)

        # Write coordinates and variables
        point_set_output_process = PointSetOutputProcessFactory(parameters, model)
        point_set_output_process.ExecuteInitialize()
        point_set_output_process.ExecuteFinalizeSolutionStep()

#        # Open the results' file
#        file_parameters = parameters["file_parameters"].Clone()
#        file_parameters["file_access_mode"].SetString("read_only")
#        if 1 < self.communicator.Size():
#            HDF5File = KratosMultiphysics.HDF5Application.HDF5FileParallel
#            file_parameters.AddString("file_driver", "mpio")
#        else:
#            HDF5File = KratosMultiphysics.HDF5Application.HDF5FileSerial
#
#        file = HDF5File(file_parameters)
#
#        # Check paths
#        self.assertTrue(file.IsGroup("test_point_set_output_main"))
#        self.assertTrue(file.IsDataSet("test_point_set_output_main/POSITION"))


    @staticmethod
    def MakeModel(number_of_elements = 100,
                  edge_length = 1.0):
        """Generate a horizontal row of square elements in the xy plane."""

        model = KratosMultiphysics.Model()
        model_part = model.CreateModelPart("main")
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)

        # Linear node/element partition for mpi (shifted for debugging purposes)
        communicator = KratosMultiphysics.ParallelEnvironment.GetDefaultDataCommunicator()
        number_of_ranks = communicator.Size()
        rank = communicator.Rank()

        chunk_size = math.ceil(number_of_elements / float(number_of_ranks))
        element_index_range = [(rank_id * chunk_size, min((rank_id+1) * chunk_size, number_of_elements)) for rank_id in range(number_of_ranks)][(rank+1) % number_of_ranks]

        # Generate nodes
        for bottom_index in range(element_index_range[0], element_index_range[1] + 1):
            bottom_node = model_part.CreateNewNode(bottom_index + 1,
                                                   bottom_index * edge_length,
                                                   0.0,
                                                   0.0)
            bottom_node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT, [bottom_index, 0.0, 0.0])
            bottom_node.SetSolutionStepValue(KratosMultiphysics.REACTION, [0.0, 0.0, bottom_index])

            top_index = bottom_index + number_of_elements + 1
            top_node = model_part.CreateNewNode(top_index + 1,
                                                bottom_index * edge_length,
                                                edge_length,
                                                0.0)
            top_node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT, [top_index, 0.0, 0.0])
            top_node.SetSolutionStepValue(KratosMultiphysics.REACTION, [0.0, 0.0, top_index])

        # Generate elements
        properties = model_part.GetProperties()[1]
        for element_index in range(*element_index_range):
            model_part.CreateNewElement(
                "Element2D4N",
                element_index,
                [element_index + 1, element_index + 2, element_index + number_of_elements + 3, element_index + number_of_elements + 2],
                properties)

        return model, model_part





if __name__ == "__main__":
    UnitTest.main()