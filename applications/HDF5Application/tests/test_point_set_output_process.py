import KratosMultiphysics
import KratosMultiphysics.kratos_utilities as KratosUtils
from KratosMultiphysics import KratosUnittest as UnitTest
from KratosMultiphysics.HDF5Application.point_set_output_process import Factory as PointSetOuputProcessFactory


class TestPointSetOutputProcess(UnitTest.TestCase):

    def setUp(self):
        self.file_name = "test_point_set_output.h5"
        #KratosUtils.DeleteFileIfExisting(self.file_name)


    def tearDown(self):
        #KratosUtils.DeleteFileIfExisting(self.file_name)
        pass


    def test_PointSetOutputProcess(self):
        number_of_elements = 3
        edge_length = 1.0

        model, model_part = self.MakeModel(
            number_of_elements = number_of_elements,
            edge_length = edge_length)

        parameters = KratosMultiphysics.Parameters("""{
            "model_part_name" : "main",
            "positions" : [],
            "output_variables" : ["DISPLACEMENT", "REACTION"],
            "file_parameters"  : {
                "file_name"    : ""
            }
        }""")
        parameters["file_parameters"]["file_name"].SetString(self.file_name)

        # Create sample points across the the element row
        positions = parameters["positions"]
        for i_point in range(2 * (number_of_elements+1) - 1):
            #position = KratosMultiphysics.Vector([i_point*edge_length/2.0, edge_length/2.0, 0.0])
            position = KratosMultiphysics.Vector([1, 2, 3])
            print(position)
            positions.Append(position)
            print(positions)

        point_set_output_process = PointSetOuputProcessFactory(parameters, model)
        point_set_output_process.ExecuteInitialize()
        point_set_output_process.ExecuteFinalizeSolutionStep()


    @staticmethod
    def MakeModel(number_of_elements = 100,
                  edge_length = 1.0):
        """Generate a horizontal row of square elements in the xy plane."""

        model = KratosMultiphysics.Model()
        model_part = model.CreateModelPart("main")
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)

        # Generate nodes
        for bottom_index in range(number_of_elements + 1):
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
        for element_index in range(number_of_elements):
            model_part.CreateNewElement(
                "Element2D4N",
                element_index,
                [element_index + 1, element_index + 2, element_index + number_of_elements + 2, element_index + number_of_elements + 1],
                properties)

        return model, model_part





if __name__ == "__main__":
    UnitTest.main()