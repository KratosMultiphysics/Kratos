# Import the necessary Kratos modules for multiphysics simulations and testing.
import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import specific applications for mapping and co-simulation.
import KratosMultiphysics.MappingApplication
import KratosMultiphysics.CoSimulationApplication as CoSimulationApplication

# Import utilities for processing JSON results.
from KratosMultiphysics.from_json_check_result_process import FromJsonCheckResultProcess
from KratosMultiphysics.json_output_process import JsonOutputProcess

# Import basic dependencies for math and file path handling.
import math
from pathlib import Path

# Function to get the full file path for a given file name.
def GetFilePath(fileName):
    return Path(__file__).parent / fileName

class Test3D1DDataTransferProcessBlock(KratosUnittest.TestCase):
    """Test class for the 3D to 1D data transfer process."""

    def setUp(self):
        """Setup function to initialize model parts for testing."""
        self.current_model = KM.Model()

        # Initialize a 3D model part called "Block".
        self.model_part_block = self.current_model.CreateModelPart("Block")
        self.model_part_block.ProcessInfo[KM.STEP] = 0
        self.model_part_block.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)
        self.model_part_block.ProcessInfo.SetValue(KM.TIME, 0.0)
        self.model_part_block.ProcessInfo.SetValue(KM.DELTA_TIME, 1.0)
        self.model_part_block.AddNodalSolutionStepVariable(KM.TEMPERATURE)
        input_mdpa = GetFilePath("3d_1d_data_transfer/3d_1d_data_transfer_mdpa/block")
        model_part_io_block = KM.ModelPartIO(input_mdpa)
        model_part_io_block.ReadModelPart(self.model_part_block)

        # Initialize a 1D model part called "Line".
        self.model_part_line = self.current_model.CreateModelPart("Line")
        self.model_part_line.ProcessInfo[KM.STEP] = 0
        self.model_part_line.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)
        self.model_part_line.ProcessInfo.SetValue(KM.TIME, 0.0)
        self.model_part_line.ProcessInfo.SetValue(KM.DELTA_TIME, 1.0)
        self.model_part_line.AddNodalSolutionStepVariable(KM.TEMPERATURE)

        # Create nodes and elements in the 1D model part.
        number_of_line_elements = 100
        initial_z_coordinate = 0.0
        end_z_coordinate = 5.0
        self.model_part_line.CreateNewNode(1, 0.5, 0.5, initial_z_coordinate)
        prop = self.model_part_line.GetProperties()[1]
        for i in range(number_of_line_elements):
            z = initial_z_coordinate + (i + 1) * (end_z_coordinate - initial_z_coordinate) / number_of_line_elements
            self.model_part_line.CreateNewNode(i + 2, 0.5, 0.5, z)
            self.model_part_line.CreateNewElement("Element3D2N", i + 1, [i + 1, i + 2], prop)

        # Set temperature values for nodes in the block and line model parts.
        for node in self.model_part_block.Nodes:
            node.SetSolutionStepValue(KM.TEMPERATURE, node.Z)  # Set temperature based on Z coordinate.
        del node

        for node in self.model_part_line.Nodes:
            node.SetSolutionStepValue(KM.TEMPERATURE, -node.Z)  # Set negative temperature based on Z coordinate.
        del node

    def test_block_from_1d_to_3d(self):
        """Test data transfer from the 1D line to the 3D block."""
        # Define parameters for the data transfer process.
        parameters = KM.Parameters("""{
            "origin_variables"         : ["TEMPERATURE"],
            "destination_variables"    : ["TEMPERATURE"]
        }""")

        # Execute the data transfer process.
        process = CoSimulationApplication.DataTransfer3D1DProcess(self.model_part_block, self.model_part_line, parameters)
        process.Execute()

        # Define result file path for the output of the transfer.
        result_file = GetFilePath("3d_1d_data_transfer/3d_1d_data_transfer_solutions/block_from_1d_to_3d")

        # # Generate output results in JSON format.
        # generate_result(self.current_model, result_file, "Block")

        # Check the results against expected outcomes.
        check_results(self.current_model, result_file, "Block")

        # Debug output for visualization.
        # debug_vtk(self.current_model, ["Block", "Line"])

    def test_block_from_3d_to_1d(self):
        """Test data transfer from the 3D block to the 1D line."""
        # Define parameters for the data transfer process.
        parameters = KM.Parameters("""{
            "origin_variables"         : ["TEMPERATURE"],
            "destination_variables"    : ["TEMPERATURE"]
        }""")

        # Execute the data transfer process with transpose setting.
        process = CoSimulationApplication.DataTransfer3D1DProcess(self.model_part_block, self.model_part_line, parameters)
        process.Set(KratosMultiphysics.Mapper.USE_TRANSPOSE)
        process.Execute()

        # Define result file path for the output of the transfer.
        result_file = GetFilePath("3d_1d_data_transfer/3d_1d_data_transfer_solutions/block_from_3d_to_1d")

        # # Generate output results in JSON format.
        # generate_result(self.current_model, result_file, "Line")

        # Check the results against expected outcomes.
        check_results(self.current_model, result_file, "Line")

        # Debug output for visualization.
        # debug_vtk(self.current_model, ["Block", "Line"])

# Test class for the toroidal data transfer process.
class Test3D1DDataTransferProcessTorus(KratosUnittest.TestCase):
    def setUp(self):
        """Setup function to initialize model parts for testing with a torus."""
        self.current_model = KM.Model()

        # Initialize a 3D model part called "Torus".
        self.model_part_torus = self.current_model.CreateModelPart("Torus")
        self.model_part_torus.ProcessInfo[KM.STEP] = 0
        self.model_part_torus.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)
        self.model_part_torus.ProcessInfo.SetValue(KM.TIME, 0.0)
        self.model_part_torus.ProcessInfo.SetValue(KM.DELTA_TIME, 1.0)
        self.model_part_torus.AddNodalSolutionStepVariable(KM.TEMPERATURE)
        input_mdpa = GetFilePath("3d_1d_data_transfer/3d_1d_data_transfer_mdpa/torus3d")
        model_part_io_torus = KM.ModelPartIO(input_mdpa)
        model_part_io_torus.ReadModelPart(self.model_part_torus)

        # Initialize a 1D model part called "Circle".
        self.model_part_circle = self.current_model.CreateModelPart("Circle")
        self.model_part_circle.ProcessInfo[KM.STEP] = 0
        self.model_part_circle.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)
        self.model_part_circle.ProcessInfo.SetValue(KM.TIME, 0.0)
        self.model_part_circle.ProcessInfo.SetValue(KM.DELTA_TIME, 1.0)
        self.model_part_circle.AddNodalSolutionStepVariable(KM.TEMPERATURE)

        # Create nodes and elements in the circular model part.
        number_of_circle_elements = 100
        radius = 1.0
        self.model_part_circle.CreateNewNode(1, radius, 0.0, 0)
        prop = self.model_part_circle.GetProperties()[1]
        for i in range(number_of_circle_elements - 1):
            angle = 2 * math.pi * (i + 1) / number_of_circle_elements
            self.model_part_circle.CreateNewNode(i + 2, radius * math.cos(angle), radius * math.sin(angle), 0)
            self.model_part_circle.CreateNewElement("Element3D2N", i + 1, [i + 1, i + 2], prop)
        self.model_part_circle.CreateNewElement("Element3D2N", number_of_circle_elements, [i + 2, 1], prop)

        # Set temperature values for nodes in the torus and circle model parts.
        for node in self.model_part_torus.Nodes:
            node.SetSolutionStepValue(KM.TEMPERATURE, -node.X + node.Y)  # Example temperature calculation.
        del node

        for node in self.model_part_circle.Nodes:
            node.SetSolutionStepValue(KM.TEMPERATURE, node.X - node.Y)  # Example temperature calculation.
        del node

    def test_torus_from_1d_to_3d(self):
        """Test data transfer from the 1D circle to the 3D torus."""
        # Define parameters for the data transfer process.
        parameters = KM.Parameters("""{
            "origin_variables"         : ["TEMPERATURE"],
            "destination_variables"    : ["TEMPERATURE"]
        }""")

        # Execute the data transfer process.
        process = CoSimulationApplication.DataTransfer3D1DProcess(self.model_part_torus, self.model_part_circle, parameters)
        process.Execute()

        # Define result file path for the output of the transfer.
        result_file = GetFilePath("3d_1d_data_transfer/3d_1d_data_transfer_solutions/torus_from_1d_to_3d")

        # # Generate output results in JSON format.
        # generate_result(self.current_model, result_file, "Torus")

        # Check the results against expected outcomes.
        check_results(self.current_model, result_file, "Torus")

        # Debug output for visualization.
        # debug_vtk(self.current_model, ["Torus", "Circle"])

    def test_torus_from_3d_to_1d(self):
        """Test data transfer from the 3D torus to the 1D circle."""
        # Define parameters for the data transfer process.
        parameters = KM.Parameters("""{
            "origin_variables"         : ["TEMPERATURE"],
            "destination_variables"    : ["TEMPERATURE"]
        }""")

        # Execute the data transfer process with transpose setting.
        process = CoSimulationApplication.DataTransfer3D1DProcess(self.model_part_torus, self.model_part_circle, parameters)
        process.Set(KratosMultiphysics.Mapper.USE_TRANSPOSE)
        process.Execute()

        # Define result file path for the output of the transfer.
        result_file = GetFilePath("3d_1d_data_transfer/3d_1d_data_transfer_solutions/torus_from_3d_to_1d")

        # # Generate output results in JSON format.
        # generate_result(self.current_model, result_file, "Circle")

        # Check the results against expected outcomes.
        check_results(self.current_model, result_file, "Circle")

        # Debug output for visualization.
        # debug_vtk(self.current_model, ["Torus", "Circle"])

def debug_vtk(model, list_names):
    """Generate VTK output for visualization of the model parts."""
    import KratosMultiphysics.vtk_output_process as vtk_output_process

    for name in list_names:
        # Define parameters for VTK output.
        vtk_output_parameters = KratosMultiphysics.Parameters("""{
            "Parameters" : {
                "model_part_name"                    : "",
                "output_precision"                   : 8,
                "output_interval"                    : 2,
                "output_sub_model_parts"             : true,
                "output_path"                        : "",
                "nodal_solution_step_data_variables" : ["TEMPERATURE"],
                "element_flags"                      : []
            }
        }""")
        vtk_output_parameters["Parameters"]["model_part_name"].SetString(name)
        vtk_output_parameters["Parameters"]["output_path"].SetString("vtk_output_" + name)

        # Initialize and execute VTK output process.
        vtk_output_process_torus = vtk_output_process.Factory(vtk_output_parameters, model)
        vtk_output_process_torus.ExecuteInitializeSolutionStep()
        vtk_output_process_torus.ExecuteFinalizeSolutionStep()
        vtk_output_process_torus.PrintOutput()

def check_results(model, input_filename, domain):
    """Check the results against expected output using JSON files."""
    # Define parameters for result checking.
    check_parameters = KM.Parameters("""{
        "check_variables"      : ["TEMPERATURE"],
        "input_file_name"      : "",
        "model_part_name"      : "",
        "historical_value"     : true,
        "time_frequency"       : 0.0
    }""")

    check_parameters["model_part_name"].SetString(domain)
    check_parameters["input_file_name"].SetString(str(input_filename) + "_data_transfer.json")

    # Execute the checking process.
    check = FromJsonCheckResultProcess(model, check_parameters)
    check.ExecuteInitialize()
    check.ExecuteBeforeSolutionLoop()
    check.ExecuteFinalizeSolutionStep()

def generate_result(model, output_filename, domain):
    """Generate output results in JSON format."""
    # Define parameters for output generation.
    out_parameters = KM.Parameters("""{
        "output_variables"     : ["TEMPERATURE"],
        "output_file_name"     : "",
        "model_part_name"      : "",
        "historical_value"     : true,
        "time_frequency"       : 0.0
    }""")

    out_parameters["model_part_name"].SetString(domain)
    out_parameters["output_file_name"].SetString(str(output_filename) + "_data_transfer.json")

    # Execute the output generation process.
    out = JsonOutputProcess(model, out_parameters)
    out.ExecuteInitialize()
    out.ExecuteBeforeSolutionLoop()
    out.ExecuteFinalizeSolutionStep()

# Main execution block to run the unit tests.
if __name__ == '__main__':
    KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.WARNING)  # Set the logging severity level.
    KratosUnittest.main()  # Execute the unit tests.