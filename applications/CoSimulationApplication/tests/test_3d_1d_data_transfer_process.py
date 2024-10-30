# Import Kratos
import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the mapping application
import KratosMultiphysics.MappingApplication

# Import the extension application
import KratosMultiphysics.CoSimulationApplication  as CoSimulationApplication

# Import the test utilities
from KratosMultiphysics.from_json_check_result_process import FromJsonCheckResultProcess
from KratosMultiphysics.json_output_process import JsonOutputProcess

# Import basic dependencies
import math
import os

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class Test3D1DDataTransferProcessBlock(KratosUnittest.TestCase):
    def setUp(self):
        self.current_model = KM.Model()

        # Import Block
        self.model_part_block = self.current_model.CreateModelPart("Block")
        self.model_part_block.ProcessInfo[KM.STEP] = 0
        self.model_part_block.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)
        self.model_part_block.ProcessInfo.SetValue(KM.TIME, 0.0)
        self.model_part_block.ProcessInfo.SetValue(KM.DELTA_TIME, 1.0)
        self.model_part_block.AddNodalSolutionStepVariable(KM.TEMPERATURE)
        input_mdpa = GetFilePath("3d_1d_data_transfer/3d_1d_data_transfer_mdpa/block")
        model_part_io_block = KM.ModelPartIO(input_mdpa)
        model_part_io_block.ReadModelPart(self.model_part_block)

        # Import Line
        self.model_part_line = self.current_model.CreateModelPart("Line")
        self.model_part_line.ProcessInfo[KM.STEP] = 0
        self.model_part_line.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)
        self.model_part_line.ProcessInfo.SetValue(KM.TIME, 0.0)
        self.model_part_line.ProcessInfo.SetValue(KM.DELTA_TIME, 1.0)
        self.model_part_line.AddNodalSolutionStepVariable(KM.TEMPERATURE)
        number_of_line_elements = 100
        initial_z_coordinate = 0.0
        end_z_coordinate = 5.0
        self.model_part_line.CreateNewNode(1, 0.5, 0.5, initial_z_coordinate)
        prop = self.model_part_line.GetProperties()[1]
        for i in range(number_of_line_elements):
            z = initial_z_coordinate + (i+1) * (end_z_coordinate - initial_z_coordinate) / number_of_line_elements
            self.model_part_line.CreateNewNode(i + 2, 0.5, 0.5, z)
            self.model_part_line.CreateNewElement("Element3D2N", i + 1, [i + 1, i + 2], prop)

        # Set values
        for node in self.model_part_block.Nodes:
            node.SetSolutionStepValue(KM.TEMPERATURE,   node.Z)
        del node

        for node in self.model_part_line.Nodes:
            node.SetSolutionStepValue(KM.TEMPERATURE, - node.Z)
        del node

    def tearDown(self):
        pass

    @KratosUnittest.skipIf(KM.IsDistributedRun(), "This test is designed for serial runs only.")
    def test_block_from_1d_to_3d(self):

        # Data transfer
        parameters = KM.Parameters("""{
            "origin_variables"         : ["TEMPERATURE"],
            "destination_variables"    : ["TEMPERATURE"]
        }""")
        process = CoSimulationApplication.DataTransfer3D1DProcess(self.model_part_block, self.model_part_line, parameters)
        process.Execute()

        result_file = GetFilePath("3d_1d_data_transfer/3d_1d_data_transfer_solutions/block_from_1d_to_3d")

        # Check results
        # generate_result(self.current_model, result_file, "Block")
        check_results(self.current_model, result_file, "Block")

        # # Debug
        debug_vtk(self.current_model, ["Block", "Line"])

    def test_block_from_3d_to_1d(self):

        # Data transfer
        parameters = KM.Parameters("""{
            "origin_variables"         : ["TEMPERATURE"],
            "destination_variables"    : ["TEMPERATURE"]
        }""")
        process = CoSimulationApplication.DataTransfer3D1DProcess(self.model_part_block, self.model_part_line, parameters)
        process.Set(KratosMultiphysics.Mapper.USE_TRANSPOSE)
        process.Execute()

        result_file = GetFilePath("3d_1d_data_transfer/3d_1d_data_transfer_solutions/block_from_3d_to_1d")

        # Check results
        # generate_result(self.current_model, result_file, "Line")
        check_results(self.current_model, result_file, "Line")

        # # Debug
        # debug_vtk(self.current_model, ["Block", "Line"])

class Test3D1DDataTransferProcessTorus(KratosUnittest.TestCase):
    def setUp(self):
        self.current_model = KM.Model()

        # Import Torus
        self.model_part_torus = self.current_model.CreateModelPart("Torus")
        self.model_part_torus.ProcessInfo[KM.STEP] = 0
        self.model_part_torus.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)
        self.model_part_torus.ProcessInfo.SetValue(KM.TIME, 0.0)
        self.model_part_torus.ProcessInfo.SetValue(KM.DELTA_TIME, 1.0)
        self.model_part_torus.AddNodalSolutionStepVariable(KM.TEMPERATURE)
        input_mdpa = GetFilePath("3d_1d_data_transfer/3d_1d_data_transfer_mdpa/torus3d")
        model_part_io_torus = KM.ModelPartIO(input_mdpa)
        model_part_io_torus.ReadModelPart(self.model_part_torus)

        # Import circle
        self.model_part_circle = self.current_model.CreateModelPart("Circle")
        self.model_part_circle.ProcessInfo[KM.STEP] = 0
        self.model_part_circle.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)
        self.model_part_circle.ProcessInfo.SetValue(KM.TIME, 0.0)
        self.model_part_circle.ProcessInfo.SetValue(KM.DELTA_TIME, 1.0)
        self.model_part_circle.AddNodalSolutionStepVariable(KM.TEMPERATURE)
        number_of_circle_elements = 100
        radius = 1.0
        self.model_part_circle.CreateNewNode(1, radius, 0.0, 0)
        prop = self.model_part_circle.GetProperties()[1]
        for i in range(number_of_circle_elements - 1):
            angle = 2 * math.pi * (i+1) / number_of_circle_elements
            self.model_part_circle.CreateNewNode(i + 2, radius*math.cos(angle), radius*math.sin(angle), 0)
            self.model_part_circle.CreateNewElement("Element3D2N", i + 1, [i + 1, i + 2], prop)
        self.model_part_circle.CreateNewElement("Element3D2N", number_of_circle_elements, [i + 2, 1], prop)

        # Set values
        for node in self.model_part_torus.Nodes:
            node.SetSolutionStepValue(KM.TEMPERATURE, - node.X + node.Y)
        del node

        for node in self.model_part_circle.Nodes:
            node.SetSolutionStepValue(KM.TEMPERATURE,   node.X - node.Y)
        del node

    def tearDown(self):
        pass

    def test_torus_from_1d_to_3d(self):

        # Data transfer
        parameters = KM.Parameters("""{
            "origin_variables"         : ["TEMPERATURE"],
            "destination_variables"    : ["TEMPERATURE"]
        }""")
        process = CoSimulationApplication.DataTransfer3D1DProcess(self.model_part_torus, self.model_part_circle, parameters)
        process.Execute()

        result_file = GetFilePath("3d_1d_data_transfer/3d_1d_data_transfer_solutions/torus_from_1d_to_3d")

        # Check results
        # generate_result(self.current_model, result_file, "Torus")
        check_results(self.current_model, result_file, "Torus")

        # # Debug
        debug_vtk(self.current_model, ["Torus", "Circle"])

    def test_torus_from_3d_to_1d(self):

        # Data transfer
        parameters = KM.Parameters("""{
            "origin_variables"         : ["TEMPERATURE"],
            "destination_variables"    : ["TEMPERATURE"]
        }""")
        process = CoSimulationApplication.DataTransfer3D1DProcess(self.model_part_torus, self.model_part_circle, parameters)
        process.Set(KratosMultiphysics.Mapper.USE_TRANSPOSE)
        process.Execute()

        result_file = GetFilePath("3d_1d_data_transfer/3d_1d_data_transfer_solutions/torus_from_3d_to_1d")

        # Check results
        # generate_result(self.current_model, result_file, "Circle")
        check_results(self.current_model, result_file, "Circle")

        # # Debug
        # debug_vtk(self.current_model, ["Torus", "Circle"])

def debug_vtk(model, list_names):
    import KratosMultiphysics.vtk_output_process as vtk_output_process

    for name in list_names:
        vtk_output_parameters = KratosMultiphysics.Parameters("""{
            "Parameters" : {
                "model_part_name"                    : "",
                "file_format"                        : "ascii",
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

        vtk_output_process_torus = vtk_output_process.Factory(vtk_output_parameters, model)
        vtk_output_process_torus.ExecuteInitializeSolutionStep()
        vtk_output_process_torus.ExecuteFinalizeSolutionStep()
        vtk_output_process_torus.PrintOutput()

def check_results(model, input_filename, domain):
    check_parameters = KM.Parameters("""
    {
        "check_variables"      : ["TEMPERATURE"],
        "input_file_name"      : "",
        "model_part_name"      : "",
        "historical_value"     : true,
        "time_frequency"       : 0.0
    }
    """)

    check_parameters["model_part_name"].SetString(domain)
    check_parameters["input_file_name"].SetString(input_filename + "_data_transfer.json")

    check = FromJsonCheckResultProcess(model, check_parameters)
    check.ExecuteInitialize()
    check.ExecuteBeforeSolutionLoop()
    check.ExecuteFinalizeSolutionStep()

def generate_result(model, output_filename, domain):
    out_parameters = KM.Parameters("""
    {
        "output_variables"     : ["TEMPERATURE"],
        "output_file_name"     : "",
        "model_part_name"      : "",
        "historical_value"     : true,
        "time_frequency"       : 0.0
    }
    """)

    out_parameters["model_part_name"].SetString(domain)
    out_parameters["output_file_name"].SetString(output_filename + "_data_transfer.json")

    out = JsonOutputProcess(model, out_parameters)
    out.ExecuteInitialize()
    out.ExecuteBeforeSolutionLoop()
    out.ExecuteFinalizeSolutionStep()

if __name__ == '__main__':
    KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.WARNING)
    KratosUnittest.main()
