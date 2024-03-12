# Import necessary modules from KratosMultiphysics and other standard libraries
import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils
from KratosMultiphysics.testing.utilities import ReadModelPart

# Import os module
import os

# Import math module
import math

def GetFilePath(fileName):
    """Construct the full file path for a given file name, located in the same directory as this script.

    Args:
        fileName (str): Name of the file to construct the path for.

    Returns:
        str: The full file path.
    """
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

def RemoveFiles(mdpa_name):
    """Remove files related to a specific model part, typically used for cleanup.

    Args:
        mdpa_name (str): Base name of the file (without extension) to be removed.
    """
    kratos_utils.DeleteFileIfExisting(mdpa_name + ".time")

def post_process_vtk(model):
    """Generate VTK files for visualization of results using the VTK output process.

    Args:
        model (KratosMultiphysics.Model): The model containing the model part to be visualized.
    """
    # Import the VTK output process and define the output parameters
    import KratosMultiphysics.vtk_output_process as vtk_output_process
    vtk_output_parameters = KratosMultiphysics.Parameters("""
    {
        "Parameters" : {
            "model_part_name"                    : "main_model_part",
            "file_format"                        : "ascii",
            "output_precision"                   : 8,
            "output_interval"                    : 2,
            "output_sub_model_parts"             : true,
            "entity_type"                        : "element",
            "nodal_data_value_variables"         : ["DISTANCE"],
            "output_path"                        : "test_vtk_output"
        }
    }""")

    # Initialize and execute the VTK output process
    process = vtk_output_process.Factory(vtk_output_parameters, model)
    process.ExecuteInitialize()
    process.ExecuteInitializeSolutionStep()
    process.ExecuteFinalizeSolutionStep()
    process.PrintOutput()

def post_process_gid(model_part):
    """Generate files for visualization using GiD by executing the GiD output process.

    Args:
        model_part (KratosMultiphysics.ModelPart): The model part to be visualized in GiD.
    """
    # Import the GiD output process and define the output configuration
    from KratosMultiphysics.gid_output_process import GiDOutputProcess
    gid_output = GiDOutputProcess(model_part,
                                "gid_output",
                                KratosMultiphysics.Parameters("""
                                    {
                                        "result_file_configuration" : {
                                            "gidpost_flags": {
                                                "GiDPostMode": "GiD_PostBinary",
                                                "WriteDeformedMeshFlag": "WriteUndeformed",
                                                "WriteConditionsFlag": "WriteConditions",
                                                "MultiFileFlag": "SingleFile"
                                            },
                                            "nodal_nonhistorical_results" : ["DISTANCE"]
                                        }
                                    }
                                    """)
                                )

    # Execute the GiD output process
    gid_output.ExecuteInitialize()
    gid_output.ExecuteBeforeSolutionLoop()
    gid_output.ExecuteInitializeSolutionStep()
    gid_output.PrintOutput()
    gid_output.ExecuteFinalizeSolutionStep()
    gid_output.ExecuteFinalize()

def CalculateAnalyticalDistance(node):
    """Calculate the analytical distance of a node from the origin.

    Args:
        node (KratosMultiphysics.Node): The node to calculate the distance for.

    Returns:
        float: The calculated distance.
    """
    distance = 0.5 - math.sqrt(node.X**2 + node.Y**2 + node.Z**2)
    return distance

class TestCalculateNodalDistanceToSkinProcessCoarseSphere(KratosUnittest.TestCase):
    """Test case for verifying the nodal distance calculation to a skin surface on a coarse sphere model."""

    @classmethod
    def setUpClass(cls):
        """Prepare the test environment with model parts and their respective .mdpa files."""
        cls.current_model = KratosMultiphysics.Model()
        cls.model_part = cls.current_model.CreateModelPart("main_model_part")
        cls.model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = 3
        cls.mdpa_name = GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/coarse_sphere")
        ReadModelPart(cls.mdpa_name, cls.model_part)
        cls.skin_model_part = cls.current_model.CreateModelPart("skin_model_part")
        cls.skin_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = 3
        cls.skin_mdpa_name = GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/coarse_sphere_skin")
        ReadModelPart(cls.skin_mdpa_name, cls.skin_model_part)

    @classmethod
    def tearDownClass(cls):
        """Clean up by removing the generated files."""
        RemoveFiles(cls.mdpa_name)
        RemoveFiles(cls.skin_mdpa_name)

    def setUp(self):
        """Setup tasks before each test method."""
        pass

    @KratosUnittest.skipIf(KratosMultiphysics.IsDistributedRun(), "This test is designed for serial runs only.")
    def test_ComputeDistanceToSkin(self):
        """Test the computation of nodal distances to the skin in the model part."""
        # Define the settings for the distance calculation process
        settings = KratosMultiphysics.Parameters("""
        {
            "distance_database"      : "nodal_non_historical",
            "distance_variable"      : "DISTANCE",
            "volume_model_part"      : "main_model_part",
            "skin_model_part"        : "skin_model_part"
        }""")
        # Execute the distance calculation process
        KratosMultiphysics.CalculateNodalDistanceToSkinProcess(self.current_model, settings).Execute()

        ## DEBUG
        # Uncomment to visualize output
        #post_process_gid(self.model_part)
        #post_process_vtk(self.current_model)

        ## Check the calculated distances against analytical values
        for node in self.model_part.Nodes:
            distance = CalculateAnalyticalDistance(node)
            solution_distance = node.GetValue(KratosMultiphysics.DISTANCE)
            self.assertAlmostEqual(distance, solution_distance, delta=3.8e-2)

if __name__ == '__main__':
    # Configure logging level and start the test runner
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()
