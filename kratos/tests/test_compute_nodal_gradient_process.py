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
            "nodal_data_value_variables"         : ["DISTANCE", "DISTANCE_GRADIENT"],
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
                                            "nodal_nonhistorical_results" : ["DISTANCE", "DISTANCE_GRADIENT"]
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

def CalculateAnalyticalNormal(node):
    """
    Calculates the normalized vector representing the normal direction at a given node.

    Parameters:
        node : Node object
            An object representing a node in the computational domain.

    Returns:
        Array3
            A KratosMultiphysics Array3 object representing the normalized normal vector at the given node.

    Note:
        This function assumes that the input node has attributes X, Y, and Z representing its coordinates.
    """
    norm = math.sqrt(node.X**2+node.Y**2+node.Z**2)
    normal = KratosMultiphysics.Array3([- node.X/norm, - node.Y/norm, - node.Z/norm])
    return normal

def CalculateAnalyticalDistanceGradientNorm(node):
    """
    Calculates the norm of the analytical distance gradient vector at a given node.

    Parameters:
        node : Node object
            An object representing a node in the computational domain.

    Returns:
        float
            The norm of the analytical distance gradient vector at the given node.

    Note:
        This function assumes that the input node has attributes X, Y, and Z representing its coordinates.
        The distance gradient is calculated assuming a constant gradient of 0.5.
    """
    return math.sqrt(node.X**2 + node.Y**2 + node.Z**2)/0.5

def CalculateNorm(array_3d_value):
    """
    Calculates the Euclidean norm of a 3D vector.

    Parameters:
        array_3d_value : list or tuple
            A list or tuple representing a 3D vector in the form [x, y, z].

    Returns:
        float
            The Euclidean norm of the input 3D vector.

    Note:
        This function assumes that the input vector has three elements corresponding to its x, y, and z components.
    """
    return math.sqrt(array_3d_value[0]**2+array_3d_value[1]**2+array_3d_value[2]**2)

class TestComputeNodalGradientProcessCoarseSphere(KratosUnittest.TestCase):
    """Test case for verifying the nodal distance gradient calculation on a coarse sphere model."""

    @classmethod
    def setUpClass(cls):
        """Prepare the test environment with model parts and their respective .mdpa files."""
        cls.current_model = KratosMultiphysics.Model()
        cls.model_part = cls.current_model.CreateModelPart("main_model_part")
        cls.model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = 3
        cls.mdpa_name = GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/coarse_sphere")
        ReadModelPart(cls.mdpa_name, cls.model_part)

    @classmethod
    def tearDownClass(cls):
        """Clean up by removing the generated files."""
        RemoveFiles(cls.mdpa_name)

    def setUp(self):
        """Setup tasks before each test method."""
        pass

    def test_ComputeNodalGradient(self):
        """Test the computation of nodal gradient in the model part."""
        # Set the values
        for node in self.model_part.Nodes:
            distance = CalculateAnalyticalDistance(node)
            node.SetValue(KratosMultiphysics.DISTANCE, distance)

        # Define the settings for the distance gradient calculation process
        settings = KratosMultiphysics.Parameters("""
        {
            "model_part_name"                : "main_model_part",
            "origin_variable"                : "DISTANCE",
            "gradient_variable"              : "DISTANCE_GRADIENT",
            "area_variable"                  : "NODAL_AREA",
            "non_historical_origin_variable" :  true
        }""")
        KratosMultiphysics.ComputeNonHistoricalNodalGradientProcess(self.current_model, settings).Execute()

        ## DEBUG
        # Uncomment to visualize output
        #post_process_gid(self.model_part)
        #post_process_vtk(self.current_model)

        ## Check the calculated distances gradients norms against analytical values
        nodes_to_test = [2, 4, 11, 26, 29, 30, 31, 35, 40, 43, 45, 50, 63, 64, 66, 68, 69, 70, 74, 76, 77, 78, 79, 81, 83]
        for node in self.model_part.Nodes:
            if node.Id in nodes_to_test:
                #distance = CalculateAnalyticalDistance(node)
                solution_gradient_distance = node.GetValue(KratosMultiphysics.DISTANCE_GRADIENT)
                solution_gradient_distance_norm = CalculateNorm(solution_gradient_distance)
                analytical_solution_gradient_distance_norm = CalculateAnalyticalDistanceGradientNorm(node)
                solution_gradient_distance_direction = solution_gradient_distance/solution_gradient_distance_norm
                normal = CalculateAnalyticalNormal(node)
                #print("Node ", node.Id, "\t", solution_gradient_distance_norm, "\t", analytical_solution_gradient_distance_norm)
                self.assertAlmostEqual(solution_gradient_distance_norm, analytical_solution_gradient_distance_norm, delta=5.0e-2)
                #print("Node ", node.Id, "\t", solution_gradient_distance_direction, "\t", normal, "\t", CalculateNorm(normal - solution_gradient_distance_direction))
                self.assertLess(CalculateNorm(normal - solution_gradient_distance_direction), 0.15)

if __name__ == '__main__':
    # Configure logging level and start the test runner
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()