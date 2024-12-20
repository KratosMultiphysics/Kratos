# Import necessary modules from KratosMultiphysics
import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils
from KratosMultiphysics.testing.utilities import ReadModelPart

# Import the pathlib module
from pathlib import Path

# To detect the os
import platform

# Import math module
import math

def debug_vtk(model_part):
    """
    Generates VTK output for the given model part for debugging purposes.

    Parameters:
    model_part (KratosMultiphysics.ModelPart): The model part to generate VTK output for.

    This function initializes and executes the VTK output process for the provided model part.
    It sets up the necessary parameters, including the model part name and the nodal data value variables.
    The VTK output process is then executed through its initialization, pre-solution loop, and solution step phases.
    """
    # Compute normals
    for cond in model_part.Conditions:
        geom = cond.GetGeometry()
        normal = geom.UnitNormal()
        cond.SetValue(KratosMultiphysics.NORMAL, normal)

    # Generate VTK output
    from KratosMultiphysics.vtk_output_process import VtkOutputProcess
    parameters = KratosMultiphysics.Parameters("""{
        "model_part_name"                : "",
        "nodal_data_value_variables"     : ["NORMAL"],
        "condition_data_value_variables" : ["NORMAL"]
    }
    """)
    parameters["model_part_name"].SetString(model_part.Name)
    output = VtkOutputProcess(model_part.GetModel(), parameters)
    output.ExecuteInitialize()
    output.ExecuteBeforeSolutionLoop()
    output.ExecuteInitializeSolutionStep()
    output.PrintOutput()
    output.ExecuteFinalizeSolutionStep()
    output.ExecuteFinalize()

def calculate_analytical_normal(node):
    """
    Calculates the analytical normal vector for a given node.

    Parameters:
    node (KratosMultiphysics.Node): The node for which to calculate the normal vector.

    Returns:
    KratosMultiphysics.Array3: The normalized vector representing the normal at the node.

    This function computes the normal vector for a node based on its coordinates (X, Y, Z).
    The normal vector is calculated by normalizing the node's position vector.
    """
    norm = math.sqrt(node.X**2 + node.Y**2 + node.Z**2)
    normal = KratosMultiphysics.Array3([node.X / norm, node.Y / norm, node.Z / norm])
    return normal

def calculate_norm(array_3d_value):
    """
    Calculates the Euclidean norm (magnitude) of a 3D vector.

    Parameters:
    array_3d_value (list or tuple of float): A 3D vector represented as a list or tuple of three float values.

    Returns:
    float: The Euclidean norm of the 3D vector.

    This function computes the Euclidean norm of a 3D vector by taking the square root of the sum of the squares
    of its components.
    """
    return math.sqrt(array_3d_value[0]**2 + array_3d_value[1]**2 + array_3d_value[2]**2)

# Define a function to get the file path given a file name
def GetFilePath(fileName):
    """
    Get the absolute file path for a given file name.

    Args:
        fileName (str): The name of the file.

    Returns:
        str: The absolute file path.
    """
    current_dir = Path(__file__).resolve().parent
    return str(current_dir / fileName)

# Define a function to remove files with the given name
def RemoveFiles(file_name):
    """
    Remove files with a given name.

    Args:
        file_name (str): The name of the file to be removed.
    """
    kratos_utils.DeleteFileIfExisting(file_name)

# Define a function to write a Kratos model part to OBJ format
def WriteModelPartToOBJ(model_part, obj_file):
    """
    Write a Kratos model part to an OBJ file.

    Args:
        model_part (KratosMultiphysics.ModelPart): The Kratos model part to be written.
        obj_file (str): The name of the OBJ file to write to.
    """
    if platform.system() == 'Windows':
        import tempfile
        with tempfile.NamedTemporaryFile(delete=False) as temp:
            obj_file = temp.name
    write_settings = KratosMultiphysics.Parameters("""{"open_mode" : "write"}""")
    obj_io = KratosMultiphysics.ObjIO(obj_file, write_settings)
    obj_io.WriteModelPart(model_part)
    return obj_file

# Define a function to read a Kratos model part from an OBJ file
def ReadModelPartFromOBJ(model_part, obj_file, decompose_quad=False):
    """
    Read a Kratos model part from an OBJ file.

    Args:
        model_part (KratosMultiphysics.ModelPart): The Kratos model part to populate with data from the OBJ file.
        data_comm (KratosMultiphysics.DataCommunicator): The data communicator.
        obj_file (str): The name of the OBJ file to read from.
    """
    read_settings = KratosMultiphysics.Parameters("""{
        "open_mode"                      : "read",
        "entity_type"                    : "condition",
        "decompose_quads_into_triangles" : false
    }""")
    read_settings["decompose_quads_into_triangles"].SetBool(decompose_quad)
    obj_io = KratosMultiphysics.ObjIO(obj_file, read_settings)
    obj_io.ReadModelPart(model_part)

class TestObjIO(KratosUnittest.TestCase):
    """
    Define a test class for testing OBJ input/output operations
    """
    @classmethod
    def setUpClass(cls):
        """
        Set up the test class before running tests.
        """
        file = Path.cwd() / "aux.obj"
        cls.obj_file = file.as_posix()

    @classmethod
    def tearDownClass(cls):
        """
        Clean up after running tests.
        """
        RemoveFiles(cls.obj_file)

    def test_ReadObjIO(self):
        """
        Test the ReadModelPart function from ObjIO
        """
        # Create a model part and set the domain size
        self.current_model = KratosMultiphysics.Model()
        self.model_part = self.current_model.CreateModelPart("Main")
        self.model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = 3

        # Read a model part from an OBJ file
        obj_name = GetFilePath("auxiliar_files_for_python_unittest/obj_files/cube.obj")
        ReadModelPartFromOBJ(self.model_part, obj_name)

        # # Debug
        # debug_vtk(self.model_part)

        # Assert that the model part has the correct number of nodes and elements
        self.assertEqual(self.model_part.NumberOfNodes(), 8)
        self.assertEqual(self.model_part.NumberOfConditions(), 6)

        # Assert that the normals are the same as the nodes coordinates
        for node in self.model_part.Nodes:
            normal = node.GetValue(KratosMultiphysics.NORMAL)
            self.assertAlmostEqual(normal[0], node.X)
            self.assertAlmostEqual(normal[1], node.Y)
            self.assertAlmostEqual(normal[2], node.Z)

    def test_ReadObjIOTriangles(self):
        """
        Test the ReadModelPart function from ObjIO
        """
        # Create a model part and set the domain size
        self.current_model = KratosMultiphysics.Model()
        self.model_part = self.current_model.CreateModelPart("Main")
        self.model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = 3

        # Read a model part from an OBJ file
        obj_name = GetFilePath("auxiliar_files_for_python_unittest/obj_files/cube.obj")
        ReadModelPartFromOBJ(self.model_part, obj_name, True)

        # # Debug
        # debug_vtk(self.model_part)

        # Assert that the model part has the correct number of nodes and elements
        self.assertEqual(self.model_part.NumberOfNodes(), 8)
        self.assertEqual(self.model_part.NumberOfConditions(), 12)

        # Assert that the normals are the same as the nodes coordinates
        for node in self.model_part.Nodes:
            normal = node.GetValue(KratosMultiphysics.NORMAL)
            self.assertAlmostEqual(normal[0], node.X)
            self.assertAlmostEqual(normal[1], node.Y)
            self.assertAlmostEqual(normal[2], node.Z)

    def test_ReadObjIOTrianglesDegenerated(self):
        """
        Test the ReadModelPart function from ObjIO
        """
        # Create a model part and set the domain size
        self.current_model = KratosMultiphysics.Model()
        self.model_part = self.current_model.CreateModelPart("Main")
        self.model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = 3

        # Read a model part from an OBJ file
        obj_name = GetFilePath("auxiliar_files_for_python_unittest/obj_files/cube_degenerated.obj")
        ReadModelPartFromOBJ(self.model_part, obj_name, True)

        clean_settings = KratosMultiphysics.Parameters("""{
            "model_part_name"    : "Main",
            "entity_type"        : "condition"
        }""")
        KratosMultiphysics.CleanUpProblematicTrianglesModeler(self.current_model, clean_settings).SetupModelPart()

        # # Debug
        # debug_vtk(self.model_part)

        # Assert that the model part has the correct number of nodes and elements
        self.assertEqual(self.model_part.NumberOfNodes(), 8)
        self.assertEqual(self.model_part.NumberOfConditions(), 12)

    def test_WriteObjIO(self):
        """
        Test the WriteModelPart function form ObjIO
        """
        # Create a model part and set the domain size
        self.current_model = KratosMultiphysics.Model()
        self.model_part = self.current_model.CreateModelPart("Main")
        self.model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = 3

        # Create a sphere model part
        mdpa_name = GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/coarse_sphere_skin")
        ReadModelPart(mdpa_name, self.model_part)

        # Compute normals
        KratosMultiphysics.VariableUtils().SetNonHistoricalVariableToZero(KratosMultiphysics.NORMAL, self.model_part.Nodes)
        KratosMultiphysics.NormalCalculationUtils().CalculateOnSimplexNonHistorical(self.model_part)

        # Compute area of the given mesh
        area_1 = 0.0
        for cond in self.model_part.Conditions:
            geom = cond.GetGeometry()
            area_1 += geom.Area()

        # Write current mesh into OBJ
        self.obj_file = WriteModelPartToOBJ(self.model_part, self.obj_file)

        # Read it back
        obj_model_part = self.current_model.CreateModelPart("OBJ")
        ReadModelPartFromOBJ(obj_model_part, self.obj_file)

        # # Debug
        # debug_vtk(obj_model_part)

        # Compute resulting area
        number_of_conditions = self.model_part.NumberOfConditions()
        area_2 = 0.0
        for elem in obj_model_part.Conditions:
            geom = elem.GetGeometry()
            area_2 += geom.Area()

        # Assert number of nodes and elements
        self.assertEqual(obj_model_part.NumberOfNodes(), self.model_part.NumberOfNodes())
        self.assertEqual(number_of_conditions, obj_model_part.NumberOfConditions())

        # Assert that the areas match approximately
        self.assertAlmostEqual(area_1, area_2)

        # Check normals results
        for node in obj_model_part.Nodes:
            normal = calculate_analytical_normal(node)
            solution_normal = node.GetValue(KratosMultiphysics.NORMAL)
            solution_normal /= calculate_norm(solution_normal)
            self.assertLess(calculate_norm(normal - solution_normal), 0.15)

if __name__ == '__main__':
    # Set the logger severity level and run the KratosUnittest
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()