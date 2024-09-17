# Import necessary modules from KratosMultiphysics
import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils
from KratosMultiphysics.testing.utilities import ReadModelPart

# Import the pathlib module
from pathlib import Path

# To detect the os
import platform

def debug_vtk(model_part):
    from KratosMultiphysics.vtk_output_process import VtkOutputProcess
    parameters = KratosMultiphysics.Parameters("""{
        "model_part_name"            : "",
        "nodal_data_value_variables" : ["NORMAL"]
    }
    """)
    parameters["model_part_name"].SetString(model_part.Name)
    output = VtkOutputProcess(model_part.GetModel(), parameters )
    output.ExecuteInitialize()
    output.ExecuteBeforeSolutionLoop()
    output.ExecuteInitializeSolutionStep()

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
def ReadModelPartFromOBJ(model_part, obj_file):
    """
    Read a Kratos model part from an OBJ file.

    Args:
        model_part (KratosMultiphysics.ModelPart): The Kratos model part to populate with data from the OBJ file.
        data_comm (KratosMultiphysics.DataCommunicator): The data communicator.
        obj_file (str): The name of the OBJ file to read from.
    """
    read_settings = KratosMultiphysics.Parameters("""{"open_mode" : "read", "entity_type" : "element"}""")
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
        self.assertEqual(self.model_part.NumberOfElements(), 6)

        # Assert that the normals are the same as the nodes coordinates
        for node in self.model_part.Nodes:
            normal = node.GetValue(KratosMultiphysics.NORMAL)
            self.assertAlmostEqual(normal[0], node.X)
            self.assertAlmostEqual(normal[1], node.Y)
            self.assertAlmostEqual(normal[2], node.Z)

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
        # debug_vtk(self.model_part)

        # Compute resulting area
        number_of_conditions = self.model_part.NumberOfConditions()
        area_2 = 0.0
        for elem in obj_model_part.Elements:
            geom = elem.GetGeometry()
            area_2 += geom.Area()

        # Assert number of nodes and elements
        self.assertEqual(obj_model_part.NumberOfNodes(), 384)
        self.assertEqual(number_of_conditions, obj_model_part.NumberOfElements())

        # Assert that the areas match approximately
        self.assertAlmostEqual(area_1, area_2)

if __name__ == '__main__':
    # Set the logger severity level and run the KratosUnittest
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()