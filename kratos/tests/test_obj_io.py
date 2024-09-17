# Import necessary modules from KratosMultiphysics
import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils
from KratosMultiphysics.testing.utilities import ReadModelPart

# Import the pathlib module
from pathlib import Path

# To detect the os
import platform

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
        KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
        cls.current_model = KratosMultiphysics.Model()
        cls.model_part = cls.current_model.CreateModelPart("Main")
        cls.model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = 3
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
        # Read a model part from an OBJ file
        obj_name = GetFilePath("auxiliar_files_for_python_unittest/obj_files/cube.obj")
        ReadModelPartFromOBJ(self.model_part, obj_name)

        # Assert that the model part has the correct number of nodes and elements
        self.assertEqual(self.model_part.NumberOfNodes(), 8)
        self.assertEqual(self.model_part.NumberOfElements(), 6)

        # Assert that the normals are the same as the nodes coordinates
        for node in self.model_part.Nodes:
            normal = node.GetValue(KratosMultiphysics.NORMAL)
            self.assertAlmostEqual(normal[0], node.X)
            self.assertAlmostEqual(normal[1], node.Y)
            self.assertAlmostEqual(normal[2], node.Z)

if __name__ == '__main__':
    # Set the logger severity level and run the KratosUnittest
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()