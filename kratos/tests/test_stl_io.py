import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils
from KratosMultiphysics.testing.utilities import ReadModelPart

# Import os module
import os

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

def RemoveFiles(mdpa_name):
    kratos_utils.DeleteFileIfExisting(mdpa_name + ".time")

class TestStlIO(KratosUnittest.TestCase):
    @classmethod
    def setUpClass(cls):
        KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
        cls.current_model = KratosMultiphysics.Model()
        cls.model_part = cls.current_model.CreateModelPart("Main")
        cls.model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = 3
        cls.stl_file = "aux.stl"

    @classmethod
    def tearDownClass(cls):
        RemoveFiles(cls.stl_file)

    # TODO: Add Read test, this test is added mainly to test Write STL files in MPI

    def test_WriteStlIO(self):
        mdpa_name = GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/coarse_sphere_skin")
        ReadModelPart(mdpa_name, self.model_part)

        # Compute area of given 
        area_1 = 0.0
        for cond in self.model_part.Conditions:
            geom = cond.GetGeometry()
            area_1 += geom.Area()

        # Write current mesh into STL
        write_settings = KratosMultiphysics.Parameters("""{"open_mode" : "write"}""")
        stl_io_1 = KratosMultiphysics.StlIO(self.stl_file, write_settings)
        stl_io_1.WriteModelPart(self.model_part)

        # Read it back
        stl_model_part = self.current_model.CreateModelPart("STL")
        read_settings = KratosMultiphysics.Parameters("""{"open_mode" : "read", "new_entity_type" : "element"}""")
        stl_io_2 = KratosMultiphysics.StlIO(self.stl_file, read_settings)
        stl_io_2.ReadModelPart(stl_model_part)

        # Compute resulting area
        area_2 = 0.0
        for elem in stl_model_part.Elements:
            geom = elem.GetGeometry()
            area_2 += geom.Area()

        self.assertAlmostEqual(area_1, area_2)

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()