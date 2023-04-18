# Importing the Kratos Library
import KratosMultiphysics as KM

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Importing the Kratos utilities
import KratosMultiphysics.kratos_utilities as kratos_utils

# Importing the read model part utility from the testing utilities
from KratosMultiphysics.testing.utilities import ReadModelPart

# Importing the OS module
import os

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

def RemoveFiles(mdpa_name):
    kratos_utils.DeleteFileIfExisting(mdpa_name + ".time")

class TestGeometricalObjectBins(KratosUnittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.current_model = KM.Model()
        cls.model_part = cls.current_model.CreateModelPart("Main")
        cls.model_part.ProcessInfo[KM.DOMAIN_SIZE] = 3
        cls.model_part.AddNodalSolutionStepVariable(KM.BULK_MODULUS)
        cls.model_part.AddNodalSolutionStepVariable(KM.NODAL_VAUX)
        cls.model_part.AddNodalSolutionStepVariable(KM.EXTERNAL_FORCES_VECTOR)
        cls.model_part.AddNodalSolutionStepVariable(KM.LOCAL_AXES_MATRIX)
        cls.mdpa_name = GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/coarse_sphere_with_conditions")
        ReadModelPart(cls.mdpa_name, cls.model_part)

    @classmethod
    def tearDownClass(cls):
        RemoveFiles(cls.mdpa_name)

    def setUp(self):
        # Create search
        self.search = KM.GeometricalObjectsBins(self.model_part.Conditions)

        # Create node for search
        self.node = self.model_part.CreateNewNode(100000, 0.0, 0.0, 0.15)

    def test_GeometricalObjectsBins_SearchNearestInRadius(self):
        radius = 0.35
        results = self.search.SearchInRadius(self.node, radius)

        self.assertEqual(len(results), 8)
        cond_id_ref = [125,78,117,18,68,1,41,119]
        for result in results:
            self.assertTrue(result.Get().Id in cond_id_ref)

if __name__ == '__main__':
    KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.WARNING)
    KratosUnittest.main()