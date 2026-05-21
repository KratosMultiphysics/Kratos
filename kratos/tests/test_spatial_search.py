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

class TestSpatialSearchSphere(KratosUnittest.TestCase):
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
        self.settings = KM.Parameters("""
        {
            "container_type"  : "",
            "search_parameters" : {
                "bucket_size"     : 4,
                "allocation_size" : 1000
            }
        }
        """)

    def test_KDTree_nodes(self):
        # Create search
        self.settings["container_type"].SetString("KDTree")
        self.search = KM.SpecializedSpatialSearch(self.settings)

        # Create node for search
        second_model_part = self.current_model.CreateModelPart("KDTree")
        second_model_part.CreateNewNode(100000, 0.0, 0.0, 0.0)
        radius_list = [0.3]
        [results, distances] = self.search.SearchNodesInRadiusExclusive(self.model_part, second_model_part.Nodes, radius_list)

        self.assertEqual(len(results), 1)
        self.assertEqual(len(results[0]), 7)
        distance_ref = [0.077385615, 0.0008331217999999999, 0.0899807529, 0.0627019979, 0.07703137859999999, 0.0789991779, 0.0708403121]
        node_id_ref = [7, 17, 18, 23, 33, 39, 44, 56]
        for distance in distances[0]:
            self.assertTrue(distance in distance_ref)
        for node in results[0]:
            self.assertTrue(node.Id in node_id_ref)

    def test_Octree_nodes(self):
        # Create search
        self.settings["container_type"].SetString("Octree")
        self.search = KM.SpecializedSpatialSearch(self.settings)

        # Create node for search
        second_model_part = self.current_model.CreateModelPart("Octree")
        second_model_part.CreateNewNode(100000, 0.0, 0.0, 0.0)
        radius_list = [0.3]
        [results, distances] = self.search.SearchNodesInRadiusExclusive(self.model_part, second_model_part.Nodes, radius_list)

        self.assertEqual(len(results), 1)
        self.assertEqual(len(results[0]), 7)
        distance_ref = [0.077385615, 0.0008331217999999999, 0.0899807529, 0.0627019979, 0.07703137859999999, 0.0789991779, 0.0708403121]
        node_id_ref = [7, 17, 18, 23, 33, 39, 44, 56]
        for distance in distances[0]:
            self.assertTrue(distance in distance_ref)
        for node in results[0]:
            self.assertTrue(node.Id in node_id_ref)

    def test_BinsStatic_nodes(self):
        # Create search
        self.settings["container_type"].SetString("BinsStatic")
        self.search = KM.SpecializedSpatialSearch(self.settings)

        # Create node for search
        second_model_part = self.current_model.CreateModelPart("BinsStatic")
        second_model_part.CreateNewNode(100000, 0.0, 0.0, 0.0)
        radius_list = [0.3]
        [results, distances] = self.search.SearchNodesInRadiusExclusive(self.model_part, second_model_part.Nodes, radius_list)

        self.assertEqual(len(results), 1)
        self.assertEqual(len(results[0]), 7)
        distance_ref = [0.077385615, 0.0008331217999999999, 0.0899807529, 0.0627019979, 0.07703137859999999, 0.0789991779, 0.0708403121]
        node_id_ref = [7, 17, 18, 23, 33, 39, 44, 56]
        for distance in distances[0]:
            self.assertTrue(distance in distance_ref)
        for node in results[0]:
            self.assertTrue(node.Id in node_id_ref)

    def test_BinsDynamic_nodes(self):
        # Create search
        self.settings["container_type"].SetString("BinsDynamic")
        self.search = KM.SpecializedSpatialSearch(self.settings)

        # Create node for search
        second_model_part = self.current_model.CreateModelPart("BinsDynamic")
        second_model_part.CreateNewNode(100000, 0.0, 0.0, 0.0)
        radius_list = [0.3]
        [results, distances] = self.search.SearchNodesInRadiusExclusive(self.model_part, second_model_part.Nodes, radius_list)

        self.assertEqual(len(results), 1)
        self.assertEqual(len(results[0]), 7)
        distance_ref = [0.077385615, 0.0008331217999999999, 0.0899807529, 0.0627019979, 0.07703137859999999, 0.0789991779, 0.0708403121]
        node_id_ref = [7, 17, 18, 23, 33, 39, 44, 56]
        for distance in distances[0]:
            self.assertTrue(distance in distance_ref)
        for node in results[0]:
            self.assertTrue(node.Id in node_id_ref)

if __name__ == '__main__':
    KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.WARNING)
    KratosUnittest.main()