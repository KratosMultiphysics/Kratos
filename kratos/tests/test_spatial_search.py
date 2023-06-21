# Importing the Kratos Library
import KratosMultiphysics as KM

# Importing the MPI
if KM.IsDistributedRun():
    import KratosMultiphysics.mpi as KratosMPI

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
        # Adding PARTITION_INDEX
        if KM.IsDistributedRun():
            cls.model_part.AddNodalSolutionStepVariable(KM.PARTITION_INDEX)
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
        
        # Get the data communicator
        self.data_comm = self.model_part.GetCommunicator().GetDataCommunicator()

    def GenerateSearch(self, container_type = "KDTree"):
        self.settings["container_type"].SetString(container_type)
        if KM.IsDistributedRun():
            raise Exception("MPI version comming in a future PR")
        else:
            self.search = KM.SpecializedSpatialSearch(self.settings)

        # Creating submodel part
        self.second_model_part = self.current_model.CreateModelPart(container_type)
        # Adding PARTITION_INDEX
        if KM.IsDistributedRun():
            self.second_model_part.AddNodalSolutionStepVariable(KM.PARTITION_INDEX)

        # Create node for search
        self.new_node_id = 100000
        self.point = KM.Point(0.0, 0.0, 0.0)
        if KM.IsDistributedRun():
            # Only added to first rank to actualy check it works in all ranks
            if self.data_comm.Rank() == 0:
                self.node = self.second_model_part.CreateNewNode(self.new_node_id, self.point.X, self.point.Y, self.point.Z)
                self.node.SetSolutionStepValue(KM.PARTITION_INDEX, 0)
            #KratosMPI.ParallelFillCommunicator(self.model_part, self.data_comm).Execute()
            KratosMPI.ParallelFillCommunicator(self.second_model_part, self.data_comm).Execute()
        else:
            self.node = self.second_model_part.CreateNewNode(self.new_node_id, self.point.X, self.point.Y, self.point.Z)

        # Reference solution
        radius_list = [0.3]
        distance_ref = [0.077385615, 0.0008331217999999999, 0.0899807529, 0.0627019979, 0.07703137859999999, 0.0789991779, 0.0708403121]
        node_id_ref = [7, 17, 18, 23, 33, 39, 44, 56]

        # Serial interface
        if not KM.IsDistributedRun():
            [results, distances] = self.search.SearchNodesInRadiusExclusive(self.model_part, self.second_model_part.Nodes, radius_list)

            # Assert results
            self.assertEqual(len(results), len(radius_list))
            self.assertEqual(len(results[0]), len(distance_ref))
            for distance in distances[0]:
                self.assertTrue(distance in distance_ref)
            for node in results[0]:
                self.assertTrue(node.Id in node_id_ref)
        
        # Parallel interface (also works in serial mode)
        results = self.search.SearchNodesInRadiusExclusive(self.model_part, self.second_model_part.Nodes, radius_list, self.data_comm)
        self.assertEqual(results.NumberOfSearchResults(), len(radius_list))
        self.assertTrue(results.HasResult(self.new_node_id))
        node_results = results[self.new_node_id]
        if not KM.IsDistributedRun():
            self.assertEqual(node_results.NumberOfLocalResults(), len(distance_ref))
        self.assertEqual(node_results.NumberOfGlobalResults(), len(distance_ref))
        ids = node_results.GetResultIndices()
        for id in ids:
            self.assertTrue(id in node_id_ref)
        distances = node_results.GetDistances()
        for distance in distances:
            self.assertTrue(distance in distance_ref)

        # New interface (works in serial and parallel mode). Interface is very similar to previous one, but it is because a node is type of point and this interface is purely for points (nodes in python)
        results = self.search.SearchNodesOverPointsInRadius(self.model_part.Nodes, self.second_model_part.Nodes, radius_list, self.data_comm)
        self.assertEqual(results.NumberOfSearchResults(), len(radius_list))
        self.assertFalse(results.HasResult(self.new_node_id)) # This is false because the hash is with coordinates
        self.assertTrue(results.HasResult(self.point))
        node_results = results[self.point]
        if not KM.IsDistributedRun():
            self.assertEqual(node_results.NumberOfLocalResults(), len(distance_ref))
        self.assertEqual(node_results.NumberOfGlobalResults(), len(distance_ref))
        ids = node_results.GetResultIndices()
        for id in ids:
            self.assertTrue(id in node_id_ref)
        distances = node_results.GetDistances()
        for distance in distances:
            self.assertTrue(distance in distance_ref)

    def test_KDTree_nodes(self):
        # Create search
        self.GenerateSearch("KDTree")

    def test_Octree_nodes(self):
        # Create search
        self.GenerateSearch("Octree")

    def test_BinsStatic_nodes(self):
        # Create search
        self.GenerateSearch("BinsStatic")

    def test_BinsDynamic_nodes(self):
        # Create search
        self.GenerateSearch("BinsDynamic")

if __name__ == '__main__':
    KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.WARNING)
    KratosUnittest.main()