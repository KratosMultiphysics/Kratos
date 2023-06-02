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

class TestGeometricalObjectBins(KratosUnittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.current_model = KM.Model()
        cls.model_part = cls.current_model.CreateModelPart("Main")
        cls.sub_model_part = cls.model_part.CreateSubModelPart("SubModelPart")
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
        # Create search
        self.data_comm = self.model_part.GetCommunicator().GetDataCommunicator()
        if KM.IsDistributedRun():
            self.search = KratosMPI.GeometricalObjectsBinsMPI(self.model_part.Conditions, self.data_comm)
        else:
            self.search = KM.GeometricalObjectsBins(self.model_part.Conditions)

        # Create node for search
        self.node_coordinates = KM.Point(0.0, 0.0, 0.15)
        if KM.IsDistributedRun():
            # Only added to first rank to actualy check it works in all ranks
            if self.data_comm.Rank() == 0:
                self.node = self.sub_model_part.CreateNewNode(100000, self.node_coordinates.X, self.node_coordinates.Y, self.node_coordinates.Z)
                self.node.SetSolutionStepValue(KM.PARTITION_INDEX, 0)
            ParallelFillCommunicator = KratosMPI.ParallelFillCommunicator(self.model_part)
            ParallelFillCommunicator.Execute()
        else:
            self.node = self.sub_model_part.CreateNewNode(100000, self.node_coordinates.X, self.node_coordinates.Y, self.node_coordinates.Z)

    def test_GeometricalObjectsBins_SearchInRadius(self):
        # Define radius
        radius = 0.35

        # Reference solution
        cond_id_ref = [125,78,117,18,68,1,41,119]

        # One node search
        results = self.search.SearchInRadius(self.node_coordinates, radius)
        self.assertEqual(results.NumberOfGlobalResults(), 8)
        ids = results.GetResultIndices()
        for id in ids:
            self.assertTrue(id in cond_id_ref)

        # Nodes array search
        results = self.search.SearchInRadius(self.sub_model_part.Nodes, radius)
        self.assertEqual(results.NumberOfPointsResults(), 1)
        node_results = results[self.node_coordinates]
        node_results.SynchronizeAll(self.data_comm)   
        self.assertEqual(node_results.NumberOfGlobalResults(), 8)
        ids = node_results.GetResultIndices()
        self.assertEqual(len(ids), 8)
        for id in ids:
            self.assertTrue(id in cond_id_ref)

    def test_GeometricalObjectsBins_SearchNearestInRadius(self):
        radius = 0.35

        # One node search
        result = self.search.SearchNearestInRadius(self.node_coordinates, radius)
        self.assertEqual(result.NumberOfGlobalResults(), 1)
        ids = result.GetResultIndices()
        self.assertTrue(1 in ids)

        # Nodes array search
        results = self.search.SearchNearestInRadius(self.sub_model_part.Nodes, radius)   
        self.assertEqual(results.NumberOfPointsResults(), 1)
        node_results = results[self.node_coordinates]
        self.assertEqual(node_results.NumberOfGlobalResults(), 1)
        # Local result
        if node_results.NumberOfLocalResults() == 1:
            self.assertEqual(node_results[0].Id, 1)
        # Global result
        ids = node_results.GetResultIndices()
        self.assertEqual(len(ids), 1)
        self.assertTrue(1 in ids)

    def test_GeometricalObjectsBins_SearchNearest(self):
        # One node search
        result = self.search.SearchNearest(self.node_coordinates)
        self.assertEqual(result.NumberOfGlobalResults(), 1)
        ids = result.GetResultIndices()
        self.assertTrue(1 in ids)

        # Nodes array search
        results = self.search.SearchNearest(self.sub_model_part.Nodes) 
        self.assertEqual(results.NumberOfPointsResults(), 1)
        node_results = results[self.node_coordinates] 
        self.assertEqual(node_results.NumberOfGlobalResults(), 1)
        # Local result
        if node_results.NumberOfLocalResults() == 1:
            self.assertEqual(node_results[0].Id, 1)
        # Global result
        ids = node_results.GetResultIndices()
        self.assertEqual(len(ids), 1)
        self.assertTrue(1 in ids)

    def test_GeometricalObjectsBins_SearchIsInside(self):
        # One node search
        result = self.search.SearchIsInside(self.node_coordinates)
        self.assertFalse(result.IsObjectFound())

        # Nodes array search
        results = self.search.SearchIsInside(self.sub_model_part.Nodes) 
        self.assertEqual(results.NumberOfPointsResults(), 1)
        node_results = results[self.node_coordinates] 
        self.assertFalse(node_results.IsObjectFound())

if __name__ == '__main__':
    KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.WARNING)
    KratosUnittest.main()