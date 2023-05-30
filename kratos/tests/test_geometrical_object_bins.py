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
            self.search = KM.GeometricalObjectsBinsMPI(self.model_part.Conditions, self.data_comm)
        else:
            self.search = KM.GeometricalObjectsBins(self.model_part.Conditions)

        # Create node for search
        if KM.IsDistributedRun():
            import KratosMultiphysics.mpi as KratosMPI
            # Only added to first rank to actualy check it works in all ranks
            if self.data_comm.Rank() == 0:
                self.node = self.sub_model_part.CreateNewNode(100000, 0.0, 0.0, 0.15)
                self.node.SetSolutionStepValue(KM.PARTITION_INDEX, 0)
            ParallelFillCommunicator = KratosMPI.ParallelFillCommunicator(self.model_part)
            ParallelFillCommunicator.Execute()
        else:
            self.node = self.sub_model_part.CreateNewNode(100000, 0.0, 0.0, 0.15)

    def test_GeometricalObjectsBins_SearchInRadius(self):
        # Define radius
        radius = 0.35

        # Reference solution
        cond_id_ref = [125,78,117,18,68,1,41,119]

        # One node search
        results = self.search.SearchInRadius(self.node, radius)
        results.SynchronizeAll(self.data_comm)
        self.assertEqual(results.NumberOfGlobalResults(), 8)
        ids = results.GetResultIndices()
        for id in ids:
            self.assertTrue(id in cond_id_ref)

        # Nodes array search
        results = self.search.SearchInRadius(self.sub_model_part.Nodes, radius)
        self.assertEqual(results.NumberOfPointsResults(), 1)
        node_results = results[self.node]
        node_results.SynchronizeAll(self.data_comm)   
        self.assertEqual(node_results.NumberOfGlobalResults(), 8)
        ids = node_results.GetResultIndices()
        for id in ids:
            self.assertTrue(id in cond_id_ref)

    def test_GeometricalObjectsBins_SearchNearestInRadius(self):
        radius = 0.35

        # One node search
        result = self.search.SearchNearestInRadius(self.node, radius)
        self.assertEqual(result.Get().Id, 1)

        # Nodes array search
        results = self.search.SearchNearestInRadius(self.sub_model_part.Nodes, radius)   
        self.assertEqual(results.NumberOfPointsResults(), 1)
        node_results = results[self.node]
        node_results.SynchronizeAll(self.data_comm)    
        self.assertEqual(node_results.NumberOfGlobalResults(), 1)
        if self.data_comm.Rank() == 0:
            self.assertEqual(node_results[0].Id, 1) # Local result
        self.assertEqual(node_results(0).Id, 1)     # Global result

    def test_GeometricalObjectsBins_SearchNearest(self):
        # One node search
        result = self.search.SearchNearest(self.node)
        self.assertEqual(result.Get().Id, 1)

        # Nodes array search
        results = self.search.SearchNearest(self.sub_model_part.Nodes) 
        self.assertEqual(results.NumberOfPointsResults(), 1)
        node_results = results[self.node]
        node_results.SynchronizeAll(self.data_comm)    
        self.assertEqual(node_results.NumberOfGlobalResults(), 1)
        if self.data_comm.Rank() == 0:
            self.assertEqual(node_results[0].Id, 1) # Local result
        self.assertEqual(node_results(0).Id, 1)     # Global result

    def test_GeometricalObjectsBins_SearchIsInside(self):
        # One node search
        result = self.search.SearchIsInside(self.node)
        self.assertEqual(result.Get(), None)

        # Nodes array search
        results = self.search.SearchIsInside(self.sub_model_part.Nodes) 
        self.assertEqual(results.NumberOfPointsResults(), 1)
        node_results = results[self.node]
        node_results.SynchronizeAll(self.data_comm)    
        self.assertEqual(node_results.NumberOfGlobalResults(), 0)

if __name__ == '__main__':
    KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.WARNING)
    KratosUnittest.main()