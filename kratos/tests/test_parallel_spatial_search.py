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

def PostProcessVTK(model, entity_type = "element", rank = 0):
    import KratosMultiphysics.vtk_output_process as vtk_output_process
    vtk_output_parameters = KM.Parameters("""{
        "Parameters" : {
            "model_part_name"                    : "Main",
            "file_format"                        : "ascii",
            "output_precision"                   : 8,
            "output_interval"                    : 2,
            "output_sub_model_parts"             : true,
            "entity_type"                        : "element",
            "nodal_solution_step_data_variables" : ["PARTITION_INDEX"],
            "nodal_data_value_variables"         : [],
            "nodal_flags"                        : [],
            "element_data_value_variables"       : [],
            "element_flags"                      : [],
            "condition_data_value_variables"     : [],
            "condition_flags"                    : [],
            "output_path"                        : "test_vtk_output"
        }
    }""")
    vtk_output_parameters["Parameters"]["entity_type"].SetString(entity_type)
    vtk_output_parameters["Parameters"]["output_path"].SetString("test_vtk_output_" + str(rank))

    process = vtk_output_process.Factory(vtk_output_parameters, model)
    process.ExecuteInitialize()
    process.ExecuteInitializeSolutionStep()
    process.ExecuteFinalizeSolutionStep()
    process.PrintOutput()

def GetFilePath(fileName):
    """
    Get the absolute path of the given file name.
    """
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

def RemoveFiles(mdpa_name):
    """
    Remove the time-related files associated with a given .mdpa file.
    """
    kratos_utils.DeleteFileIfExisting(mdpa_name + ".time")

def CreateSearch(model_part, data_comm, search_type, entity_type = "Conditions"):
    """
    Create an appropriate search wrapper based on the given search type.
    """
    if search_type == "Octree":
        if entity_type == "Conditions":
            search = KM.ParallelSpatialSearchOCTreeCondition(model_part.Conditions, data_comm)
        elif entity_type == "Elements":
            search = KM.ParallelSpatialSearchOCTreeElement(model_part.Elements, data_comm)
        elif entity_type == "Nodes":
            search = KM.ParallelSpatialSearchOCTreeNode(model_part.Nodes, data_comm)
        else:
            raise Exception("Invalid entity type: " + entity_type)
    elif search_type == "KDTree":
        if entity_type == "Conditions":
            search = KM.ParallelSpatialSearchKDTreeCondition(model_part.Conditions, data_comm)
        elif entity_type == "Elements":
            search = KM.ParallelSpatialSearchKDTreeElement(model_part.Elements, data_comm)
        elif entity_type == "Nodes":
            search = KM.ParallelSpatialSearchKDTreeNode(model_part.Nodes, data_comm)
        else:
            raise Exception("Invalid entity type: " + entity_type)
    elif search_type == "StaticBinsTree":
        if entity_type == "Conditions":
            search = KM.ParallelSpatialSearchStaticBinsTreeCondition(model_part.Conditions, data_comm)
        elif entity_type == "Elements":
            search = KM.ParallelSpatialSearchStaticBinsTreeElement(model_part.Elements, data_comm)
        elif entity_type == "Nodes":
            search = KM.ParallelSpatialSearchStaticBinsTreeNode(model_part.Nodes, data_comm)
        else:
            raise Exception("Invalid entity type: " + entity_type)
    elif search_type == "DynamicBins":
        if entity_type == "Conditions":
            search = KM.ParallelSpatialSearchDynamicBinsCondition(model_part.Conditions, data_comm)
        elif entity_type == "Elements":
            search = KM.ParallelSpatialSearchDynamicBinsElement(model_part.Elements, data_comm)
        elif entity_type == "Nodes":
            search = KM.ParallelSpatialSearchDynamicBinsNode(model_part.Nodes, data_comm)
        else:
            raise Exception("Invalid entity type: " + entity_type)
    elif search_type == "GeometricalObjectBins":
        if entity_type == "Conditions":
            search = KM.ParallelSpatialSearchGeometricalObjectBins(model_part.Conditions, data_comm)
        elif entity_type == "Elements":
            search = KM.ParallelSpatialSearchGeometricalObjectBins(model_part.Elements, data_comm)
        else:
            raise Exception("Invalid entity type: " + entity_type)
    else:
        raise Exception("Invalid search type: " + search_type)

    return search

class TestParallelSpatialSearch(KratosUnittest.TestCase):
    """
    A test class for the ParallelSpatialSearchGeometricalObjectBins class.

    This class is used to test the functionality of the GeometricalObjectBins search wrapper
    for searching and querying geometric objects within a specified radius.

    Attributes:
        current_model (KM.Model): The Kratos Model used for testing.
        model_part (KM.ModelPart): The main model part for testing.
        sub_model_part (KM.ModelPart): The sub-model part for testing.
        mdpa_name (str): The name of the MDPA file used for testing.
        data_comm (KM.DataCommunicator): The data communicator for parallel execution.
        node_id (int): The ID of the test node.
        node (KM.Node): The test node used for searching.
        search (KM.ParallelSpatialSearchGeometricalObjectBins): The search wrapper instance for testing.

    Methods:
        setUpClass(cls):
            Set up the model part and its related entities for the test.

        tearDownClass(cls):
            Clean up after all tests are run.

        setUp(self):
            Set up for each individual test.

        test_ParallelSpatialSearchGeometricalObjectBins_SearchInRadius(self):
            Test for the 'SearchInRadius' method of the GeometricalObjectBins search wrapper.

        test_ParallelSpatialSearchGeometricalObjectBins_SearchNearestInRadius(self):
            Test for the 'SearchNearestInRadius' method of the GeometricalObjectBins search wrapper.

        test_ParallelSpatialSearchGeometricalObjectBins_SearchNearest(self):
            Test for the 'SearchNearest' method of the GeometricalObjectBins search wrapper.

        test_ParallelSpatialSearchGeometricalObjectBins_SearchIsInside(self):
            Test for the 'SearchIsInside' method of the GeometricalObjectBins search wrapper.
    """
    @classmethod
    def setUpClass(cls):
        """
        Setting up the model part and its related entities for the test.
        """
        # Create model and model parts
        cls.current_model = KM.Model()
        cls.model_part = cls.current_model.CreateModelPart("Main")
        cls.sub_model_part = cls.model_part.CreateSubModelPart("SubModelPart")

        # Define model part variables
        cls.model_part.ProcessInfo[KM.DOMAIN_SIZE] = 3
        cls.model_part.AddNodalSolutionStepVariable(KM.BULK_MODULUS)
        cls.model_part.AddNodalSolutionStepVariable(KM.NODAL_VAUX)
        cls.model_part.AddNodalSolutionStepVariable(KM.EXTERNAL_FORCES_VECTOR)
        cls.model_part.AddNodalSolutionStepVariable(KM.LOCAL_AXES_MATRIX)

        # If distributed run, add PARTITION_INDEX to model part
        if KM.IsDistributedRun():
            cls.model_part.AddNodalSolutionStepVariable(KM.PARTITION_INDEX)

        # Define mdpa file name and read it into the model part
        cls.mdpa_name = GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/coarse_sphere_with_conditions")
        ReadModelPart(cls.mdpa_name, cls.model_part)

    @classmethod
    def tearDownClass(cls):
        """
        Cleanup after all tests are run.
        """
        RemoveFiles(cls.mdpa_name)

    def setUp(self):
        """
        Setup for each individual test.
        """
        # Create search
        self.data_comm = self.model_part.GetCommunicator().GetDataCommunicator()

        # Create node for search
        self.node_id = 100
        if KM.IsDistributedRun():
            # Only added to first rank to actually check it works in all ranks
            if self.data_comm.Rank() == 0:
                self.node = self.sub_model_part.CreateNewNode(self.node_id, 0.0, 0.0, 0.15)
                self.node.SetSolutionStepValue(KM.PARTITION_INDEX, 0)
            ParallelFillCommunicator = KratosMPI.ParallelFillCommunicator(self.model_part, self.data_comm)
            ParallelFillCommunicator.Execute()
        else:
            self.node = self.sub_model_part.CreateNewNode(self.node_id, 0.0, 0.0, 0.15)

    def _SearchInRadius(self, search_type = "GeometricalObjectBins"):
        """
        Test for the 'SearchInRadius' method
        """
        # Create search
        self.search = CreateSearch(self.model_part, self.data_comm, search_type)

        # Define radius
        radius = 0.35

        # Reference solution
        cond_id_ref = [125,78,117,18,68,1,41,119]

        # Nodes array search
        results = self.search.SearchInRadius(self.sub_model_part.Nodes, radius)
        # Only in partitions were results are found
        number_search_results = results.NumberOfSearchResults()
        check = self.data_comm.SumAll(number_search_results)
        self.assertGreater(check, 0)
        if number_search_results > 0:
            self.assertEqual(number_search_results, 1)
            node_results = results[0]
            self.assertEqual(node_results.NumberOfGlobalResults(), 8)
            ids = results.GetResultIndices()[0]
            self.assertEqual(len(ids), 8)
            for id in ids:
                self.assertTrue(id in cond_id_ref)

    def _SearchNearestInRadius(self, search_type = "GeometricalObjectBins"):
        """
        Test for the 'SearchNearestInRadius' method
        """
        # Create search
        self.search = CreateSearch(self.model_part, self.data_comm, search_type)

        # Define radius
        radius = 0.35

        # Nodes array search
        results = self.search.SearchNearestInRadius(self.sub_model_part.Nodes, radius)
        # Only in partitions were results are found
        number_search_results = results.NumberOfSearchResults()
        check = self.data_comm.SumAll(number_search_results)
        self.assertGreater(check, 0)
        if number_search_results > 0:
            self.assertEqual(number_search_results, 1)
            node_results = results[0]
            self.assertEqual(node_results.NumberOfGlobalResults(), 1)
            # Local result
            if node_results.NumberOfLocalResults() == 1:
                self.assertTrue(node_results[0].IsObjectFound())
            # Global result
            ids = results.GetDistances()[0]
            self.assertEqual(len(ids), 1)
            # Distance is different depending on the search type used (GeometricalObjectBins or not)
            if search_type == "GeometricalObjectBins":
                self.assertAlmostEqual(ids[0], 0.326335, 6)
            else:
                self.assertAlmostEqual(ids[0], 0.106584, 6)

    def _SearchNearest(self, search_type = "GeometricalObjectBins"):
        """
        Test for the 'SearchNearest' method
        """
        # Create search
        self.search = CreateSearch(self.model_part, self.data_comm, search_type)

        # Nodes array search
        results = self.search.SearchNearest(self.sub_model_part.Nodes)
        # Only in partitions were results are found
        number_search_results = results.NumberOfSearchResults()
        check = self.data_comm.SumAll(number_search_results)
        self.assertGreater(check, 0)
        if number_search_results > 0:
            self.assertEqual(number_search_results, 1)
            node_results = results[0]
            self.assertEqual(node_results.NumberOfGlobalResults(), 1)
            # Local result
            if node_results.NumberOfLocalResults() == 1:
                self.assertTrue(node_results[0].IsObjectFound())
            # Global result
            ids = results.GetDistances()[0]
            self.assertEqual(len(ids), 1)
            # Distance is different depending on the search type used (GeometricalObjectBins or not)
            if search_type == "GeometricalObjectBins":
                self.assertAlmostEqual(ids[0], 0.326335, 6)
            else:
                self.assertAlmostEqual(ids[0], 0.106584, 6)

    def test_ParallelSpatialSearchGeometricalObjectBins_SearchInRadius(self):
        """
        Test for the 'SearchInRadius' method of the GeometricalObjectBins search wrapper.
        """
        self._SearchInRadius("GeometricalObjectBins")

    def test_ParallelSpatialSearchGeometricalObjectBins_SearchNearestInRadius(self):
        """
        Test for the 'SearchNearestInRadius' method of the GeometricalObjectBins search wrapper.
        """
        self._SearchNearestInRadius("GeometricalObjectBins")

    def test_ParallelSpatialSearchGeometricalObjectBins_SearchNearest(self):
        """
        Test for the 'SearchNearest' method of the GeometricalObjectBins search wrapper.
        """
        self._SearchNearest("GeometricalObjectBins")

    def test_ParallelSpatialSearchGeometricalObjectBins_SearchIsInside(self):
        """
        Test for the 'SearchIsInside' method of the GeometricalObjectBins search wrapper.
        """
        # Create search
        self.search = CreateSearch(self.model_part, self.data_comm, "GeometricalObjectBins")

        # Nodes array search
        results = self.search.SearchIsInside(self.sub_model_part.Nodes)
        # Only in partitions were results are found
        number_search_results = results.NumberOfSearchResults()
        check = self.data_comm.SumAll(number_search_results)
        self.assertGreater(check, 0)
        if number_search_results > 0:
            self.assertEqual(number_search_results, 1)
            node_results = results[0]
            self.assertFalse(node_results.IsObjectFound())

    def test_ParallelSpatialSearchKDTree_SearchInRadius(self):
        """
        Test for the 'SearchInRadius' method of the KDTree search wrapper.
        """
        self._SearchInRadius("KDTree")

    def test_ParallelSpatialSearchKDTree_SearchNearestInRadius(self):
        """
        Test for the 'SearchNearestInRadius' method of the KDTree search wrapper.
        """
        self._SearchNearestInRadius("KDTree")

    def test_ParallelSpatialSearchKDTree_SearchNearest(self):
        """
        Test for the 'SearchNearest' method of the KDTree search wrapper.
        """
        self._SearchNearest("KDTree")

    def test_ParallelSpatialSearchOctree_SearchInRadius(self):
        """
        Test for the 'SearchInRadius' method of the Octree search wrapper.
        """
        self._SearchInRadius("Octree")

    def test_ParallelSpatialSearchOctree_SearchNearestInRadius(self):
        """
        Test for the 'SearchNearestInRadius' method of the Octree search wrapper.
        """
        self._SearchNearestInRadius("Octree")

    def test_ParallelSpatialSearchOctree_SearchNearest(self):
        """
        Test for the 'SearchNearest' method of the Octree search wrapper.
        """
        self._SearchNearest("Octree")

    def test_ParallelSpatialSearchDynamicBins_SearchInRadius(self):
        """
        Test for the 'SearchInRadius' method of the DynamicBins search wrapper.
        """
        self._SearchInRadius("DynamicBins")

    def test_ParallelSpatialSearchDynamicBins_SearchNearestInRadius(self):
        """
        Test for the 'SearchNearestInRadius' method of the DynamicBins search wrapper.
        """
        self._SearchNearestInRadius("DynamicBins")

    def test_ParallelSpatialSearchDynamicBins_SearchNearest(self):
        """
        Test for the 'SearchNearest' method of the DynamicBins search wrapper.
        """
        self._SearchNearest("DynamicBins")

    def test_ParallelSpatialSearchStaticBinsTree_SearchInRadius(self):
        """
        Test for the 'SearchInRadius' method of the StaticBinsTree search wrapper.
        """
        self._SearchInRadius("StaticBinsTree")

    def test_ParallelSpatialSearchStaticBinsTree_SearchNearestInRadius(self):
        """
        Test for the 'SearchNearestInRadius' method of the StaticBinsTree search wrapper.
        """
        self._SearchNearestInRadius("StaticBinsTree")

    def test_ParallelSpatialSearchStaticBinsTree_SearchNearest(self):
        """
        Test for the 'SearchNearest' method of the StaticBinsTree search wrapper.
        """
        self._SearchNearest("StaticBinsTree")

class TestParallelSpatialSearchSmallSquare(KratosUnittest.TestCase):
    """

    """
    @classmethod
    def setUpClass(cls):
        """
        Setting up the model part and its related entities for the test.
        """
        # Create model and model parts
        cls.current_model = KM.Model()
        cls.model_part = cls.current_model.CreateModelPart("Main")
        cls.sub_model_part = cls.model_part.CreateSubModelPart("SubModelPart")

        # Define model part variables
        cls.model_part.ProcessInfo[KM.DOMAIN_SIZE] = 3

        # If distributed run, add PARTITION_INDEX to model part
        if KM.IsDistributedRun():
            cls.model_part.AddNodalSolutionStepVariable(KM.PARTITION_INDEX)

        # Define mdpa file name and read it into the model part
        cls.mdpa_name = GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/small_square")
        ReadModelPart(cls.mdpa_name, cls.model_part)

    @classmethod
    def tearDownClass(cls):
        """
        Cleanup after all tests are run.
        """
        RemoveFiles(cls.mdpa_name)

    def setUp(self):
        """
        Setup for each individual test.
        """
        # Create search
        self.data_comm = self.model_part.GetCommunicator().GetDataCommunicator()

        # Adding nodes in the border x = 0.0
        for node in self.model_part.Nodes:
            if node.X == 0.0:
                self.sub_model_part.AddNode(node)

        # # Debug mesh
        # PostProcessVTK(self.current_model, "element", self.data_comm.Rank())

    def _SearchInRadius(self, entity_type = "Nodes"):
        """
        Test for the 'SearchInRadius' method
        """
        # Create search
        self.search = CreateSearch(self.model_part, self.data_comm, "KDTree", entity_type)

        # Define radius
        radius = 0.5

        # Reference solution
        if entity_type == "Nodes":
            id_ref = {
                1 :  [1, 2, 3, 4, 5],
                2 :  [1, 2, 3, 4, 5, 6, 9],
                6 :  [2, 4, 5, 6, 9, 12],
                12 : [5, 6, 9, 12, 14, 16, 20],
                16 : [9, 12, 14, 16, 20],
            }
        elif entity_type == "Elements":
            id_ref = {
                1 :  [4, 17, 18, 29, 30],
                2 :  [4, 17, 18, 20, 24, 27, 28, 29, 30, 31],
                6 :  [1, 4, 11, 18, 24, 25, 26, 27, 28, 29, 30, 31],
                12 : [1, 9, 11, 24, 25, 26, 27, 28, 29, 30],
                16 : [1, 9, 11, 24, 25, 27],
            }

        # Nodes array search
        results = self.search.SearchInRadius(self.sub_model_part.Nodes, radius)
        # Only in partitions were results are found
        number_search_results = results.NumberOfSearchResults()
        check = self.data_comm.SumAll(number_search_results)
        self.assertGreater(check, 0)
        if number_search_results > 0:
            self.assertEqual(number_search_results, 5)
            all_ids = results.GetResultIndices()
            for i in range(5):
                node_results = results[i]
                global_id = node_results.GetGlobalIndex()
                ids = all_ids[i]
                number_of_global_results = node_results.NumberOfGlobalResults()
                if global_id > 0: # Solution defined in this rank
                    ref_ids = id_ref[global_id]
                    self.assertEqual(number_of_global_results, len(ref_ids))
                    self.assertEqual(len(ids), len(ref_ids))
                    for id in ids:
                        self.assertTrue(id in ref_ids)

    def test_ParallelSpatialSearchKDTree_SearchInRadius_Nodes(self):
        """
        Test for the 'SearchInRadius' method of the KDTree search wrapper for nodes.
        """
        self._SearchInRadius("Nodes")

    def test_ParallelSpatialSearchKDTree_SearchInRadius_Elements(self):
        """
        Test for the 'SearchInRadius' method of the KDTree search wrapper for elements.
        """
        self._SearchInRadius("Elements")

if __name__ == '__main__':
    # Set logging severity and start the unittest
    KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.WARNING)
    KratosUnittest.main()