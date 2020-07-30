import KratosMultiphysics
import KratosMultiphysics.mpi as KratosMPI
from KratosMultiphysics.mpi import distributed_import_model_part_utility
import KratosMultiphysics.kratos_utilities as kratos_utilities
import KratosMultiphysics.KratosUnittest as UnitTest
import sys
import os

dependencies_are_available = kratos_utilities.CheckIfApplicationsAvailable("MetisApplication")

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

# Use cPickle on Python 2.7 (Note that only the cPickle module is supported on Python 2.7)
# Source: https://pybind11.readthedocs.io/en/stable/advanced/classes.html
pickle_message = ""
try:
    import cPickle as pickle
    have_pickle_module = True
except ImportError:
    if sys.version_info > (3, 0):
        try:
            import pickle
            have_pickle_module = True
        except ImportError:
            have_pickle_module = False
            pickle_message = "No pickle module found"
    else:
        have_pickle_module = False
        pickle_message = "No valid pickle module found"

def executeComputeArea_Task(pickled_model):
    # Unpickling model
    serialized_model = pickle.loads(pickled_model)

    # Unserializing model
    deserialized_model = KratosMultiphysics.Model()
    serialized_model.Load("ModelSerialization", deserialized_model)

    model_part = deserialized_model.GetModelPart("ModelPart").GetRootModelPart()
    ## Construct and execute the Parallel fill communicator
    ParallelFillCommunicator = KratosMPI.ParallelFillCommunicator(model_part)
    ParallelFillCommunicator.Execute()
    communicator = model_part.GetCommunicator()

    # Computing local areas
    for node in model_part.Nodes:
        node.SetValue(KratosMultiphysics.NODAL_AREA, 0.0)

    for elem in model_part.Elements:
        for node in elem.GetNodes():
            current_nodal_area = node.GetValue(KratosMultiphysics.NODAL_AREA)
            current_nodal_area += 1/3*elem.GetGeometry().Area()
            node.SetValue(KratosMultiphysics.NODAL_AREA, current_nodal_area)

    # Assembling nodal values
    communicator.AssembleNonHistoricalData(KratosMultiphysics.NODAL_AREA)

    # Computing sum of total area to check results.
    local_sum = sum(node.GetValue(KratosMultiphysics.NODAL_AREA) for node in model_part.Nodes if node.GetSolutionStepValue(KratosMultiphysics.PARTITION_INDEX) == communicator.GetDataCommunicator().Rank())
    total_sum = communicator.GetDataCommunicator().SumAll(local_sum)

    return total_sum

class TestMPISerializer(UnitTest.TestCase):

    def setUp(self):
        communicator = KratosMultiphysics.DataCommunicator.GetDefault()
        rank = communicator.Rank()

        current_model = KratosMultiphysics.Model()
        model_part_to_read = current_model.CreateModelPart("ModelPart")
        model_part_to_read.AddNodalSolutionStepVariable(KratosMultiphysics.PARTITION_INDEX)

        importer_settings = KratosMultiphysics.Parameters("""{
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": "",
                "partition_in_memory" : true
            },
            "echo_level" : 0
        }""")
        mdpa_name = GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/test_mpi_serializer")
        importer_settings["model_import_settings"]["input_filename"].SetString(mdpa_name)

        model_part_import_util = distributed_import_model_part_utility.DistributedImportModelPartUtility(model_part_to_read, importer_settings)
        model_part_import_util.ImportModelPart()

        ## Serialize each model_part
        serialized_model = KratosMultiphysics.MpiSerializer()
        serialized_model.Save("ModelSerialization",current_model)

        self.pickled_data = pickle.dumps(serialized_model, protocol=2) # Second argument is the protocol and is NECESSARY (according to pybind11 docs)

    def tearDown(self):
        communicator = KratosMultiphysics.DataCommunicator.GetDefault()
        rank = communicator.Rank()
        kratos_utilities.DeleteFileIfExisting("test_mpi_serializer.time")
        kratos_utilities.DeleteFileIfExisting("test_mpi_serializer_"+str(rank)+".mdpa")
        kratos_utilities.DeleteFileIfExisting("test_mpi_serializer_"+str(rank)+".time")

    @UnitTest.skipUnless(have_pickle_module, "Pickle module error: " + pickle_message)
    @UnitTest.skipUnless(dependencies_are_available,"MetisApplication is not available")
    def testCalculateNodalArea(self):
        total_nodal_area = executeComputeArea_Task(self.pickled_data)
        # Check result
        self.assertAlmostEqual(total_nodal_area, 1.0)

if __name__ == "__main__":
    UnitTest.main()
