import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.mpi.distributed_import_model_part_utility import DistributedImportModelPartUtility
from KratosMultiphysics.mpi import DataCommunicatorFactory
from KratosMultiphysics.kratos_utilities import DeleteDirectoryIfExisting

import os

def GetFilePath(fileName):
    return os.path.dirname(os.path.realpath(__file__)) + "/" + fileName

@KratosUnittest.skipIfApplicationsNotAvailable("MetisApplication")
class TestDistributedImportModelPartUtility(KratosUnittest.TestCase):

    def tearDown(self):
        # Remove the Metis partitioning files
        # this can only be done after all processes arrived here!
        # In these tests not all ranks participate in the test (=> SubComm),
        # hence they could remove the files for the other ranks before they could read!
        KM.Testing.GetDefaultDataCommunicator().Barrier()

        DeleteDirectoryIfExisting("test_mpi_communicator_partitioned")

        # next test can only start after all the processes arrived here, otherwise race conditions with deleting the files can occur
        KM.Testing.GetDefaultDataCommunicator().Barrier()

    def __execute_test(self, in_memory, all_ranks):
        settings = KM.Parameters("""{
            "model_import_settings" : {
                "input_type" : "mdpa"
            },
            "echo_level" : 0
        }""")
        settings["model_import_settings"].AddEmptyValue("input_filename").SetString(GetFilePath("test_files/mdpa_files/test_mpi_communicator"))
        settings["model_import_settings"].AddEmptyValue("partition_in_memory").SetBool(in_memory)

        data_comm_name = "World"
        if not all_ranks:
            default_data_comm = KM.ParallelEnvironment.GetDefaultDataCommunicator()
            size = default_data_comm.Size()
            if size < 3:
                self.skipTest("This test needs at least 3 mpi processes")

            ranks = [i for i in range(1,size)]
            data_comm_name = "AllExceptFirst"
            sub_comm = DataCommunicatorFactory.CreateFromRanksAndRegister(default_data_comm, ranks, data_comm_name)
            self.addCleanup(KM.ParallelEnvironment.UnregisterDataCommunicator, data_comm_name)

            settings["model_import_settings"].AddEmptyValue("data_communicator_name").SetString(data_comm_name)

            if default_data_comm.Rank() == 0:
                self.assertFalse(sub_comm.IsDefinedOnThisRank())
            else:
                self.assertTrue(sub_comm.IsDefinedOnThisRank())

            if not sub_comm.IsDefinedOnThisRank():
                # this rank does not participate
                return

        current_model = KM.Model()
        model_part = current_model.CreateModelPart("main_model_part")
        model_part.AddNodalSolutionStepVariable(KM.PARTITION_INDEX)
        model_part.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(KM.VISCOSITY)

        import_util = DistributedImportModelPartUtility(model_part, settings)
        import_util.ImportModelPart()
        import_util.CreateCommunicators()

        # check main ModelPart
        self.assertTrue(model_part.IsDistributed())

        self.assertEqual(model_part.GetCommunicator().GlobalNumberOfNodes(), 9)
        self.assertEqual(model_part.GetCommunicator().GlobalNumberOfElements(), 8)
        self.assertEqual(model_part.GetCommunicator().GlobalNumberOfConditions(), 8)

        self.assertEqual(model_part.NumberOfSubModelParts(), 1)
        self.assertTrue(model_part.HasSubModelPart("Skin"))

        # Check SubModelPart
        smp = model_part.GetSubModelPart("Skin")
        self.assertTrue(smp.IsDistributed())

        self.assertEqual(smp.GetCommunicator().GlobalNumberOfNodes(), 8)
        self.assertEqual(smp.GetCommunicator().GlobalNumberOfElements(), 0)
        self.assertEqual(smp.GetCommunicator().GlobalNumberOfConditions(), 8)

        self.assertEqual(smp.NumberOfSubModelParts(), 1)
        self.assertTrue(smp.HasSubModelPart("Top_side_1"))

        # Check SubSubModelPart
        sub_smp = smp.GetSubModelPart("Top_side_1")
        self.assertTrue(sub_smp.IsDistributed())

        self.assertEqual(sub_smp.GetCommunicator().GlobalNumberOfNodes(), 2)
        self.assertEqual(sub_smp.GetCommunicator().GlobalNumberOfElements(), 0)
        self.assertEqual(sub_smp.GetCommunicator().GlobalNumberOfConditions(), 1)

        self.assertEqual(sub_smp.NumberOfSubModelParts(), 0)


    def test_with_all_ranks(self):
        self.__execute_test(in_memory=False, all_ranks=True)

    def test_with_sub_comm(self):
        self.__execute_test(in_memory=False, all_ranks=False)

    def test_in_memory_with_all_ranks(self):
        self.__execute_test(in_memory=True, all_ranks=True)

    def test_in_memory_with_sub_comm(self):
        self.__execute_test(in_memory=True, all_ranks=False)


if __name__ == '__main__':
    KratosUnittest.main()
