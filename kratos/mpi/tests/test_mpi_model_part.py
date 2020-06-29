import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.mpi as KratosMPI
import KratosMultiphysics.MetisApplication as KratosMetis # TODO refactor this by using the "distributed_import_model_part_utility"
import KratosMultiphysics.kratos_utilities as kratos_utilities

import os

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)


class TestMPICommunicator(KratosUnittest.TestCase):

    def setUp(self):
        self.communicator = KratosMultiphysics.DataCommunicator.GetDefault()

    def tearDown(self):
        rank = self.communicator.Rank()
        if rank == 0:
            kratos_utilities.DeleteFileIfExisting("test_mpi_communicator.time")
        kratos_utilities.DeleteFileIfExisting("test_mpi_communicator_"+str(rank)+".mdpa")
        kratos_utilities.DeleteFileIfExisting("test_mpi_communicator_"+str(rank)+".time")
        self.communicator.Barrier()

    def _read_model_part_mpi(self,main_model_part):

        if self.communicator.Size() == 1:
            self.skipTest("Test can be run only using more than one mpi process")

        ## Add variables to the model part
        main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DENSITY)
        main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VISCOSITY)
        main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PARTITION_INDEX)

        ## Serial partition of the original .mdpa file
        input_filename = GetFilePath("test_mpi_communicator")
        if self.communicator.Rank() == 0 :

            # Original .mdpa file reading
            model_part_io = KratosMultiphysics.ModelPartIO(input_filename)

            # Partition of the original .mdpa file
            number_of_partitions = self.communicator.Size() # Number of partitions equals the number of processors
            domain_size = main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
            verbosity = 0
            sync_conditions = True # Make sure that the condition goes to the same partition as the element is a face of

            partitioner = KratosMetis.MetisDivideHeterogeneousInputProcess(model_part_io, number_of_partitions , domain_size, verbosity, sync_conditions)
            partitioner.Execute()

            KratosMultiphysics.Logger.PrintInfo("TestMPICommunicator","Metis divide finished.")

        self.communicator.Barrier()

        ## Read the partitioned .mdpa files
        mpi_input_filename = input_filename + "_" + str(self.communicator.Rank())
        model_part_io = KratosMultiphysics.ModelPartIO(mpi_input_filename)
        model_part_io.ReadModelPart(main_model_part)

        ## Construct and execute the Parallel fill communicator
        ParallelFillCommunicator = KratosMPI.ParallelFillCommunicator(main_model_part.GetRootModelPart())
        ParallelFillCommunicator.Execute()

        ## Check submodelpart of each main_model_part of each processor
        self.assertTrue(main_model_part.HasSubModelPart("Skin"))
        skin_sub_model_part = main_model_part.GetSubModelPart("Skin")

    def test_remove_nodes_parallel_interfaces(self):
        current_model = KratosMultiphysics.Model()
        main_model_part = current_model.CreateModelPart("MainModelPart")
        main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 2)

        self._read_model_part_mpi(main_model_part)

        for node in main_model_part.Nodes:
            if node.Id % 2:
                node.Set(KratosMultiphysics.TO_ERASE)
            node.SetSolutionStepValue(KratosMultiphysics.DENSITY, self.communicator.Rank())

        main_model_part.RemoveNodesFromAllLevels(KratosMultiphysics.TO_ERASE)

        # Try to synchronize after deleting
        main_model_part.GetCommunicator().SynchronizeVariable(KratosMultiphysics.DENSITY)

        # Local nodes have no nodes with odd id:
        for node in main_model_part.Nodes:
            self.assertEqual(node.Id % 2, 0)

        # Local itnerfaces have no nodes with odd id:
        for node in main_model_part.GetCommunicator().LocalMesh().Nodes:
            self.assertEqual(node.Id % 2, 0)
            self.assertEqual(node.GetSolutionStepValue(KratosMultiphysics.DENSITY), node.GetSolutionStepValue(KratosMultiphysics.PARTITION_INDEX))

        # Ghost itnerfaces have no nodes with odd id:
        for node in main_model_part.GetCommunicator().GhostMesh().Nodes:
            self.assertEqual(node.Id % 2, 0)
            self.assertEqual(node.GetSolutionStepValue(KratosMultiphysics.DENSITY), node.GetSolutionStepValue(KratosMultiphysics.PARTITION_INDEX))

        # Local itnerfaces have no nodes with odd id:
        for node in main_model_part.GetCommunicator().InterfaceMesh().Nodes:
            self.assertEqual(node.Id % 2, 0)
            self.assertEqual(node.GetSolutionStepValue(KratosMultiphysics.DENSITY), node.GetSolutionStepValue(KratosMultiphysics.PARTITION_INDEX))

        # Local itnerfaces have no nodes with odd id:
        for i in range(0, main_model_part.GetCommunicator().GetNumberOfColors()):
            for node in main_model_part.GetCommunicator().LocalMesh(i).Nodes:
                self.assertEqual(node.Id % 2, 0)
                self.assertEqual(node.GetSolutionStepValue(KratosMultiphysics.DENSITY), node.GetSolutionStepValue(KratosMultiphysics.PARTITION_INDEX))

            for node in main_model_part.GetCommunicator().GhostMesh(i).Nodes:
                self.assertEqual(node.Id % 2, 0)
                self.assertEqual(node.GetSolutionStepValue(KratosMultiphysics.DENSITY), node.GetSolutionStepValue(KratosMultiphysics.PARTITION_INDEX))

            for node in main_model_part.GetCommunicator().InterfaceMesh(i).Nodes:
                self.assertEqual(node.Id % 2, 0)
                self.assertEqual(node.GetSolutionStepValue(KratosMultiphysics.DENSITY), node.GetSolutionStepValue(KratosMultiphysics.PARTITION_INDEX))


if __name__ == '__main__':
    KratosUnittest.main()
