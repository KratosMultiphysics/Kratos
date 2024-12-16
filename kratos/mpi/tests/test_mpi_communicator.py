import KratosMultiphysics
import KratosMultiphysics.mpi as KratosMPI
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.testing.utilities import ReadModelPart
from KratosMultiphysics.kratos_utilities import DeleteDirectoryIfExisting

import os

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)
class TestMPICommunicator(KratosUnittest.TestCase):

    def tearDown(self):
        DeleteDirectoryIfExisting("test_mpi_communicator_partitioned")

    def _read_model_part_mpi(self,main_model_part):

        if KratosMultiphysics.Testing.GetDefaultDataCommunicator().Size() == 1:
            self.skipTest("Test can be run only using more than one mpi process")

        ## Add variables to the model part
        main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DENSITY)
        main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VISCOSITY)
        main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PARTITION_INDEX)

        ReadModelPart(GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/test_mpi_communicator"), main_model_part)

        ## Check submodelpart of each main_model_part of each processor
        self.assertTrue(main_model_part.HasSubModelPart("Skin"))

    def test_assemble_variable_in_model_part(self):
        current_model = KratosMultiphysics.Model()
        main_model_part = current_model.CreateModelPart("MainModelPart")
        main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 2)

        self._read_model_part_mpi(main_model_part)

        ## Initialize DENSITY variable
        for node in main_model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.DENSITY, 0, 0.0)
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X, 0, 0.0)
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y, 0, 0.0)
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z, 0, 0.0)

        ## Fill the DENSITY value in each processor
        for elem in main_model_part.Elements:
            for node in elem.GetNodes():
                d = node.GetSolutionStepValue(KratosMultiphysics.DENSITY)
                d += 1
                node.SetSolutionStepValue(KratosMultiphysics.DENSITY,0,d)
                node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X,0,d)
                node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y,0,d)
                node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z,0,d)

        ## Communicate between processors
        main_model_part.GetCommunicator().AssembleCurrentData(KratosMultiphysics.DENSITY)
        main_model_part.GetCommunicator().AssembleCurrentData(KratosMultiphysics.DISPLACEMENT)

        ## Expected nodal values after communication
        expected_values = {1 : 2,
                           2 : 3,
                           3 : 1,
                           4 : 3,
                           5 : 6,
                           6 : 3,
                           7 : 1,
                           8 : 3,
                           9 : 2}

        ## Check the obtained values
        for node in main_model_part.Nodes:
            self.assertEqual(node.GetSolutionStepValue(KratosMultiphysics.DENSITY), expected_values[node.Id])
            self.assertEqual(node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X), expected_values[node.Id])
            self.assertEqual(node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y), expected_values[node.Id])
            self.assertEqual(node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z), expected_values[node.Id])


    def test_assemble_variable_in_sub_model_part(self):
        current_model = KratosMultiphysics.Model()
        main_model_part = current_model.CreateModelPart("MainModelPart")
        main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 2)

        self._read_model_part_mpi(main_model_part)

        submodelpart = main_model_part.GetSubModelPart("Skin")

        # self.communicator.Barrier()

        # Check partitioning
        #~ for i in range(KratosMPI.mpi.size):
            #~ if(KratosMPI.mpi.rank == i):
                #~ print(" ")
                #~ print("Submodelpart elements of processor id = ", i)
                #~ for elem in submodelpart.Elements:
                    #~ print(elem.Id)
                #~ print("Submodelpart conditions of processor id = ", i)
                #~ for elem in submodelpart.Conditions:
                    #~ print(elem.Id)
                #~ print("Submodelpart local nodes of processor id = ", i)
                #~ for node in submodelpart.GetCommunicator().LocalMesh().Nodes:
                    #~ print(node.Id)
                #~ print("Submodelpart ghost nodes of processor id = ", i)
                #~ for node in submodelpart.GetCommunicator().GhostMesh().Nodes:
                    #~ print(node.Id)
                #~ print("Submodelpart interface nodes of processor id = ", i)
                #~ for node in submodelpart.GetCommunicator().InterfaceMesh().Nodes:
                    #~ print(node.Id)
                #~ print(" ")
            #~ KratosMPI.mpi.world.barrier()

        for node in submodelpart.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.DENSITY, 0, 0.0)
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X, 0, 0.0)
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y, 0, 0.0)
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z, 0, 0.0)

        for condition in submodelpart.Conditions:
            for node in condition.GetNodes():
                d = node.GetSolutionStepValue(KratosMultiphysics.DENSITY)
                d += 1
                node.SetSolutionStepValue(KratosMultiphysics.DENSITY,0,d)
                node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X, 0, d)
                node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y, 0, d)
                node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z, 0, d)

        submodelpart.GetCommunicator().AssembleCurrentData(KratosMultiphysics.DENSITY)
        submodelpart.GetCommunicator().AssembleCurrentData(KratosMultiphysics.DISPLACEMENT)

        for node in submodelpart.Nodes:
            self.assertEqual(node.GetSolutionStepValue(KratosMultiphysics.DENSITY), 2.0)
            self.assertEqual(node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X), 2.0)
            self.assertEqual(node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y), 2.0)
            self.assertEqual(node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z), 2.0)

    def testCommunicatorRankSize(self):
        current_model = KratosMultiphysics.Model()
        main_model_part = current_model.CreateModelPart("MainModelPart")
        main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 2)

        self._read_model_part_mpi(main_model_part)

        self.world = KratosMultiphysics.ParallelEnvironment.GetDataCommunicator("World")
        self.assertEqual(self.world.Rank(), main_model_part.GetCommunicator().MyPID())
        self.assertEqual(self.world.Size(), main_model_part.GetCommunicator().TotalProcesses())

        submodelpart = main_model_part.GetSubModelPart("Skin")
        self.assertEqual(self.world.Rank(), submodelpart.GetCommunicator().MyPID())
        self.assertEqual(self.world.Size(), submodelpart.GetCommunicator().TotalProcesses())

    def testCommunicatorReduction(self):
        current_model = KratosMultiphysics.Model()
        main_model_part = current_model.CreateModelPart("MainModelPart")
        main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 2)

        self._read_model_part_mpi(main_model_part)

        self.world = KratosMultiphysics.ParallelEnvironment.GetDataCommunicator("World")

        comm = main_model_part.GetCommunicator().GetDataCommunicator()

        self.assertEqual(comm.SumAll(1), self.world.Size())
        self.assertEqual(comm.SumAll(2.0), 2.0*self.world.Size())
        self.assertEqual(comm.MinAll(comm.Rank()), self.world.MinAll(self.world.Rank()))
        self.assertEqual(comm.MaxAll(comm.Rank()), self.world.MaxAll(self.world.Rank()))
        self.assertEqual(comm.ScanSum(1), self.world.Rank()+1)

    def testCommunicatorReductionSerial(self):
        current_model = KratosMultiphysics.Model()
        main_model_part = current_model.CreateModelPart("MainModelPart")

        # this one is not set, so it should be a serial Communicator
        comm = main_model_part.GetCommunicator().GetDataCommunicator()

        self.assertEqual(comm.SumAll(1), 1)
        self.assertEqual(comm.SumAll(2.0), 2.0)
        self.assertEqual(comm.MinAll(comm.Rank()), 0)
        self.assertEqual(comm.MaxAll(comm.Rank()), 0)
        self.assertEqual(comm.ScanSum(1), 1)

    def test_GlobalNumberOf_Methods(self):
        current_model = KratosMultiphysics.Model()
        main_model_part = current_model.CreateModelPart("MainModelPart")
        main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 2)

        self._read_model_part_mpi(main_model_part)

        # Adding new nodes for the constraint in rank 0 always
        data_comm = main_model_part.GetCommunicator().GetDataCommunicator()
        if data_comm.Rank() == 0:
            node100 = main_model_part.CreateNewNode(100, 0.00000, 1.00000, 0.00000)
            node100.SetSolutionStepValue(KratosMultiphysics.PARTITION_INDEX, 0)
            node101 = main_model_part.CreateNewNode(101, 0.00000, 0.50000, 0.00000)
            node101.SetSolutionStepValue(KratosMultiphysics.PARTITION_INDEX, 0)

        # Adding Dofs
        dofs_list = ["DISPLACEMENT_X", "DISPLACEMENT_Y", "DISPLACEMENT_Z"]
        KratosMultiphysics.VariableUtils.AddDofsList(dofs_list, main_model_part)

        # Create constraint in rank 0 always
        if data_comm.Rank() == 0:
            main_model_part.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", 1, node100, KratosMultiphysics.DISPLACEMENT_X, node101, KratosMultiphysics.DISPLACEMENT_X, 1.0, 0)

        ParallelFillCommunicator = KratosMPI.ParallelFillCommunicator(main_model_part)
        ParallelFillCommunicator.Execute()

        main_comm = main_model_part.GetCommunicator()

        self.assertEqual(main_comm.GlobalNumberOfNodes(), 11)
        self.assertEqual(main_comm.GlobalNumberOfElements(), 8)
        self.assertEqual(main_comm.GlobalNumberOfConditions(), 8)
        self.assertEqual(main_comm.GlobalNumberOfMasterSlaveConstraints(), 1)

        sub_comm = main_model_part.GetSubModelPart("Skin").GetCommunicator()

        self.assertEqual(sub_comm.GlobalNumberOfNodes(), 8)
        self.assertEqual(sub_comm.GlobalNumberOfElements(), 0)
        self.assertEqual(sub_comm.GlobalNumberOfConditions(), 8)
        self.assertEqual(sub_comm.GlobalNumberOfMasterSlaveConstraints(), 0)

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()
