import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest

import os

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

def ReadDistributedModelPart(model_part, mdpa_file_name):
    from KratosMultiphysics.mpi import distributed_import_model_part_utility
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PARTITION_INDEX)
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DENSITY)
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VISCOSITY)
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PARTITION_INDEX)

    importer_settings = KratosMultiphysics.Parameters("""{
        "model_import_settings": {
            "input_type": "mdpa",
            "input_filename": \"""" + mdpa_file_name + """\",
            "partition_in_memory" : true
        },
        "echo_level" : 0
    }""")

    model_part_import_util = distributed_import_model_part_utility.DistributedImportModelPartUtility(model_part, importer_settings)
    model_part_import_util.ImportModelPart()
    model_part_import_util.CreateCommunicators()

class TestMPIModelPart(KratosUnittest.TestCase):

    def setUp(self):
        self.communicator = KratosMultiphysics.DataCommunicator.GetDefault()


    def test_remove_nodes_parallel_interfaces(self):
        current_model = KratosMultiphysics.Model()
        main_model_part = current_model.CreateModelPart("MainModelPart")
        main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 2)

        ReadDistributedModelPart(main_model_part, GetFilePath("test_mpi_communicator"))

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
