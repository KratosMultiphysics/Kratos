import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.testing.utilities import ReadDistributedModelPart

import os

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)


class TestMPIModelPart(KratosUnittest.TestCase):

    def test_remove_nodes_parallel_interfaces(self):
        current_model = KratosMultiphysics.Model()
        main_model_part = current_model.CreateModelPart("MainModelPart")
        main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DENSITY)
        main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VISCOSITY)
        main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)

        main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 2)

        ReadDistributedModelPart(GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/test_mpi_communicator"), main_model_part)

        for node in main_model_part.Nodes:
            if node.Id % 2:
                node.Set(KratosMultiphysics.TO_ERASE)
            node.SetSolutionStepValue(KratosMultiphysics.DENSITY, main_model_part.GetCommunicator().MyPID())

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
