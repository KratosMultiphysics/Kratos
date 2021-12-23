import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.mpi import DistributedModelPartInitializer

import os

def GetFilePath(fileName):
    return os.path.dirname(os.path.realpath(__file__)) + "/" + fileName

class TestDistributedModelPartInitializer(KratosUnittest.TestCase):

    def test_with_real_model(self):
        current_model = KM.Model()
        model_part = current_model.CreateModelPart("main_model_part")
        model_part.AddNodalSolutionStepVariable(KM.PARTITION_INDEX)
        model_part.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(KM.VISCOSITY)

        if KM.DataCommunicator.GetDefault().Rank() == 0:
            KM.ModelPartIO(GetFilePath("test_mpi_communicator")).ReadModelPart(model_part)

        DistributedModelPartInitializer(model_part, 0).Execute()

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


if __name__ == '__main__':
    KratosUnittest.main()
