from __future__ import print_function, absolute_import, division

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest

import sys

class TestCommunicator(KratosUnittest.TestCase):

    def setUp(self):
        self.model = KratosMultiphysics.Model()

    def testCommunicatorSetUp(self):
        model_part = self.model.CreateModelPart("Test_model_part",1)
        self.assertTrue("Communicator" in model_part.GetCommunicator().__str__())

    def testCommunicatorMyPID(self):
        model_part = self.model.CreateModelPart("Test_model_part",1)
        self.assertEqual(model_part.GetCommunicator().MyPID(), 0)

    def testCommunicatorTotalProcesses(self):
        model_part = self.model.CreateModelPart("Test_model_part",1)
        self.assertEqual(model_part.GetCommunicator().TotalProcesses(), 1)

    def testCommunicatorGetNumberOfColors(self):
        model_part = self.model.CreateModelPart("Test_model_part",1)
        self.assertEqual(model_part.GetCommunicator().GetNumberOfColors(), 1)

    def testCommunicatorSumAll(self):
        model_part = self.model.CreateModelPart("Test_model_part",1)
        a = 0
        self.assertEqual(model_part.GetCommunicator().SumAll(a), a)
        b = 0.0
        self.assertEqual(model_part.GetCommunicator().SumAll(b), b)
        c = KratosMultiphysics.Array3()
        c[0] = 0.0
        c[1] = 0.0
        c[2] = 0.0
        d = model_part.GetCommunicator().SumAll(c)
        for i in range(3):
            self.assertEqual(d[i], 0.0)

    def testCommunicatorMinAll(self):
        model_part = self.model.CreateModelPart("Test_model_part",1)
        a = 0
        self.assertEqual(model_part.GetCommunicator().MinAll(a), a)
        b = 0.0
        self.assertEqual(model_part.GetCommunicator().MinAll(b), b)

    def testCommunicatorMaxAll(self):
        model_part = self.model.CreateModelPart("Test_model_part",1)
        a = 0
        self.assertEqual(model_part.GetCommunicator().MaxAll(a), a)
        b = 0.0
        self.assertEqual(model_part.GetCommunicator().MaxAll(b), b)

    def testCommunicatorSynchronizeVariable(self):
        model_part = self.model.CreateModelPart("Test_model_part",1)
        self.assertEqual(model_part.GetCommunicator().SynchronizeVariable(KratosMultiphysics.IS_RESTARTED), True)
        self.assertEqual(model_part.GetCommunicator().SynchronizeVariable(KratosMultiphysics.DOMAIN_SIZE), True)
        self.assertEqual(model_part.GetCommunicator().SynchronizeVariable(KratosMultiphysics.TIME), True)
        self.assertEqual(model_part.GetCommunicator().SynchronizeVariable(KratosMultiphysics.VECTOR_LAGRANGE_MULTIPLIER), True)
        self.assertEqual(model_part.GetCommunicator().SynchronizeVariable(KratosMultiphysics.RECOVERED_STRESS), True)
        self.assertEqual(model_part.GetCommunicator().SynchronizeVariable(KratosMultiphysics.GREEN_LAGRANGE_STRAIN_TENSOR), True)

    def testCommunicatorSynchronizeNonHistoricalVariable(self):
        model_part = self.model.CreateModelPart("Test_model_part",1)
        self.assertEqual(model_part.GetCommunicator().SynchronizeNonHistoricalVariable(KratosMultiphysics.IS_RESTARTED), True)
        self.assertEqual(model_part.GetCommunicator().SynchronizeNonHistoricalVariable(KratosMultiphysics.DOMAIN_SIZE), True)
        self.assertEqual(model_part.GetCommunicator().SynchronizeNonHistoricalVariable(KratosMultiphysics.TIME), True)
        self.assertEqual(model_part.GetCommunicator().SynchronizeNonHistoricalVariable(KratosMultiphysics.VECTOR_LAGRANGE_MULTIPLIER), True)
        self.assertEqual(model_part.GetCommunicator().SynchronizeNonHistoricalVariable(KratosMultiphysics.RECOVERED_STRESS), True)
        self.assertEqual(model_part.GetCommunicator().SynchronizeNonHistoricalVariable(KratosMultiphysics.GREEN_LAGRANGE_STRAIN_TENSOR), True)

    def testCommunicatorAssembleCurrentData(self):
        model_part = self.model.CreateModelPart("Test_model_part",1)
        self.assertEqual(model_part.GetCommunicator().AssembleCurrentData(KratosMultiphysics.IS_RESTARTED), True)
        self.assertEqual(model_part.GetCommunicator().AssembleCurrentData(KratosMultiphysics.DOMAIN_SIZE), True)
        self.assertEqual(model_part.GetCommunicator().AssembleCurrentData(KratosMultiphysics.TIME), True)
        self.assertEqual(model_part.GetCommunicator().AssembleCurrentData(KratosMultiphysics.VECTOR_LAGRANGE_MULTIPLIER), True)
        self.assertEqual(model_part.GetCommunicator().AssembleCurrentData(KratosMultiphysics.RECOVERED_STRESS), True)
        self.assertEqual(model_part.GetCommunicator().AssembleCurrentData(KratosMultiphysics.GREEN_LAGRANGE_STRAIN_TENSOR), True)

    def testCommunicatorAssembleNonHistoricalData(self):
        model_part = self.model.CreateModelPart("Test_model_part",1)
        self.assertEqual(model_part.GetCommunicator().AssembleNonHistoricalData(KratosMultiphysics.IS_RESTARTED), True)
        self.assertEqual(model_part.GetCommunicator().AssembleNonHistoricalData(KratosMultiphysics.DOMAIN_SIZE), True)
        self.assertEqual(model_part.GetCommunicator().AssembleNonHistoricalData(KratosMultiphysics.TIME), True)
        self.assertEqual(model_part.GetCommunicator().AssembleNonHistoricalData(KratosMultiphysics.VECTOR_LAGRANGE_MULTIPLIER), True)
        self.assertEqual(model_part.GetCommunicator().AssembleNonHistoricalData(KratosMultiphysics.RECOVERED_STRESS), True)
        self.assertEqual(model_part.GetCommunicator().AssembleNonHistoricalData(KratosMultiphysics.GREEN_LAGRANGE_STRAIN_TENSOR), True)

if __name__ == "__main__":
    KratosUnittest.main()
