from __future__ import print_function, absolute_import, division

import KratosMultiphysics
import KratosMultiphysics.mpi as MPI #TODO: do not import the so directly (but I need a nice Python module first)

import KratosMultiphysics.KratosUnittest as KratosUnittest

class TestMPICommunicatorSetUp(KratosUnittest.TestCase):

    def setUp(self):
        self.model = KratosMultiphysics.Model()

    def testMPICommunicatorSetUp(self):
        model_part = self.model.CreateModelPart("Test_model_part",1)

        # I'd do more complex tests, but this one should work in serial too (JC)
        self.assertNotRegex(model_part.GetCommunicator().__str__(), "MPICommunicator")

        MPI.ModelPartCommunicatorUtilities.SetMPICommunicator(model_part)

        self.assertRegex(model_part.GetCommunicator().__str__(), "MPICommunicator")

    def testMPICommunicatorMyPID(self):
        model_part = self.model.CreateModelPart("Test_model_part",1)
        MPI.ModelPartCommunicatorUtilities.SetMPICommunicator(model_part)
        #self.assertEqual(model_part.GetCommunicator().MyPID(), 0)

    def testMPICommunicatorTotalProcesses(self):
        model_part = self.model.CreateModelPart("Test_model_part",1)
        MPI.ModelPartCommunicatorUtilities.SetMPICommunicator(model_part)
        #self.assertEqual(model_part.GetCommunicator().TotalProcesses(), 1)

    def testMPICommunicatorGetNumberOfColors(self):
        model_part = self.model.CreateModelPart("Test_model_part",1)
        MPI.ModelPartCommunicatorUtilities.SetMPICommunicator(model_part)
        #self.assertEqual(model_part.GetCommunicator().GetNumberOfColors(), 1)

    def testMPICommunicatorSumAll(self):
        model_part = self.model.CreateModelPart("Test_model_part",1)
        MPI.ModelPartCommunicatorUtilities.SetMPICommunicator(model_part)
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

    def testMPICommunicatorMinAll(self):
        model_part = self.model.CreateModelPart("Test_model_part",1)
        MPI.ModelPartCommunicatorUtilities.SetMPICommunicator(model_part)
        a = 0
        self.assertEqual(model_part.GetCommunicator().MinAll(a), a)
        b = 0.0
        self.assertEqual(model_part.GetCommunicator().MinAll(b), b)

    def testMPICommunicatorMaxAll(self):
        model_part = self.model.CreateModelPart("Test_model_part",1)
        MPI.ModelPartCommunicatorUtilities.SetMPICommunicator(model_part)
        a = 0
        self.assertEqual(model_part.GetCommunicator().MaxAll(a), a)
        b = 0.0
        self.assertEqual(model_part.GetCommunicator().MaxAll(b), b)

    def testMPICommunicatorSynchronizeVariable(self):
        model_part = self.model.CreateModelPart("Test_model_part",1)
        MPI.ModelPartCommunicatorUtilities.SetMPICommunicator(model_part)
        self.assertEqual(model_part.GetCommunicator().SynchronizeVariable(KratosMultiphysics.IS_RESTARTED), True)
        self.assertEqual(model_part.GetCommunicator().SynchronizeVariable(KratosMultiphysics.DOMAIN_SIZE), True)
        self.assertEqual(model_part.GetCommunicator().SynchronizeVariable(KratosMultiphysics.TIME), True)
        self.assertEqual(model_part.GetCommunicator().SynchronizeVariable(KratosMultiphysics.VECTOR_LAGRANGE_MULTIPLIER), True)
        self.assertEqual(model_part.GetCommunicator().SynchronizeVariable(KratosMultiphysics.RECOVERED_STRESS), True)
        self.assertEqual(model_part.GetCommunicator().SynchronizeVariable(KratosMultiphysics.GREEN_LAGRANGE_STRAIN_TENSOR), True)

    def testMPICommunicatorSynchronizeNonHistoricalVariable(self):
        model_part = self.model.CreateModelPart("Test_model_part",1)
        MPI.ModelPartCommunicatorUtilities.SetMPICommunicator(model_part)
        self.assertEqual(model_part.GetCommunicator().SynchronizeNonHistoricalVariable(KratosMultiphysics.IS_RESTARTED), True)
        self.assertEqual(model_part.GetCommunicator().SynchronizeNonHistoricalVariable(KratosMultiphysics.DOMAIN_SIZE), True)
        self.assertEqual(model_part.GetCommunicator().SynchronizeNonHistoricalVariable(KratosMultiphysics.TIME), True)
        self.assertEqual(model_part.GetCommunicator().SynchronizeNonHistoricalVariable(KratosMultiphysics.VECTOR_LAGRANGE_MULTIPLIER), True)
        self.assertEqual(model_part.GetCommunicator().SynchronizeNonHistoricalVariable(KratosMultiphysics.RECOVERED_STRESS), True)
        self.assertEqual(model_part.GetCommunicator().SynchronizeNonHistoricalVariable(KratosMultiphysics.GREEN_LAGRANGE_STRAIN_TENSOR), True)

    def testMPICommunicatorAssembleCurrentData(self):
        model_part = self.model.CreateModelPart("Test_model_part",1)
        MPI.ModelPartCommunicatorUtilities.SetMPICommunicator(model_part)
        self.assertEqual(model_part.GetCommunicator().AssembleCurrentData(KratosMultiphysics.IS_RESTARTED), True)
        self.assertEqual(model_part.GetCommunicator().AssembleCurrentData(KratosMultiphysics.DOMAIN_SIZE), True)
        self.assertEqual(model_part.GetCommunicator().AssembleCurrentData(KratosMultiphysics.TIME), True)
        self.assertEqual(model_part.GetCommunicator().AssembleCurrentData(KratosMultiphysics.VECTOR_LAGRANGE_MULTIPLIER), True)
        self.assertEqual(model_part.GetCommunicator().AssembleCurrentData(KratosMultiphysics.RECOVERED_STRESS), True)
        self.assertEqual(model_part.GetCommunicator().AssembleCurrentData(KratosMultiphysics.GREEN_LAGRANGE_STRAIN_TENSOR), True)

    def testMPICommunicatorAssembleNonHistoricalData(self):
        model_part = self.model.CreateModelPart("Test_model_part",1)
        MPI.ModelPartCommunicatorUtilities.SetMPICommunicator(model_part)
        self.assertEqual(model_part.GetCommunicator().AssembleNonHistoricalData(KratosMultiphysics.IS_RESTARTED), True)
        self.assertEqual(model_part.GetCommunicator().AssembleNonHistoricalData(KratosMultiphysics.DOMAIN_SIZE), True)
        self.assertEqual(model_part.GetCommunicator().AssembleNonHistoricalData(KratosMultiphysics.TIME), True)
        self.assertEqual(model_part.GetCommunicator().AssembleNonHistoricalData(KratosMultiphysics.VECTOR_LAGRANGE_MULTIPLIER), True)
        self.assertEqual(model_part.GetCommunicator().AssembleNonHistoricalData(KratosMultiphysics.RECOVERED_STRESS), True)
        self.assertEqual(model_part.GetCommunicator().AssembleNonHistoricalData(KratosMultiphysics.GREEN_LAGRANGE_STRAIN_TENSOR), True)

if __name__ == "__main__":
    KratosUnittest.main()
