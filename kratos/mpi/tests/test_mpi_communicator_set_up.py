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

if __name__ == "__main__":
    KratosUnittest.main()
