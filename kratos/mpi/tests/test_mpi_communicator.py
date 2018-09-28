import KratosMultiphysics as Kratos
import KratosMultiphysics.mpi as MPI #TODO: do not import the so directly (but I need a nice Python module first)

import KratosMultiphysics.KratosUnittest as UnitTest

class MPICommunicatorTest(UnitTest.TestCase):

    def testMPICommunicatorSetup(self):
        model_part = Kratos.ModelPart("Test_model_part")

        # I'd do more complex tests, but this one should work in serial too (JC)
        self.assertNotEqual(model_part.GetCommunicator().__str__(), "MPICommunicator")

        MPI.ModelPartCommunicatorUtilities.SetMPICommunicator(model_part)

        self.assertEqual(model_part.GetCommunicator().__str__(), "MPICommunicator")

if __name__ == "__main__":
    UnitTest.main()
