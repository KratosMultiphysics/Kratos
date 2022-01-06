import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics
import KratosMultiphysics.TrilinosApplication as KratosTrilinos

class TestTrilinosMatrix(KratosUnittest.TestCase):

    def test_resize(self):
        comm = KratosTrilinos.CreateEpetraCommunicator(KratosMultiphysics.DataCommunicator.GetDefault())
        space = KratosTrilinos.TrilinosSparseSpace()

        pb = space.CreateEmptyVectorPointer(comm)
        space.ResizeVector(pb,2)

        n = space.Size(pb.GetReference())
        self.assertEqual(n,2)

if __name__ == '__main__':
    KratosUnittest.main()