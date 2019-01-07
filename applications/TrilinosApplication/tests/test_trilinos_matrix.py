from __future__ import print_function, absolute_import, division

import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics
import KratosMultiphysics.mpi as KratosMPI
import KratosMultiphysics.TrilinosApplication as KratosTrilinos


class TestTrilinosMatrix(KratosUnittest.TestCase):

    def test_resize(self):
        comm = KratosMultiphysics.TrilinosApplication.CreateCommunicator()
        space = KratosMultiphysics.TrilinosApplication.TrilinosSparseSpace()
        
        pb = space.CreateEmptyVectorPointer(comm)
        space.ResizeVector(pb,2)

        n = space.Size(pb.GetReference())
        self.assertEqual(n,2)
        
if __name__ == '__main__':
    KratosUnittest.main()