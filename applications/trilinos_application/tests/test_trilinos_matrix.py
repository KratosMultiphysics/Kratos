from __future__ import print_function, absolute_import, division

import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics
import KratosMultiphysics.mpi as KratosMPI
import KratosMultiphysics.TrilinosApplication as KratosTrilinos


class TestMPICommunicator(KratosUnittest.TestCase):

    def test_resize(self):
        comm = KratosMultiphysics.TrilinosApplication.CreateCommunicator()
        space = KratosMultiphysics.TrilinosApplication.TrilinosSparseSpace()
        
        print(dir(space))
        print(space)
        pb = space.CreateEmptyVectorPointer(comm)
        print("000")
        space.ResizeVector(pb,2)
        print("111")
        rb = pb.GetReference()
        print("222")
        n = space.Size(rb)
        self.assertEqual(n,0)
        
if __name__ == '__main__':
    KratosUnittest.main()