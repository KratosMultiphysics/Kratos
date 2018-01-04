from __future__ import print_function, absolute_import, division

import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics import *
import math


class TestMatrixInterface(KratosUnittest.TestCase):
        
    def test_range(self):
        a = Matrix(2,3)
        
        self.assertEqual(2,a.Size1())
        self.assertEqual(3,a.Size2())
        
        for i in range(a.Size1()):
            for j in range(a.Size2()):
                a[i,j] = i
        
        print(a)
        
        
if __name__ == '__main__':
    KratosUnittest.main()
