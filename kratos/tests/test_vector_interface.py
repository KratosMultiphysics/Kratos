from __future__ import print_function, absolute_import, division

import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics import *
import math


class TestVectorInterface(KratosUnittest.TestCase):

    def test_range(self):
        a = Vector(4)

        for i in range(len(a)):
            a[i] = i

        for i in range(4):
            self.assertEqual(a[i],i)

    def test_sum(self):
        a = Vector(3)
        a.fill(1.0)
        b = Vector(3)
        b.fill(2.0)

        c = a+b
        for it in c:
            self.assertEqual(it,3.0)

        d = b - a;
        for it in d:
            self.assertEqual(it,1.0)

        b -= a
        for i in range(len(b)):
            self.assertEqual(b[i], d[i])

    def test_list_construction(self):
        a = Vector([1,2,3])
        self.assertEqual(a[0],1)
        self.assertEqual(a[1],2)
        self.assertEqual(a[2],3)

    def test_scalar_op(self):
        e = Vector(2)
        e[0] = 0.0
        e[1] = 0.0

        e += 2

        for it in e:
            self.assertEqual(it,2.0)

        e -= 1.0
        for it in e:
            self.assertEqual(it,1.0 )

        e *= 5.0
        for it in e:
            self.assertEqual(it,5.0 )

        e /= 5.0
        for it in e:
            self.assertEqual(it,1.0 )

    def test_slicing(self):
        v = Vector(5)
        for i in range(len(v)):
            v[i] = i

        slice1 = v[1:4] # simple slice get
        self.assertEqual(len(slice1),3)
        for i,value in enumerate(slice1):
            self.assertEqual(value, i+1 )

        slice2 = v[2:] # open slice get
        self.assertEqual(len(slice2),3)
        for i,value in enumerate(slice2):
            self.assertEqual(value, i+2 )

        slice3 = v[:3] # open slice get
        self.assertEqual(len(slice3),3)
        for i,value in enumerate(slice3):
            self.assertEqual(value, i )

        # slice set
        v[0:3] = [2*i for i in range(3)]
        expected = [0,2,4,3,4]
        for i in range(5):
            self.assertEqual(v[i],expected[i])

        # modifying the vector through the slice
        # (numpy-like behaviour)
        s = v[0:3]

        # double assignment
        s[0] = 99
        self.assertEqual(v[0],99)

        # vector assignment
        w = Vector(3)
        w[:] = [10.,20.,30.]
        s[:] = w

        expected = [10,20,30,3,4]
        for i in range(5):
            self.assertEqual(v[i],expected[i])

        # slice assignment
        w = Vector(5)
        w[:] = [-1.,-2.,-3.,-4.,-5.]
        s[:] = w[1:4]

        expected = [-2.,-3.,-4.,3,4]
        for i in range(5):
            self.assertEqual(v[i],expected[i])

    def test_array_construction(self):
        a = Array3(5.0)
        for it in a:
            self.assertEqual(it,5.0)

        b = Vector(a)
        for it in b:
            self.assertEqual(it,5.0)

    def test_iter_read(self):
        a = Vector(4)

        counter = 0
        for i in range(len(a)):
            a[i] = i

        for it in a:
            self.assertEqual(it,counter)
            counter += 1

    #it is not possible to fix the next test. = is always type assignement in python
    #@KratosUnittest.expectedFailure
    #def test_iter_write(self):
        #a = Vector(4)

        #for it in a:
            #it = 2.0 #this simply does not work in python

        #for it in a:
            #self.assertEqual(it,2.0)



if __name__ == '__main__':
    KratosUnittest.main()
