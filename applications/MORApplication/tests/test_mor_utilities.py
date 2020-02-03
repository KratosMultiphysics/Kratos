
from __future__ import print_function, absolute_import, division
import KratosMultiphysics

import KratosMultiphysics.MORApplication as KratosMOR
import KratosMultiphysics.KratosUnittest as KratosUnittest

from math import sqrt

class TestMORUtilities(KratosUnittest.TestCase):
    def test_pair_complex_conjugates(self):
        # mixed complex and real
        v = KratosMultiphysics.ComplexVector([1.-4.j, 1.+4.j, -10.1, -1.+2.j, 1.+9.j, 1.-9.j, 12., -1.-2.j, 3.+2.j, 3.-2.j, 1.-4.j, 1.+4.j])
        KratosMOR.PairComplexConjugates(v)

        v_exp = KratosMultiphysics.ComplexVector([-1.-2.j, -1.+2.j, 1.-4.j, 1.+4.j, 1.-4.j, 1.+4.j, 1.-9.j, 1.+9.j, 3.-2.j, 3.+2.j, -10.1, 12])
        for i in range(11):
            self.assertAlmostEqual(v[i], v_exp[i], 7)

        # only imaginary
        v = KratosMultiphysics.ComplexVector([-5.j, 5.j, 8.j, -8.j])
        KratosMOR.PairComplexConjugates(v)

        v_exp = KratosMultiphysics.ComplexVector([-5.j, 5.j, -8.j, 8.j])
        for i in range(3):
            self.assertAlmostEqual(v[i], v_exp[i], 7)

        # only real
        v = KratosMultiphysics.ComplexVector([-5., 5., 8., -8.])
        KratosMOR.PairComplexConjugates(v)

        v_exp = KratosMultiphysics.ComplexVector([-8., -5., 5., 8.])
        for i in range(3):
            self.assertAlmostEqual(v[i], v_exp[i], 7)

    def test_generalized_eigenvalue_utility(self):
        A0 = KratosMultiphysics.Matrix(4,4,0)
        A0[0,0] = -7
        A0[0,1] = 2
        A0[1,0] = 2
        A0[2,0] = 4
        A0[0,2] = 4
        A0[1,1] = -4
        A0[1,2] = 2
        A0[2,1] = 2
        A0[2,2] = -9
        A0[2,3] = 3
        A0[3,2] = 3
        A0[3,3] = -3
        A1 = KratosMultiphysics.Matrix(4,4,0)
        A1[0,0] = .4
        A1[0,2] = -.3
        A1[2,0] = -.3
        A1[2,2] = .5
        A1[2,3] = -.2
        A1[3,2] = -.2
        A1[3,3] = .2
        A2 = KratosMultiphysics.Matrix(4,4,0)
        A2[0,0] = 3
        A2[1,1] = 1
        A2[2,2] = 3
        A2[3,3] = 1
        e = KratosMultiphysics.ComplexVector(8)

        KratosMOR.ComputePolynomialEigenvalues(A0,A1,A2,e)

        e_exp = KratosMultiphysics.ComplexVector([-2.44984944370563, -2.15361619803731, -1.62477834052925, 1.47524114347567, \
            0.335294429778544, 2.03635097664370, 2.22790873204791, -0.346551299673631])
        for i in range(7):
            self.assertAlmostEqual(e[i], e_exp[i], 7)

    def test_complex_generalized_eigenvalue_utility(self):
        A = KratosMultiphysics.ComplexMatrix(2,2,0)
        A[0,0] = 1j/sqrt(2)
        A[1,1] = 1j
        B = KratosMultiphysics.ComplexMatrix(2,2,0)
        B[0,1] = 1
        B[1,0] = -1/sqrt(2)
        e = KratosMultiphysics.ComplexVector(2)
        
        KratosMOR.ComputeGeneralizedEigenvalues(A,B,e)

        e_exp = KratosMultiphysics.ComplexVector([1, -1])
        for i in range(len(e)):
            self.assertAlmostEqual(e[i], e_exp[i], 7)

if __name__ == '__main__':
    KratosUnittest.main()