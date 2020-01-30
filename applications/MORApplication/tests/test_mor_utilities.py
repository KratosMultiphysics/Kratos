
from __future__ import print_function, absolute_import, division
import KratosMultiphysics

import KratosMultiphysics.MORApplication as KratosMOR
import KratosMultiphysics.KratosUnittest as KratosUnittest

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

if __name__ == '__main__':
    KratosUnittest.main()