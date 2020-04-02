from __future__ import print_function, absolute_import, division
import KratosMultiphysics

import KratosMultiphysics.MORApplication as MOR
import KratosMultiphysics.KratosUnittest as KratosUnittest
import math

class TestAcousticRobinConditionLine(KratosUnittest.TestCase):

    def test_acoustic_robin_condition_2d2n(self):
        # create model
        model = KratosMultiphysics.Model()
        mp = model.CreateModelPart("domain")
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)

        # create nodes
        mp.CreateNewNode(1,0.0,0.0,0.0)
        mp.CreateNewNode(2,1.0,1.0,0.0)
        mp.GetProperties()[1].SetValue(MOR.ADMITTANCE, 2.41e-3)

        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.PRESSURE, mp)

        cond = mp.CreateNewCondition("AcousticRobinCondition2D2N", 1, [1,2], mp.GetProperties()[1])

        # compute condition
        lhs = KratosMultiphysics.Matrix(0,0)
        cond.CalculateDampingMatrix(lhs,mp.ProcessInfo)

        # assert
        self.assertAlmostEqual(lhs[0,0], 0.001136084895106)
        self.assertAlmostEqual(lhs[1,1], 0.001136084895106)
        self.assertAlmostEqual(lhs[0,1], 5.680424475531930e-04)
        self.assertAlmostEqual(lhs[1,0], 5.680424475531930e-04)


if __name__ == '__main__':
    KratosUnittest.main()
