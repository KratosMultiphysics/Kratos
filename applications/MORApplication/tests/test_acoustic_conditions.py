from __future__ import print_function, absolute_import, division
import KratosMultiphysics

import KratosMultiphysics.MORApplication as MOR
import KratosMultiphysics.KratosUnittest as KratosUnittest
import math

class TestAcousticConditions(KratosUnittest.TestCase):

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

    def test_acoustic_load_condition_2d2n(self):
        # create model
        model = KratosMultiphysics.Model()
        mp = model.CreateModelPart("domain")
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)
        mp.AddNodalSolutionStepVariable(MOR.ACOUSTIC_LOAD)

        # create nodes
        n1 = mp.CreateNewNode(1,0.0,0.0,0.0)
        mp.CreateNewNode(2,1.0,1.0,0.0)
        mp.GetProperties()[1].SetValue(MOR.ACOUSTIC_LOAD, 2.41e-3)

        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.PRESSURE, mp)
        mp.ProcessInfo[MOR.FREQUENCY] = 12
        cond = mp.CreateNewCondition("AcousticLoadCondition2D2N", 1, [1,2], mp.GetProperties()[1])

        # compute condition
        lhs = KratosMultiphysics.Matrix(0,0)
        rhs = KratosMultiphysics.Vector(0)
        cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)

        # assert
        self.assertAlmostEqual(rhs[0], rhs[1])
        self.assertAlmostEqual(rhs[0], 0.245394337342979)

        # apply different load on nodes
        n1.SetSolutionStepValue(MOR.ACOUSTIC_LOAD, 1e-2)
        rhs = KratosMultiphysics.Vector(0)
        cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)

        self.assertAlmostEqual(rhs[0], 1.263628102251608)
        self.assertAlmostEqual(rhs[1], 0.245394337342979)

    def test_acoustic_load_condition_2d1n(self):
        # parameters
        load = 2.41e-3
        freq = 12

        # create model
        model = KratosMultiphysics.Model()
        mp = model.CreateModelPart("domain")
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)
        mp.AddNodalSolutionStepVariable(MOR.ACOUSTIC_LOAD)

        # create nodes
        n1 = mp.CreateNewNode(1,0.0,0.0,0.0)
        mp.GetProperties()[1].SetValue(MOR.ACOUSTIC_LOAD, load)

        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.PRESSURE, mp)
        mp.ProcessInfo[MOR.FREQUENCY] = freq
        cond = mp.CreateNewCondition("AcousticLoadCondition2D1N", 1, [1], mp.GetProperties()[1])

        # compute condition
        lhs = KratosMultiphysics.Matrix(0,0)
        rhs = KratosMultiphysics.Vector(0)
        cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)

        # assert
        self.assertAlmostEqual(rhs[0], load*freq**2)

if __name__ == '__main__':
    KratosUnittest.main()
