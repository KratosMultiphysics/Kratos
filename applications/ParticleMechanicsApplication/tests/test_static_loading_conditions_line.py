import KratosMultiphysics

import KratosMultiphysics.MPMApplication as KratosMPM
import KratosMultiphysics.KratosUnittest as KratosUnittest
from math import sqrt


class TestStaticLoadingConditionsLine(KratosUnittest.TestCase):

    def test_MPMGridLineLoadCondition2D2N(self):
        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("solid_part")
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.POSITIVE_FACE_PRESSURE)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.NEGATIVE_FACE_PRESSURE)

        # Create nodes
        mp.CreateNewNode(1,0.0,0.0,0.0)
        mp.CreateNewNode(2,sqrt(2),sqrt(2.0),0.0)
        length = 2.0

        # Ensure that the property 1 is created
        mp.GetProperties()[1]

        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z,mp)

        cond = mp.CreateNewCondition("MPMGridLineLoadCondition2D2N", 1, [1,2], mp.GetProperties()[1])

        lhs = KratosMultiphysics.Matrix(0,0)
        rhs = KratosMultiphysics.Vector(0)

        # First we apply a constant LINE_LOAD to theh condition
        load_on_cond = KratosMultiphysics.Vector(3)
        load_on_cond[0] = 1.0
        load_on_cond[1] = 2.0
        load_on_cond[2] = 0.0 # Note that this is a 2D condition
        cond.SetValue(KratosMPM.LINE_LOAD,load_on_cond)
        cond.CalculateLocalSystem(lhs,rhs,mp.ProcessInfo)
        self.assertEqual(rhs[0],0.5*length)
        self.assertEqual(rhs[1],1.0*length)
        self.assertEqual(rhs[2],0.5*length)
        self.assertEqual(rhs[3],1.0*length)

        cond.SetValue(KratosMultiphysics.POSITIVE_FACE_PRESSURE,10.0)
        cond.SetValue(KratosMultiphysics.NEGATIVE_FACE_PRESSURE,5.0)
        cond.CalculateLocalSystem(lhs,rhs,mp.ProcessInfo)
        self.assertAlmostEqual(rhs[0],0.5*(1+5*sqrt(2.0)/2.0)*length)
        self.assertAlmostEqual(rhs[1],0.5*(2.0-5*sqrt(2.0)/2.0)*length)
        self.assertAlmostEqual(rhs[2],0.5*(1+5*sqrt(2.0)/2.0)*length)
        self.assertAlmostEqual(rhs[3], 0.5*(2.0-5*sqrt(2.0)/2.0)*length)

        ## Finally we apply TO THE NODES, a linearly varying POSITIVE_FACE_PRESSURE ranging from -100.0 to -200.0
        mp.Nodes[1].SetSolutionStepValue(KratosMultiphysics.POSITIVE_FACE_PRESSURE,0,-100.0)
        mp.Nodes[2].SetSolutionStepValue(KratosMultiphysics.POSITIVE_FACE_PRESSURE,0,-200.0)
        cond.CalculateLocalSystem(lhs,rhs,mp.ProcessInfo)
        reference_res = [-89.7453702522736,92.7453702522736,-113.31559629182519,116.3155962918252]
        self.assertAlmostEqual(rhs[0],reference_res[0])
        self.assertAlmostEqual(rhs[1],reference_res[1])
        self.assertAlmostEqual(rhs[2],reference_res[2])
        self.assertAlmostEqual(rhs[3],reference_res[3])

    def test_MPMGridLineLoadCondition2D2NAngle(self):
        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("solid_part")
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.POSITIVE_FACE_PRESSURE)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.NEGATIVE_FACE_PRESSURE)

        # Create nodes
        mp.CreateNewNode(1,0.0,0.0,0.0)
        mp.CreateNewNode(2,1.0,0.0,0.0)
        mp.CreateNewNode(3,1.0,1.0,0.0)
        length = 1.0

        # Ensure that the property 1 is created
        mp.GetProperties()[1]

        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z,mp)

        cond1 = mp.CreateNewCondition("MPMGridLineLoadCondition2D2N", 1, [1,2], mp.GetProperties()[1])
        cond2 = mp.CreateNewCondition("MPMGridLineLoadCondition2D2N", 2, [2,3], mp.GetProperties()[1])

        rhs = KratosMultiphysics.Vector(6)
        rhs[0] = 0.0
        rhs[1] = 0.0
        rhs[2] = 0.0
        rhs[3] = 0.0
        rhs[4] = 0.0
        rhs[5] = 0.0

        # First we apply a constant LINE_LOAD to then condition
        load_on_cond = KratosMultiphysics.Vector(3)
        load_on_cond[0] = 1.0
        load_on_cond[1] = 0.0
        load_on_cond[2] = 0.0 # note that this is a 2D condition
        cond1.SetValue(KratosMPM.LINE_LOAD,load_on_cond)
        load_on_cond[0] =  0.0
        load_on_cond[1] = -1.0
        cond2.SetValue(KratosMPM.LINE_LOAD,load_on_cond)

        linear_solver = KratosMultiphysics.SkylineLUFactorizationSolver()
        builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(linear_solver)
        scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()

        builder_and_solver.SetEchoLevel(0)

        builder_and_solver.Check(mp)
        builder_and_solver.SetUpDofSet(scheme, mp)
        builder_and_solver.SetUpSystem(mp)
        builder_and_solver.BuildRHS(scheme, mp, rhs)

        self.assertEqual(rhs[0], 0.5*length)
        self.assertEqual(rhs[1], 0.0*length)
        self.assertEqual(rhs[2], 0.5*length)
        self.assertEqual(rhs[3],-0.5*length)
        self.assertEqual(rhs[4], 0.0*length)
        self.assertEqual(rhs[5],-0.5*length)

if __name__ == '__main__':
    KratosUnittest.main()
