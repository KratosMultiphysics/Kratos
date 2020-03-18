from __future__ import print_function, absolute_import, division
import KratosMultiphysics

import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest
import math

class TestLoadingConditionsLine(KratosUnittest.TestCase):

    def _LineLoadCondition3D2NRotDof(self, prefix = ""):
        current_model = KratosMultiphysics.Model()
        dim = 3
        mp = current_model.CreateModelPart("solid_part")
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.ROTATION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION_MOMENT)

        #create nodes
        mp.CreateNewNode(1,0.0,0.0,0.0)
        mp.CreateNewNode(2,1.0,1.0,0.0)
        length = math.sqrt(2)

        #ensure that the property 1 is created
        mp.GetProperties()[1]

        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z,mp)

        # NOTE: If adding rotations on the future
        #KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_X, KratosMultiphysics.REACTION_MOMENT_X,mp)
        #KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_Y, KratosMultiphysics.REACTION_MOMENT_Y,mp)
        #KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_Z, KratosMultiphysics.REACTION_MOMENT_Z,mp)

        cond = mp.CreateNewCondition(prefix + "LineLoadCondition3D2N", 1, [1,2], mp.GetProperties()[1])

        cond.SetValue(KratosMultiphysics.LOCAL_AXIS_2, [-1.0, 1.0, 0.0])

        lhs = KratosMultiphysics.Matrix(0,0)
        rhs = KratosMultiphysics.Vector(0)

        #first we apply a constant LINE_LOAD to theh condition
        Line_Load_i = 10000.00/math.sqrt(2) #apply a 45 degrees load

        load_on_cond = KratosMultiphysics.Vector(3)
        load_on_cond[0] = 0.00
        load_on_cond[1] = -Line_Load_i
        load_on_cond[2] = -Line_Load_i
        cond.SetValue(StructuralMechanicsApplication.LINE_LOAD,load_on_cond)
        cond.CalculateLocalSystem(lhs,rhs,mp.ProcessInfo)

        Nodal_Transversal_Forces = Line_Load_i*length/2.00

        self.assertEqual(rhs[0],0.00)
        self.assertEqual(rhs[1],-Nodal_Transversal_Forces)
        self.assertEqual(rhs[2],-Nodal_Transversal_Forces)
        self.assertEqual(rhs[3],0.00)
        self.assertEqual(rhs[4],-Nodal_Transversal_Forces)
        self.assertEqual(rhs[5],-Nodal_Transversal_Forces)

        # NOTE: If adding rotations on the future
        #Nodal_Moments = Line_Load_i*length*length/12.00/math.sqrt(2)

        #self.assertEqual(rhs[0],0.00)
        #self.assertEqual(rhs[1],-Nodal_Transversal_Forces)
        #self.assertEqual(rhs[2],-Nodal_Transversal_Forces)
        #self.assertAlmostEqual(rhs[3],-Nodal_Moments)
        #self.assertAlmostEqual(rhs[4],Nodal_Moments)
        #self.assertAlmostEqual(rhs[5],-Nodal_Moments)
        #self.assertEqual(rhs[6],0.00)
        #self.assertEqual(rhs[7],-Nodal_Transversal_Forces)
        #self.assertEqual(rhs[8],-Nodal_Transversal_Forces)
        #self.assertAlmostEqual(rhs[9],Nodal_Moments)
        #self.assertAlmostEqual(rhs[10],-Nodal_Moments)
        #self.assertAlmostEqual(rhs[11],Nodal_Moments)

    def _LineLoadCondition2D2N(self, prefix = ""):
        current_model = KratosMultiphysics.Model()
        dim = 2
        mp = current_model.CreateModelPart("solid_part")
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.POSITIVE_FACE_PRESSURE)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.NEGATIVE_FACE_PRESSURE)

        #create nodes
        mp.CreateNewNode(1,0.0,0.0,0.0)
        mp.CreateNewNode(2,math.sqrt(2),math.sqrt(2.0),0.0)
        lenght = 2.0

        #ensure that the property 1 is created
        mp.GetProperties()[1]

        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z,mp)

        cond = mp.CreateNewCondition(prefix + "LineLoadCondition2D2N", 1, [1,2], mp.GetProperties()[1])

        lhs = KratosMultiphysics.Matrix(0,0)
        rhs = KratosMultiphysics.Vector(0)

        #first we apply a constant LINE_LOAD to theh condition
        load_on_cond = KratosMultiphysics.Vector(3)
        load_on_cond[0] = 1.0
        load_on_cond[1] = 2.0
        load_on_cond[2] = 0.0 #note that this is a 2D condition
        cond.SetValue(StructuralMechanicsApplication.LINE_LOAD,load_on_cond)
        cond.CalculateLocalSystem(lhs,rhs,mp.ProcessInfo)
        self.assertEqual(rhs[0],0.5*lenght)
        self.assertEqual(rhs[1],1.0*lenght)
        self.assertEqual(rhs[2],0.5*lenght)
        self.assertEqual(rhs[3],1.0*lenght)

        ##now the condition is something like
        ##    2
        ##  +/
        ##  /-
        ## 1
        #with the + and - sides as indicated. (positive face is the one from which the normal to the face goes out)
        ##applying a NEGATIVE_FACE_PRESSURE implies applying a pressure on the - face, which hence goes in the direction of the normal (-1,1)
        ##applying a POSITIVE_FACE_PRESSURE implies applying a distributed force to the + face, which corresponds to a force going in the direciton of (1,-1), that is, the opposite to the normal
        ##here
        ##we add to it a NEGATIVE_FACE_PRESSURE of 10 corresponding to a a distributed force (10*cos(45),-10*sin(45))
        ##togheter with a POSITIVE_FACE_PRESSURE of 5 corresponding to a LINE_LOAD (-5*cos(45),5*cos(45))
        ##together with the previousy applied LINE_LOAD this gives an equivalent load of (1+5*cos(45),2.0-5*sin(45))
        cond.SetValue(KratosMultiphysics.POSITIVE_FACE_PRESSURE,10.0)
        cond.SetValue(KratosMultiphysics.NEGATIVE_FACE_PRESSURE,5.0)
        cond.CalculateLocalSystem(lhs,rhs,mp.ProcessInfo)
        self.assertAlmostEqual(rhs[0],0.5*(1+5*math.sqrt(2.0)/2.0)*lenght)
        self.assertAlmostEqual(rhs[1],0.5*(2.0-5*math.sqrt(2.0)/2.0)*lenght)
        self.assertAlmostEqual(rhs[2],0.5*(1+5*math.sqrt(2.0)/2.0)*lenght)
        self.assertAlmostEqual(rhs[3], 0.5*(2.0-5*math.sqrt(2.0)/2.0)*lenght)


        ##finally we apply TO THE NODES, a linearly varying POSITIVE_FACE_PRESSURE ranging from -100.0 to -200.0
        mp.Nodes[1].SetSolutionStepValue(KratosMultiphysics.POSITIVE_FACE_PRESSURE,0,-100.0)
        mp.Nodes[2].SetSolutionStepValue(KratosMultiphysics.POSITIVE_FACE_PRESSURE,0,-200.0)
        cond.CalculateLocalSystem(lhs,rhs,mp.ProcessInfo)
        reference_res = [-89.7453702522736,92.7453702522736,-113.31559629182519,116.3155962918252]
        self.assertAlmostEqual(rhs[0],reference_res[0])
        self.assertAlmostEqual(rhs[1],reference_res[1])
        self.assertAlmostEqual(rhs[2],reference_res[2])
        self.assertAlmostEqual(rhs[3],reference_res[3])

    def _LineLoadCondition2D3N(self, prefix = ""):
        current_model = KratosMultiphysics.Model()
        dim = 2
        mp = current_model.CreateModelPart("solid_part")
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.POSITIVE_FACE_PRESSURE)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.NEGATIVE_FACE_PRESSURE)

        #create nodes
        mp.CreateNewNode(1,0.0,0.0,0.0)
        mp.CreateNewNode(2,math.sqrt(2),math.sqrt(2.0),0.0)
        mp.CreateNewNode(3,0.5 * math.sqrt(2),0.5 * math.sqrt(2.0),0.0)
        lenght = 2.0

        #ensure that the property 1 is created
        mp.GetProperties()[1]

        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z,mp)

        cond = mp.CreateNewCondition(prefix + "LineLoadCondition2D3N", 1, [1,2,3], mp.GetProperties()[1])

        lhs = KratosMultiphysics.Matrix(0,0)
        rhs = KratosMultiphysics.Vector(0)

        #first we apply a constant LINE_LOAD to theh condition
        load_on_cond = KratosMultiphysics.Vector(3)
        load_on_cond[0] = 1.0
        load_on_cond[1] = 2.0
        load_on_cond[2] = 0.0 #note that this is a 2D condition
        cond.SetValue(StructuralMechanicsApplication.LINE_LOAD,load_on_cond)
        cond.CalculateLocalSystem(lhs,rhs,mp.ProcessInfo)

        self.assertAlmostEqual(rhs[0],1.0/6.0*lenght)
        self.assertAlmostEqual(rhs[1],1.0/3.0*lenght)
        self.assertAlmostEqual(rhs[2],1.0/6.0*lenght)
        self.assertAlmostEqual(rhs[3],1.0/3.0*lenght)
        self.assertAlmostEqual(rhs[4],2.0/3.0*lenght)
        self.assertAlmostEqual(rhs[5],4.0/3.0*lenght)

        ##now the condition is something like
        ##    2
        ##  +/
        ##  /-
        ## 1
        #with the + and - sides as indicated. (positive face is the one from which the normal to the face goes out)
        ##applying a NEGATIVE_FACE_PRESSURE implies applying a pressure on the - face, which hence goes in the direction of the normal (-1,1)
        ##applying a POSITIVE_FACE_PRESSURE implies applying a distributed force to the + face, which corresponds to a force going in the direciton of (1,-1), that is, the opposite to the normal
        ##here
        ##we add to it a NEGATIVE_FACE_PRESSURE of 10 corresponding to a a distributed force (10*cos(45),-10*sin(45))
        ##togheter with a POSITIVE_FACE_PRESSURE of 5 corresponding to a LINE_LOAD (-5*cos(45),5*cos(45))
        ##together with the previousy applied LINE_LOAD this gives an equivalent load of (1+5*cos(45),2.0-5*sin(45))
        cond.SetValue(KratosMultiphysics.POSITIVE_FACE_PRESSURE,10.0)
        cond.SetValue(KratosMultiphysics.NEGATIVE_FACE_PRESSURE,5.0)
        cond.CalculateLocalSystem(lhs,rhs,mp.ProcessInfo)

        self.assertAlmostEqual(rhs[0],1.0/6.0*(1+5*math.sqrt(2.0)/2.0)*lenght)
        self.assertAlmostEqual(rhs[1],1.0/6.0*(2.0-5*math.sqrt(2.0)/2.0)*lenght)
        self.assertAlmostEqual(rhs[2],1.0/6.0*(1+5*math.sqrt(2.0)/2.0)*lenght)
        self.assertAlmostEqual(rhs[3],1.0/6.0*(2.0-5*math.sqrt(2.0)/2.0)*lenght)
        self.assertAlmostEqual(rhs[4],2.0/3.0*(1+5*math.sqrt(2.0)/2.0)*lenght)
        self.assertAlmostEqual(rhs[5],2.0/3.0*(2.0-5*math.sqrt(2.0)/2.0)*lenght)

        ## Finally we apply TO THE NODES, a linearly varying POSITIVE_FACE_PRESSURE ranging from -100.0 to -200.0
        mp.Nodes[1].SetSolutionStepValue(KratosMultiphysics.POSITIVE_FACE_PRESSURE,0,-100.0)
        mp.Nodes[2].SetSolutionStepValue(KratosMultiphysics.POSITIVE_FACE_PRESSURE,0,-200.0)
        mp.Nodes[3].SetSolutionStepValue(KratosMultiphysics.POSITIVE_FACE_PRESSURE,0,-150.0)
        cond.CalculateLocalSystem(lhs,rhs,mp.ProcessInfo)
        reference_res = [-22.0584,23.0584,-45.6286,46.6286,-135.374,139.374]

        self.assertAlmostEqual(rhs[0],reference_res[0], 4)
        self.assertAlmostEqual(rhs[1],reference_res[1], 4)
        self.assertAlmostEqual(rhs[2],reference_res[2], 4)
        self.assertAlmostEqual(rhs[3],reference_res[3], 4)
        self.assertAlmostEqual(rhs[4],reference_res[4], 4)
        self.assertAlmostEqual(rhs[5],reference_res[5], 4)

    def _LineLoadCondition2D2NAngle(self, prefix = ""):
        current_model = KratosMultiphysics.Model()
        dim = 2
        mp = current_model.CreateModelPart("solid_part")
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.POSITIVE_FACE_PRESSURE)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.NEGATIVE_FACE_PRESSURE)

        #create nodes
        mp.CreateNewNode(1,0.0,0.0,0.0)
        mp.CreateNewNode(2,1.0,0.0,0.0)
        mp.CreateNewNode(3,1.0,1.0,0.0)
        lenght = 1.0

        #ensure that the property 1 is created
        mp.GetProperties()[1]

        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z,mp)

        cond1 = mp.CreateNewCondition(prefix + "LineLoadCondition2D2N", 1, [1,2], mp.GetProperties()[1])
        cond2 = mp.CreateNewCondition(prefix + "LineLoadCondition2D2N", 2, [2,3], mp.GetProperties()[1])

        rhs = KratosMultiphysics.Vector(6)
        rhs[0] = 0.0
        rhs[1] = 0.0
        rhs[2] = 0.0
        rhs[3] = 0.0
        rhs[4] = 0.0
        rhs[5] = 0.0

        # First we apply a constant LINE_LOAD to theh condition
        load_on_cond = KratosMultiphysics.Vector(3)
        load_on_cond[0] = 1.0
        load_on_cond[1] = 0.0
        load_on_cond[2] = 0.0 #note that this is a 2D condition
        cond1.SetValue(StructuralMechanicsApplication.LINE_LOAD,load_on_cond)
        load_on_cond[0] =  0.0
        load_on_cond[1] = -1.0
        cond2.SetValue(StructuralMechanicsApplication.LINE_LOAD,load_on_cond)

        linear_solver = KratosMultiphysics.SkylineLUFactorizationSolver()
        builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(linear_solver)
        scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()

        builder_and_solver.SetEchoLevel(0)

        builder_and_solver.Check(mp)
        builder_and_solver.SetUpDofSet(scheme, mp);
        builder_and_solver.SetUpSystem(mp);
        builder_and_solver.BuildRHS(scheme, mp, rhs)

        self.assertEqual(rhs[0], 0.5*lenght)
        self.assertEqual(rhs[1], 0.0*lenght)
        self.assertEqual(rhs[2], 0.5*lenght)
        self.assertEqual(rhs[3],-0.5*lenght)
        self.assertEqual(rhs[4], 0.0*lenght)
        self.assertEqual(rhs[5],-0.5*lenght)

    def test_SDLineLoadCondition3D2NRotDof(self):
        self._LineLoadCondition3D2NRotDof("SmallDisplacement")

    def test_LineLoadCondition3D2NRotDof(self):
        self._LineLoadCondition3D2NRotDof()

    def test_SDLineLoadCondition2D2N(self):
        self._LineLoadCondition2D2N("SmallDisplacement")

    def test_LineLoadCondition2D2N(self):
        self._LineLoadCondition2D2N()

    def test_SDLineLoadCondition2D3N(self):
        self._LineLoadCondition2D3N("SmallDisplacement")

    def test_LineLoadCondition2D3N(self):
        self._LineLoadCondition2D3N()

    def test_SDLineLoadCondition2D2NAngle(self):
        self._LineLoadCondition2D2NAngle("SmallDisplacement")

    def test_LineLoadCondition2D2NAngle(self):
        self._LineLoadCondition2D2NAngle()

if __name__ == '__main__':
    KratosUnittest.main()
