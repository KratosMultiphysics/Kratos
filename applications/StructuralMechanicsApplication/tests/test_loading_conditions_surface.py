from __future__ import print_function, absolute_import, division
import KratosMultiphysics

import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest
import math


class TestLoadingConditionsSurface(KratosUnittest.TestCase):
        
    def _SimplestSurfaceLoadCondition3D4N(self, prefix = ""):
        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("solid_part")
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)

        #create nodes
        mp.CreateNewNode(1,0.0,0.0,0.0)
        mp.CreateNewNode(2,1.0,0.0,0.0)
        mp.CreateNewNode(3,1.0,1.0,0.0)
        mp.CreateNewNode(4,0.0,1.0,0.0)

        #ensure that the property 1 is created
        mp.GetProperties()[1]

        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z,mp)

        cond = mp.CreateNewCondition(prefix + "SurfaceLoadCondition3D4N", 1, [1,2,3,4], mp.GetProperties()[1])

        #cond.Check()

        lhs = KratosMultiphysics.Matrix(0,0)
        rhs = KratosMultiphysics.Vector(0)

        #first we apply a constant SURFACE_LOAD to theh condition
        load_on_cond = KratosMultiphysics.Vector(3)
        load_on_cond[0] =  0.0
        load_on_cond[1] =  0.0
        load_on_cond[2] = -1.0
        cond.SetValue(StructuralMechanicsApplication.SURFACE_LOAD,load_on_cond)
        cond.CalculateLocalSystem(lhs,rhs,mp.ProcessInfo)

        reference_res = [0,0,-0.25,0,0,-0.25,0,0,-0.25,0,0,-0.25]

        for i in range(len(rhs)):
            self.assertAlmostEqual(rhs[i],reference_res[i])

    def _SimpleOrientationSurfaceLoadCondition3D4N(self, prefix = ""):
        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("solid_part")
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)

        #create nodes
        mp.CreateNewNode(1,4.0,6.0,4.7)
        mp.CreateNewNode(2,4.7,6.0,4.7)
        mp.CreateNewNode(3,4.0,6.0,4.0)
        mp.CreateNewNode(4,4.7,6.0,4.0)

        #ensure that the property 1 is created
        mp.GetProperties()[1]

        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z,mp)

        cond = mp.CreateNewCondition(prefix + "SurfaceLoadCondition3D4N", 1, [1,2,4,3], mp.GetProperties()[1])

        #cond.Check()

        lhs = KratosMultiphysics.Matrix(0,0)
        rhs = KratosMultiphysics.Vector(0)

        #first we apply a constant SURFACE_LOAD to theh condition
        load_on_cond = KratosMultiphysics.Vector(3)
        load_on_cond[0] =  0.0
        load_on_cond[1] = -1.0
        load_on_cond[2] =  0.0
        cond.SetValue(StructuralMechanicsApplication.SURFACE_LOAD,load_on_cond)
        cond.CalculateLocalSystem(lhs,rhs,mp.ProcessInfo)

        reference_res = [0,-0.1225,0,0,-0.1225,0,0,-0.1225,0,0,-0.1225,0]

        for i in range(len(rhs)):
            self.assertAlmostEqual(rhs[i],reference_res[i])

    def _SimpleSurfaceLoadCondition3D4N(self, prefix = ""):
        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("solid_part")
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)

        #create nodes
        mp.CreateNewNode(1,0.0,0.0,0.0)
        mp.CreateNewNode(2,2.0,0.0,0.0)
        mp.CreateNewNode(3,2.0,1.0,0.0)
        mp.CreateNewNode(4,0.0,1.0,0.0)

        #ensure that the property 1 is created
        mp.GetProperties()[1]

        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z,mp)

        cond = mp.CreateNewCondition(prefix + "SurfaceLoadCondition3D4N", 1, [1,2,3,4], mp.GetProperties()[1])

        #cond.Check()

        lhs = KratosMultiphysics.Matrix(0,0)
        rhs = KratosMultiphysics.Vector(0)

        #first we apply a constant SURFACE_LOAD to theh condition
        load_on_cond = KratosMultiphysics.Vector(3)
        load_on_cond[0] =  0.0
        load_on_cond[1] =  0.0
        load_on_cond[2] = -1.0
        cond.SetValue(StructuralMechanicsApplication.SURFACE_LOAD,load_on_cond)
        cond.CalculateLocalSystem(lhs,rhs,mp.ProcessInfo)

        reference_res = [0,0,-0.5,0,0,-0.5,0,0,-0.5,0,0,-0.5]

        for i in range(len(rhs)):
            self.assertAlmostEqual(rhs[i],reference_res[i])

    def _SurfaceLoadCondition3D4N(self, prefix = ""):
        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("solid_part")
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.POSITIVE_FACE_PRESSURE)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.NEGATIVE_FACE_PRESSURE)

        #create nodes
        mp.CreateNewNode(1,0.0,0.0,0.0)
        mp.CreateNewNode(2,0.0,0.0,1.0)
        mp.CreateNewNode(3,math.sqrt(2),math.sqrt(2.0),0.0)
        mp.CreateNewNode(4,math.sqrt(2),math.sqrt(2.0),1.0)
        lenght = 2.0
        width = 1.0

        #ensure that the property 1 is created
        mp.GetProperties()[1]

        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z,mp)

        cond = mp.CreateNewCondition(prefix + "SurfaceLoadCondition3D4N", 1, [1,2,4,3], mp.GetProperties()[1])

        #cond.Check()

        lhs = KratosMultiphysics.Matrix(0,0)
        rhs = KratosMultiphysics.Vector(0)

        #first we apply a constant SURFACE_LOAD to theh condition
        load_on_cond = KratosMultiphysics.Vector(3)
        load_on_cond[0] = 1.0
        load_on_cond[1] = 2.0
        load_on_cond[2] = 0.0
        cond.SetValue(StructuralMechanicsApplication.SURFACE_LOAD,load_on_cond)
        cond.CalculateLocalSystem(lhs,rhs,mp.ProcessInfo)
        self.assertAlmostEqual(rhs[0],0.25*lenght); self.assertAlmostEqual(rhs[1],0.5*1.0*lenght); self.assertAlmostEqual(rhs[2],0.0)
        self.assertAlmostEqual(rhs[3],0.25*lenght); self.assertAlmostEqual(rhs[4],0.5*1.0*lenght); self.assertAlmostEqual(rhs[5],0.0)
        self.assertAlmostEqual(rhs[6],0.25*lenght); self.assertAlmostEqual(rhs[7],0.5*1.0*lenght); self.assertAlmostEqual(rhs[8],0.0)
        self.assertAlmostEqual(rhs[9],0.25*lenght); self.assertAlmostEqual(rhs[10],0.5*1.0*lenght); self.assertAlmostEqual(rhs[11],0.0)

        ##now the condition is something like
        ##    2
        ##  +/
        ##  /-
        ## 1
        #with the + and - sides as indicated. (positive face is the one from which the normal to the face goes out)
        ##applying a NEGATIVE_FACE_PRESSURE implies applying a pressure on the - face, which hence goes in the direction of the normal (-1,1)
        ##applying a POSITIVE_FACE_PRESSURE implies applying a distributed force to the + face, which corresponds to a force going in the direciton of (1,-1), that is, the opposite to the normal
        ##here
        ##we add to it a NEGATIVE_FACE_PRESSURE of 10 corresponding to a a distributed force (10*cos(45),-10*sin(45))*width
        ##togheter with a POSITIVE_FACE_PRESSURE of 5 corresponding to a LINE_LOAD (-5*cos(45),5*cos(45))*width
        ##together with the previousy applied LINE_LOAD this gives an equivalent load of 0.5*(1+5*cos(45),2.0-5*sin(45))
        cond.SetValue(KratosMultiphysics.POSITIVE_FACE_PRESSURE,10.0)
        cond.SetValue(KratosMultiphysics.NEGATIVE_FACE_PRESSURE,5.0)
        cond.CalculateLocalSystem(lhs,rhs,mp.ProcessInfo)
        fxnodal = (1.0 + 5.0*math.sqrt(2.0)/2.0) * lenght*width / 4.0
        fynodal = (2.0 - 5.0*math.sqrt(2.0)/2.0) * lenght*width / 4.0
        for i in range(4):
            base = i*3
            self.assertAlmostEqual(rhs[base+0],fxnodal)
            self.assertAlmostEqual(rhs[base+1],fynodal)
            self.assertAlmostEqual(rhs[base+2],0.0)

        ##finally we apply TO THE NODES, a linearly varying POSITIVE_FACE_PRESSURE ranging from -100.0 to -200.0
        mp.Nodes[1].SetSolutionStepValue(KratosMultiphysics.POSITIVE_FACE_PRESSURE,0,-100.0)
        mp.Nodes[2].SetSolutionStepValue(KratosMultiphysics.POSITIVE_FACE_PRESSURE,0,-100.0)
        mp.Nodes[3].SetSolutionStepValue(KratosMultiphysics.POSITIVE_FACE_PRESSURE,0,-200.0)
        mp.Nodes[4].SetSolutionStepValue(KratosMultiphysics.POSITIVE_FACE_PRESSURE,0,-200.0)
        cond.CalculateLocalSystem(lhs,rhs,mp.ProcessInfo)
        reference_res = [-44.872685126136794,46.372685126136815,0.0,-44.87268512613681,46.3726851261368,0.0,-56.657798145912594,58.1577981459126,0.0,-56.65779814591261,58.157798145912615,0.0]
        for i in range(len(rhs)):
            self.assertAlmostEqual(rhs[i],reference_res[i])

    def test_SDSimplestSurfaceLoadCondition3D4N(self):
        self._SimplestSurfaceLoadCondition3D4N("SmallDisplacement")

    def test_SimplestSurfaceLoadCondition3D4N(self):
        self._SimplestSurfaceLoadCondition3D4N()

    def test_SDSimpleOrientationSurfaceLoadCondition3D4N(self):
        self._SimpleOrientationSurfaceLoadCondition3D4N("SmallDisplacement")

    def test_SimpleOrientationSurfaceLoadCondition3D4N(self):
        self._SimpleOrientationSurfaceLoadCondition3D4N()

    def test_SDSimpleSurfaceLoadCondition3D4N(self):
        self._SimpleSurfaceLoadCondition3D4N("SmallDisplacement")

    def test_SimpleSurfaceLoadCondition3D4N(self):
        self._SimpleSurfaceLoadCondition3D4N()

    def test_SDSurfaceLoadCondition3D4N(self):
        self._SurfaceLoadCondition3D4N("SmallDisplacement")

    def test_SurfaceLoadCondition3D4N(self):
        self._SurfaceLoadCondition3D4N()

if __name__ == '__main__':
    KratosUnittest.main()
