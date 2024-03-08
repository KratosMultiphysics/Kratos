import KratosMultiphysics

import KratosMultiphysics.MPMApplication as KratosMPM
import KratosMultiphysics.KratosUnittest as KratosUnittest
from math import sqrt


class TestStaticLoadingConditionsSurface(KratosUnittest.TestCase):

    def test_MPMGridSurfaceLoadCondition3D4N(self):
        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("solid_part")
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.POSITIVE_FACE_PRESSURE)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.NEGATIVE_FACE_PRESSURE)

        # Create nodes
        mp.CreateNewNode(1,0.0,0.0,0.0)
        mp.CreateNewNode(2,0.0,0.0,1.0)
        mp.CreateNewNode(3,sqrt(2),sqrt(2.0),0.0)
        mp.CreateNewNode(4,sqrt(2),sqrt(2.0),1.0)
        length = 2.0
        width = 1.0

        # Ensure that the property 1 is created
        mp.GetProperties()[1]

        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z,mp)


        cond = mp.CreateNewCondition("MPMGridSurfaceLoadCondition3D4N", 1, [1,2,4,3], mp.GetProperties()[1])

        lhs = KratosMultiphysics.Matrix(0,0)
        rhs = KratosMultiphysics.Vector(0)

        # First we apply a constant SURFACE_LOAD to then condition
        load_on_cond = KratosMultiphysics.Vector(3)
        load_on_cond[0] = 1.0
        load_on_cond[1] = 2.0
        load_on_cond[2] = 0.0 # Note that this is a 2D condition
        cond.SetValue(KratosMPM.SURFACE_LOAD,load_on_cond)
        cond.CalculateLocalSystem(lhs,rhs,mp.ProcessInfo)
        self.assertAlmostEqual(rhs[0],0.25*length); self.assertAlmostEqual(rhs[1],0.5*1.0*length); self.assertAlmostEqual(rhs[2],0.0)
        self.assertAlmostEqual(rhs[3],0.25*length); self.assertAlmostEqual(rhs[4],0.5*1.0*length); self.assertAlmostEqual(rhs[5],0.0)
        self.assertAlmostEqual(rhs[6],0.25*length); self.assertAlmostEqual(rhs[7],0.5*1.0*length); self.assertAlmostEqual(rhs[8],0.0)
        self.assertAlmostEqual(rhs[9],0.25*length); self.assertAlmostEqual(rhs[10],0.5*1.0*length); self.assertAlmostEqual(rhs[11],0.0)

        cond.SetValue(KratosMultiphysics.POSITIVE_FACE_PRESSURE,10.0)
        cond.SetValue(KratosMultiphysics.NEGATIVE_FACE_PRESSURE,5.0)
        cond.CalculateLocalSystem(lhs,rhs,mp.ProcessInfo)
        nodal_fx = (1.0 + 5.0*sqrt(2.0)/2.0) * length*width / 4.0
        nodal_fy = (2.0 - 5.0*sqrt(2.0)/2.0) * length*width / 4.0
        for i in range(4):
            base = i*3
            self.assertAlmostEqual(rhs[base+0],nodal_fx)
            self.assertAlmostEqual(rhs[base+1],nodal_fy)
            self.assertAlmostEqual(rhs[base+2],0.0)

        ## Finally we apply TO THE NODES, a linearly varying POSITIVE_FACE_PRESSURE ranging from -100.0 to -200.0
        mp.Nodes[1].SetSolutionStepValue(KratosMultiphysics.POSITIVE_FACE_PRESSURE,0,-100.0)
        mp.Nodes[2].SetSolutionStepValue(KratosMultiphysics.POSITIVE_FACE_PRESSURE,0,-100.0)
        mp.Nodes[3].SetSolutionStepValue(KratosMultiphysics.POSITIVE_FACE_PRESSURE,0,-200.0)
        mp.Nodes[4].SetSolutionStepValue(KratosMultiphysics.POSITIVE_FACE_PRESSURE,0,-200.0)
        cond.CalculateLocalSystem(lhs,rhs,mp.ProcessInfo)
        reference_res = [-44.872685126136794,46.372685126136815,0.0,-44.87268512613681,46.3726851261368,0.0,-56.657798145912594,58.1577981459126,0.0,-56.65779814591261,58.157798145912615,0.0]
        for i in range(len(rhs)):
            self.assertAlmostEqual(rhs[i],reference_res[i])

if __name__ == '__main__':
    KratosUnittest.main()
