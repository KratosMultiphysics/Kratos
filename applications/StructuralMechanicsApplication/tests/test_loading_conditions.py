from __future__ import print_function, absolute_import, division
import KratosMultiphysics 

import KratosMultiphysics.StructuralMechanicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest
import math


class TestLoadingConditions(KratosUnittest.TestCase):
    def setUp(self):
        pass
    
    

    def test_LineLoadCondition2D2N(self):
        dim = 2
        mp = KratosMultiphysics.ModelPart("solid_part")
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

        for node in mp.Nodes:
            node.AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X)
            node.AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y)
            node.AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z)
            
        cond = mp.CreateNewCondition("LineLoadCondition2D2N", 1, [1,2], mp.GetProperties()[1])
        
        lhs = KratosMultiphysics.Matrix(0,0)
        rhs = KratosMultiphysics.Vector(0)
        
        #first we apply a constant LINE_LOAD to theh condition 
        load_on_cond = KratosMultiphysics.Vector(3)
        load_on_cond[0] = 1.0
        load_on_cond[1] = 2.0
        load_on_cond[2] = 0.0 #note that this is a 2D condition
        cond.SetValue(KratosMultiphysics.StructuralMechanicsApplication.LINE_LOAD,load_on_cond)
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
 
if __name__ == '__main__':
    KratosUnittest.main()
