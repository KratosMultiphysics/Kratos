from __future__ import print_function, absolute_import, division

import KratosMultiphysics 
import KratosMultiphysics.SolidMechanicsApplication as SolidMechanicsApplication
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.ContactStructuralMechanicsApplication as ContactStructuralMechanicsApplication

import KratosMultiphysics.KratosUnittest as KratosUnittest

class TestModelPart(KratosUnittest.TestCase):
    
    #def __init__(self):
        ## TODO: Add if necessary
    
    ## Test exact integration in 2D
    # LINE
    def testLineExactIntegration(self):
        model_part = KratosMultiphysics.ModelPart("Main")
        model_part.SetBufferSize(3)
        model_part.AddProperties(KratosMultiphysics.Properties(1))
        
        normal = KratosMultiphysics.Vector(3)
        normal[0] = 0.0
        normal[1] = 1.0
        normal[2] = 0.0
        
        # Line 1
        model_part.CreateNewNode(1, 0.00,0.00,0.00)
        model_part.GetNode(1).SetValue(KratosMultiphysics.NORMAL, normal)
        model_part.CreateNewNode(2, 1.00,0.00,0.00)
        model_part.GetNode(2).SetValue(KratosMultiphysics.NORMAL, normal)
        
        cond1 = model_part.CreateNewCondition("Condition2D2N", 1, [1,2], model_part.GetProperties()[1])
        cond1.SetValue(KratosMultiphysics.NORMAL, normal)
        
        # Creating the utility:
        ExactIntegration = ContactStructuralMechanicsApplication.ExactMortarIntegrationUtility2D2N()

        # Line 2
        normal[1] = -1.0
        model_part.CreateNewNode(3, 0.50,0.001,0.00)
        model_part.GetNode(3).SetValue(KratosMultiphysics.NORMAL, normal)
        model_part.CreateNewNode(4, 1.50,0.001,0.00)
        model_part.GetNode(4).SetValue(KratosMultiphysics.NORMAL, normal)
        
        cond2 = model_part.CreateNewCondition("Condition2D2N", 2, [3,4], model_part.GetProperties()[1])
        cond2.SetValue(KratosMultiphysics.NORMAL, normal)
        
        MatrixSolution = KratosMultiphysics.Matrix()
        
        solution = ExactIntegration.TestGetExactIntegration(cond1, cond2, MatrixSolution)
        
        self.assertTrue(solution)
        self.assertEqual(MatrixSolution[0,0], 0.5)
        self.assertEqual(MatrixSolution[0,1], 1)
        
        ## Debug
        #if (solution == True):
            #print("Integration accomplished", MatrixSolution)
        
    ## Test exact integration in 3D
    # TRIANGLE
    def testTriangleExactIntegration(self):
        model_part = KratosMultiphysics.ModelPart("Main")
        model_part.SetBufferSize(3)
        model_part.AddProperties(KratosMultiphysics.Properties(1))
        
        normal = KratosMultiphysics.Vector(3)
        normal[0] = 0.0
        normal[1] = 0.0
        normal[2] = 1.0
        
        # Triangle 1
        model_part.CreateNewNode(1, 0.00,0.00,0.00)
        model_part.GetNode(1).SetValue(KratosMultiphysics.NORMAL, normal)
        model_part.CreateNewNode(2, 1.00,0.00,0.00)
        model_part.GetNode(2).SetValue(KratosMultiphysics.NORMAL, normal)
        model_part.CreateNewNode(3, 0.00,1.00,0.00)
        model_part.GetNode(3).SetValue(KratosMultiphysics.NORMAL, normal)
        
        cond1 = model_part.CreateNewCondition("SurfaceCondition3D3N", 1, [1,2,3], model_part.GetProperties()[1])
        cond1.SetValue(KratosMultiphysics.NORMAL, normal)
        
        # Creating the utility:
        ExactIntegration = ContactStructuralMechanicsApplication.ExactMortarIntegrationUtility3D3N()

        # Triangle 2
        normal[2] = -1.0
        model_part.CreateNewNode(4, 0.00,0.00,0.01)
        model_part.GetNode(4).SetValue(KratosMultiphysics.NORMAL, normal)
        model_part.CreateNewNode(5, 1.00,0.00,0.01)
        model_part.GetNode(5).SetValue(KratosMultiphysics.NORMAL, normal)
        model_part.CreateNewNode(6, 1.00,1.00,0.01)
        model_part.GetNode(6).SetValue(KratosMultiphysics.NORMAL, normal)
        
        cond2 = model_part.CreateNewCondition("SurfaceCondition3D3N", 2, [4,5,6], model_part.GetProperties()[1])
        cond2.SetValue(KratosMultiphysics.NORMAL, normal)
        
        MatrixSolution = KratosMultiphysics.Matrix()
        ExactIntegration.TestGetExactIntegration(cond1, cond2, MatrixSolution)
        
        ## Debug
        #if (solution == True):
            #print("Integration accomplished", MatrixSolution)

    # QUADRILATERAL
    def testQuadrilateralExactIntegration(self):
        model_part = KratosMultiphysics.ModelPart("Main")
        model_part.SetBufferSize(3)
        model_part.AddProperties(KratosMultiphysics.Properties(1))
        
        normal = KratosMultiphysics.Vector(3)
        normal[0] = 0.0
        normal[1] = 0.0
        normal[2] = 1.0
        
        # Quadrilateral 1
        model_part.CreateNewNode(1, 0.00,0.00,0.00)
        model_part.GetNode(1).SetValue(KratosMultiphysics.NORMAL, normal)
        model_part.CreateNewNode(2, 1.00,0.00,0.00)
        model_part.GetNode(2).SetValue(KratosMultiphysics.NORMAL, normal)
        model_part.CreateNewNode(3, 1.00,1.00,0.00)
        model_part.GetNode(3).SetValue(KratosMultiphysics.NORMAL, normal)
        model_part.CreateNewNode(4, 0.00,1.00,0.00)
        model_part.GetNode(4).SetValue(KratosMultiphysics.NORMAL, normal)
        
        cond1 = model_part.CreateNewCondition("Condition3D4N", 1, [1,2,3,4], model_part.GetProperties()[1])
        cond1.SetValue(KratosMultiphysics.NORMAL, normal)

        # Creating the utility:
        ExactIntegration = ContactStructuralMechanicsApplication.ExactMortarIntegrationUtility3D4N()

        # Quadrilateral 2
        normal[2] = -1.0
        model_part.CreateNewNode(5, 0.50,0.50,0.01)
        model_part.GetNode(5).SetValue(KratosMultiphysics.NORMAL, normal)
        model_part.CreateNewNode(6, 1.50,0.50,0.01)
        model_part.GetNode(6).SetValue(KratosMultiphysics.NORMAL, normal)
        model_part.CreateNewNode(7, 1.50,1.50,0.01)
        model_part.GetNode(7).SetValue(KratosMultiphysics.NORMAL, normal)
        model_part.CreateNewNode(8, 0.50,1.50,0.01)
        model_part.GetNode(8).SetValue(KratosMultiphysics.NORMAL, normal)
        
        cond2 = model_part.CreateNewCondition("Condition3D4N", 2, [5,6,7,8], model_part.GetProperties()[1])
        cond2.SetValue(KratosMultiphysics.NORMAL, normal)
        
        MatrixSolution = KratosMultiphysics.Matrix()
        #ExactIntegration.TestGetExactIntegration(cond1, cond2, MatrixSolution)
        
        ## Debug
        #if (solution == True):
            #print("Integration accomplished", MatrixSolution)
        
if __name__ == '__main__':
    KratosUnittest.main()
        
