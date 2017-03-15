from __future__ import print_function, absolute_import, division

import KratosMultiphysics 
import KratosMultiphysics.SolidMechanicsApplication as SolidMechanicsApplication
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.ContactStructuralMechanicsApplication as ContactStructuralMechanicsApplication

import KratosMultiphysics.KratosUnittest as KratosUnittest
    
## Test exact integration in 2D
# LINE
class TestLineExactIntegration1(KratosUnittest.TestCase):
    def test_execution(self):
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
        ExactIntegration = ContactStructuralMechanicsApplication.ExactMortarIntegrationUtility2D2N(1)

        # Line 2
        normal[1] = -1.0
        model_part.CreateNewNode(3, 0.00,0.001,0.00)
        model_part.GetNode(3).SetValue(KratosMultiphysics.NORMAL, normal)
        model_part.CreateNewNode(4, 1.00,0.001,0.00)
        model_part.GetNode(4).SetValue(KratosMultiphysics.NORMAL, normal)
        
        cond2 = model_part.CreateNewCondition("Condition2D2N", 2, [3,4], model_part.GetProperties()[1])
        cond2.SetValue(KratosMultiphysics.NORMAL, normal)
        
        MatrixSolution = KratosMultiphysics.Matrix()
        
        solution = ExactIntegration.TestGetExactIntegration(cond1, cond2, MatrixSolution)
        
        ## Debug
        #if (solution == True):
            #print("Integration accomplished", MatrixSolution)
        
        self.assertTrue(solution)
        self.assertEqual(MatrixSolution[0,0], 0)
        self.assertEqual(MatrixSolution[0,1], 2)
    
class TestLineExactIntegration2(KratosUnittest.TestCase):
    def test_execution(self):
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
        ExactIntegration = ContactStructuralMechanicsApplication.ExactMortarIntegrationUtility2D2N(1)

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
        
        ## Debug
        #if (solution == True):
            #print("Integration accomplished", MatrixSolution)
        
        self.assertTrue(solution)
        self.assertEqual(MatrixSolution[0,0], 0.5)
        self.assertEqual(MatrixSolution[0,1], 1)
    
## Test exact integration in 3D
# TRIANGLE
class TestTriangleExactIntegration1(KratosUnittest.TestCase):
    def test_execution(self):
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
        ExactIntegration = ContactStructuralMechanicsApplication.ExactMortarIntegrationUtility3D3N(2)

        # Triangle 2
        normal[2] = -1.0
        model_part.CreateNewNode(4, 0.00,0.00,0.01)
        model_part.GetNode(4).SetValue(KratosMultiphysics.NORMAL, normal)
        model_part.CreateNewNode(5, 1.00,0.00,0.01)
        model_part.GetNode(5).SetValue(KratosMultiphysics.NORMAL, normal)
        model_part.CreateNewNode(6, 0.00,1.00,0.01)
        model_part.GetNode(6).SetValue(KratosMultiphysics.NORMAL, normal)
        
        cond2 = model_part.CreateNewCondition("SurfaceCondition3D3N", 2, [4,5,6], model_part.GetProperties()[1])
        cond2.SetValue(KratosMultiphysics.NORMAL, normal)
        
        MatrixSolution = KratosMultiphysics.Matrix()
        solution = ExactIntegration.TestGetExactIntegration(cond1, cond2, MatrixSolution)
        
        ## Debug
        #if (solution == True):
            #print("Integration accomplished", MatrixSolution)
        
        self.assertTrue(solution)
        self.assertAlmostEqual(MatrixSolution[0,0], 1.0/6.0)
        self.assertAlmostEqual(MatrixSolution[0,1], 1.0/6.0)
        self.assertAlmostEqual(MatrixSolution[0,2], 1.0/6.0)
        self.assertAlmostEqual(MatrixSolution[1,0], 4.0/6.0)
        self.assertAlmostEqual(MatrixSolution[1,1], 1.0/6.0)
        self.assertAlmostEqual(MatrixSolution[1,2], 1.0/6.0)
        self.assertAlmostEqual(MatrixSolution[2,0], 1.0/6.0)
        self.assertAlmostEqual(MatrixSolution[2,1], 4.0/6.0)
        self.assertAlmostEqual(MatrixSolution[2,2], 1.0/6.0)
        
class TestTriangleExactIntegration2(KratosUnittest.TestCase):
    def test_execution(self):
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
        ExactIntegration = ContactStructuralMechanicsApplication.ExactMortarIntegrationUtility3D3N(2)

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
        solution = ExactIntegration.TestGetExactIntegration(cond1, cond2, MatrixSolution)
        
        ## Debug
        #if (solution == True):
            #print("Integration accomplished", MatrixSolution)
        
        self.assertTrue(solution)
        self.assertAlmostEqual(MatrixSolution[0,0], 0.25)
        self.assertAlmostEqual(MatrixSolution[0,1], 1.0/12.0)
        self.assertAlmostEqual(MatrixSolution[0,2], 1.0/12.0)
        self.assertAlmostEqual(MatrixSolution[1,0], 0.75)
        self.assertAlmostEqual(MatrixSolution[1,1], 1.0/12.0)
        self.assertAlmostEqual(MatrixSolution[1,2], 1.0/12.0)
        self.assertAlmostEqual(MatrixSolution[2,0], 0.5)
        self.assertAlmostEqual(MatrixSolution[2,1], 1.0/3.0)
        self.assertAlmostEqual(MatrixSolution[2,2], 1.0/12.0)
        
class TestTriangleExactIntegration3(KratosUnittest.TestCase):
    def test_execution(self):
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
        model_part.CreateNewNode(4, 1.00,1.00,0.00)
        model_part.GetNode(3).SetValue(KratosMultiphysics.NORMAL, normal)
        
        cond1 = model_part.CreateNewCondition("SurfaceCondition3D3N", 1, [1,2,3], model_part.GetProperties()[1])
        cond2 = model_part.CreateNewCondition("SurfaceCondition3D3N", 2, [2,3,4], model_part.GetProperties()[1])
        cond1.SetValue(KratosMultiphysics.NORMAL, normal)
        cond2.SetValue(KratosMultiphysics.NORMAL, normal)
        
        # Creating the utility:
        ExactIntegration = ContactStructuralMechanicsApplication.ExactMortarIntegrationUtility3D3N(2)

        # Triangle 2
        normal[2] = -1.0
        model_part.CreateNewNode(5, 0.00,0.00,0.01)
        model_part.GetNode(5).SetValue(KratosMultiphysics.NORMAL, normal)
        model_part.CreateNewNode(6, 1.00,0.00,0.01)
        model_part.GetNode(6).SetValue(KratosMultiphysics.NORMAL, normal)
        model_part.CreateNewNode(7, 0.00,1.00,0.01)
        model_part.GetNode(7).SetValue(KratosMultiphysics.NORMAL, normal)
        model_part.CreateNewNode(8, 1.00,1.00,0.01)
        model_part.GetNode(7).SetValue(KratosMultiphysics.NORMAL, normal)
        
        cond3 = model_part.CreateNewCondition("SurfaceCondition3D3N", 3, [5,6,7], model_part.GetProperties()[1])
        cond4 = model_part.CreateNewCondition("SurfaceCondition3D3N", 4, [6,7,8], model_part.GetProperties()[1])
        cond3.SetValue(KratosMultiphysics.NORMAL, normal)
        cond4.SetValue(KratosMultiphysics.NORMAL, normal)
        
        MatrixSolution = KratosMultiphysics.Matrix()
        
        
        solution = ExactIntegration.TestGetExactIntegration(cond1, cond3, MatrixSolution)
        
        # Debug
        if (solution == True):
            print("\n\nFirst Integration accomplished", MatrixSolution)
        
        solution = ExactIntegration.TestGetExactIntegration(cond1, cond4, MatrixSolution)
        
        # Debug
        if (solution == True):
            print("\n\nSecond Integration accomplished", MatrixSolution)
            
        solution = ExactIntegration.TestGetExactIntegration(cond2, cond3, MatrixSolution)
        
        # Debug
        if (solution == True):
            print("\n\nThird Integration accomplished", MatrixSolution)
        
        solution = ExactIntegration.TestGetExactIntegration(cond2, cond4, MatrixSolution)
        
        # Debug
        if (solution == True):
            print("\n\nFourth Integration accomplished", MatrixSolution)
        
        #self.assertTrue(solution)
        #self.assertAlmostEqual(MatrixSolution[0,0], 0.25)
        #self.assertAlmostEqual(MatrixSolution[0,1], 1.0/12.0)
        #self.assertAlmostEqual(MatrixSolution[0,2], 1.0/12.0)
        #self.assertAlmostEqual(MatrixSolution[1,0], 0.75)
        #self.assertAlmostEqual(MatrixSolution[1,1], 1.0/12.0)
        #self.assertAlmostEqual(MatrixSolution[1,2], 1.0/12.0)
        #self.assertAlmostEqual(MatrixSolution[2,0], 0.5)
        #self.assertAlmostEqual(MatrixSolution[2,1], 1.0/3.0)
        #self.assertAlmostEqual(MatrixSolution[2,2], 1.0/12.0)

# QUADRILATERAL
class TestQuadrilateralExactIntegration1(KratosUnittest.TestCase):
    def test_execution(self):
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
        ExactIntegration = ContactStructuralMechanicsApplication.ExactMortarIntegrationUtility3D4N(2)

        # Quadrilateral 2
        normal[2] = -1.0
        model_part.CreateNewNode(5, 0.00,0.00,0.01)
        model_part.GetNode(5).SetValue(KratosMultiphysics.NORMAL, normal)
        model_part.CreateNewNode(6, 1.00,0.00,0.01)
        model_part.GetNode(6).SetValue(KratosMultiphysics.NORMAL, normal)
        model_part.CreateNewNode(7, 1.00,1.00,0.01)
        model_part.GetNode(7).SetValue(KratosMultiphysics.NORMAL, normal)
        model_part.CreateNewNode(8, 0.00,1.00,0.01)
        model_part.GetNode(8).SetValue(KratosMultiphysics.NORMAL, normal)
        
        cond2 = model_part.CreateNewCondition("Condition3D4N", 2, [5,6,7,8], model_part.GetProperties()[1])
        cond2.SetValue(KratosMultiphysics.NORMAL, normal)
        
        MatrixSolution = KratosMultiphysics.Matrix()
        solution = ExactIntegration.TestGetExactIntegration(cond1, cond2, MatrixSolution)
        
        ## Debug
        #if (solution == True):
            #print("Integration accomplished", MatrixSolution)
        
        self.assertTrue(solution)
        self.assertAlmostEqual(MatrixSolution[0,0], -1.0/3.0)
        self.assertAlmostEqual(MatrixSolution[0,1], -2.0/3.0)
        self.assertAlmostEqual(MatrixSolution[0,2],  2.0/3.0)
        self.assertAlmostEqual(MatrixSolution[1,0],  2.0/3.0)
        self.assertAlmostEqual(MatrixSolution[1,1], -2.0/3.0)
        self.assertAlmostEqual(MatrixSolution[1,2],  2.0/3.0)
        self.assertAlmostEqual(MatrixSolution[2,0],  2.0/3.0)
        self.assertAlmostEqual(MatrixSolution[2,1],  1.0/3.0)
        self.assertAlmostEqual(MatrixSolution[2,2],  2.0/3.0)
        self.assertAlmostEqual(MatrixSolution[3,0], -2.0/3.0)
        self.assertAlmostEqual(MatrixSolution[3,1], -1.0/3.0)
        self.assertAlmostEqual(MatrixSolution[3,2],  2.0/3.0)
        self.assertAlmostEqual(MatrixSolution[4,0],  1.0/3.0)
        self.assertAlmostEqual(MatrixSolution[4,1],  2.0/3.0)
        self.assertAlmostEqual(MatrixSolution[4,2],  2.0/3.0)
        self.assertAlmostEqual(MatrixSolution[5,0], -2.0/3.0)
        self.assertAlmostEqual(MatrixSolution[5,1],  2.0/3.0)
        self.assertAlmostEqual(MatrixSolution[5,2],  2.0/3.0)
        
class TestQuadrilateralExactIntegration2(KratosUnittest.TestCase):
    def test_execution(self):
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
        solution = ExactIntegration.TestGetExactIntegration(cond1, cond2, MatrixSolution)
        
        ## Debug
        #if (solution == True):
            #print("Integration accomplished", MatrixSolution)
        
        self.assertTrue(solution)
        self.assertAlmostEqual(MatrixSolution[0,0], 2.0/6.0)
        self.assertAlmostEqual(MatrixSolution[0,1], 1.0/6.0)
        self.assertAlmostEqual(MatrixSolution[0,2], 1.0/6.0)
        self.assertAlmostEqual(MatrixSolution[1,0], 5.0/6.0)
        self.assertAlmostEqual(MatrixSolution[1,1], 1.0/6.0)
        self.assertAlmostEqual(MatrixSolution[1,2], 1.0/6.0)
        self.assertAlmostEqual(MatrixSolution[2,0], 5.0/6.0)
        self.assertAlmostEqual(MatrixSolution[2,1], 4.0/6.0)
        self.assertAlmostEqual(MatrixSolution[2,2], 1.0/6.0)
        self.assertAlmostEqual(MatrixSolution[3,0], 1.0/6.0)
        self.assertAlmostEqual(MatrixSolution[3,1], 2.0/6.0)
        self.assertAlmostEqual(MatrixSolution[3,2], 1.0/6.0)
        self.assertAlmostEqual(MatrixSolution[4,0], 4.0/6.0)
        self.assertAlmostEqual(MatrixSolution[4,1], 5.0/6.0)
        self.assertAlmostEqual(MatrixSolution[4,2], 1.0/6.0)
        self.assertAlmostEqual(MatrixSolution[5,0], 1.0/6.0)
        self.assertAlmostEqual(MatrixSolution[5,1], 5.0/6.0)
        self.assertAlmostEqual(MatrixSolution[5,2], 1.0/6.0)
            
if __name__ == '__main__':
   TestLineExactIntegration1().test_execution()
   TestLineExactIntegration2().test_execution()
   TestTriangleExactIntegration1().test_execution()
   TestTriangleExactIntegration2().test_execution()
   TestTriangleExactIntegration3().test_execution()
   TestQuadrilateralExactIntegration1().test_execution()
   TestQuadrilateralExactIntegration2().test_execution()
        
