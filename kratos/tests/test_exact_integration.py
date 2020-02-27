from __future__ import print_function, absolute_import, division

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest


class TestExactIntegration(KratosUnittest.TestCase):

    def setUp(self):
        pass
    

    # Test exact integration in 2D
    # LINE
    def test_line_exact_integration_1(self):
        current_model = KratosMultiphysics.Model()

        model_part = current_model.CreateModelPart("Main")
        model_part.SetBufferSize(3)
        model_part.AddProperties(KratosMultiphysics.Properties(1))
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)

        normal = KratosMultiphysics.Vector(3)

        # Line 1
        model_part.CreateNewNode(1, 0.00, 0.00, 0.00)
        model_part.CreateNewNode(2, 1.00, 0.00, 0.00)

        cond1 = model_part.CreateNewCondition("LineCondition2D2N", 1, [1, 2], model_part.GetProperties()[1])
        normal = cond1.GetGeometry().UnitNormal()
        cond1.SetValue(KratosMultiphysics.NORMAL, normal)
        model_part.GetNode(1).SetSolutionStepValue(KratosMultiphysics.NORMAL, normal)
        model_part.GetNode(2).SetSolutionStepValue(KratosMultiphysics.NORMAL, normal)

        # Creating the utility:
        exact_integration = KratosMultiphysics.ExactMortarIntegrationUtility2D2N(2)

        # Line 2
        model_part.CreateNewNode(3, 0.00, 0.001, 0.00)
        model_part.CreateNewNode(4, 1.00, 0.001, 0.00)

        cond2 = model_part.CreateNewCondition("LineCondition2D2N", 2, [3, 4], model_part.GetProperties()[1])
        normal = cond2.GetGeometry().UnitNormal()
        cond2.SetValue(KratosMultiphysics.NORMAL, normal)
        model_part.GetNode(3).SetSolutionStepValue(KratosMultiphysics.NORMAL, normal)
        model_part.GetNode(4).SetSolutionStepValue(KratosMultiphysics.NORMAL, normal)

        matrix_solution = KratosMultiphysics.Matrix()

        solution = exact_integration.TestGetExactIntegration(cond1, cond2, matrix_solution)

        # Debug
        #if (solution == True):
            #print("Integration accomplished", matrix_solution)

        self.assertTrue(solution)
        self.assertAlmostEqual(matrix_solution[0, 0], -0.57735026918963)
        self.assertAlmostEqual(matrix_solution[0, 1], 1.0)
        self.assertAlmostEqual(matrix_solution[1, 0],  0.57735026918963)
        self.assertAlmostEqual(matrix_solution[1, 1], 1.0)

    def test_line_exact_integration_2(self):
        current_model = KratosMultiphysics.Model()

        model_part = current_model.CreateModelPart("Main")
        model_part.SetBufferSize(3)
        model_part.AddProperties(KratosMultiphysics.Properties(1))
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)

        normal = KratosMultiphysics.Vector(3)

        # Line 1
        model_part.CreateNewNode(1, 0.00, 0.00, 0.00)
        model_part.CreateNewNode(2, 1.00, 0.00, 0.00)

        cond1 = model_part.CreateNewCondition("LineCondition2D2N", 1, [1, 2], model_part.GetProperties()[1])
        normal = cond1.GetGeometry().UnitNormal()
        cond1.SetValue(KratosMultiphysics.NORMAL, normal)
        model_part.GetNode(1).SetSolutionStepValue(KratosMultiphysics.NORMAL, normal)
        model_part.GetNode(2).SetSolutionStepValue(KratosMultiphysics.NORMAL, normal)

        # Creating the utility:
        exact_integration = KratosMultiphysics.ExactMortarIntegrationUtility2D2N(2)

        # Line 2
        model_part.CreateNewNode(3, 0.50, 0.001, 0.00)
        model_part.CreateNewNode(4, 1.50, 0.001, 0.00)

        cond2 = model_part.CreateNewCondition("LineCondition2D2N", 2, [3, 4], model_part.GetProperties()[1])
        normal = cond2.GetGeometry().UnitNormal()
        cond2.SetValue(KratosMultiphysics.NORMAL, normal)
        model_part.GetNode(3).SetSolutionStepValue(KratosMultiphysics.NORMAL, normal)
        model_part.GetNode(4).SetSolutionStepValue(KratosMultiphysics.NORMAL, normal)

        matrix_solution = KratosMultiphysics.Matrix()

        solution = exact_integration.TestGetExactIntegration(cond1, cond2, matrix_solution)

        # Debug
        #if (solution == True):
            #print("Integration accomplished", matrix_solution)

        self.assertTrue(solution)
        self.assertAlmostEqual(matrix_solution[0, 0], 0.21132486540517492)
        self.assertAlmostEqual(matrix_solution[0, 1], 0.5)
        self.assertAlmostEqual(matrix_solution[1, 0], 0.7886751345947951)
        self.assertAlmostEqual(matrix_solution[1, 1], 0.5)

    def test_line_exact_integration_3(self):
        current_model = KratosMultiphysics.Model()

        model_part = current_model.CreateModelPart("Main")
        model_part.SetBufferSize(3)
        model_part.AddProperties(KratosMultiphysics.Properties(1))
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)

        normal = KratosMultiphysics.Vector(3)

        # Line 1
        model_part.CreateNewNode(1, 0.00, -0.5, 0.00)
        model_part.CreateNewNode(2, 1.00,  0.5, 0.00)

        cond1 = model_part.CreateNewCondition("LineCondition2D2N", 1, [1, 2], model_part.GetProperties()[1])

        normal = cond1.GetGeometry().UnitNormal()
        cond1.SetValue(KratosMultiphysics.NORMAL, normal)
        model_part.GetNode(1).SetSolutionStepValue(KratosMultiphysics.NORMAL, normal)
        model_part.GetNode(2).SetSolutionStepValue(KratosMultiphysics.NORMAL, normal)

        # Creating the utility:
        exact_integration = KratosMultiphysics.ExactMortarIntegrationUtility2D2N(2)

        # Line 2
        model_part.CreateNewNode(3, 0.0, 0.5, 0.00)
        model_part.CreateNewNode(4, 1.0, 0.5, 0.00)

        cond2 = model_part.CreateNewCondition("LineCondition2D2N", 2, [3, 4], model_part.GetProperties()[1])

        normal = cond2.GetGeometry().UnitNormal()
        cond2.SetValue(KratosMultiphysics.NORMAL, normal)
        model_part.GetNode(3).SetSolutionStepValue(KratosMultiphysics.NORMAL, normal)
        model_part.GetNode(4).SetSolutionStepValue(KratosMultiphysics.NORMAL, normal)

        matrix_solution = KratosMultiphysics.Matrix()

        solution = exact_integration.TestGetExactIntegration(cond1, cond2, matrix_solution)

        # Debug
        #if (solution == True):
            #print("Integration accomplished", matrix_solution)

        self.assertTrue(solution)
        self.assertAlmostEqual(matrix_solution[0, 0], 0.21132486540517847)
        self.assertAlmostEqual(matrix_solution[1, 0], 0.7886751345948)
        self.assertAlmostEqual(matrix_solution[0, 1], 1.0 / (2.0**0.5))
        self.assertAlmostEqual(matrix_solution[1, 1], 1.0 / (2.0**0.5))

    # Test exact integration in 3D
    # TRIANGLE
    def test_triangle_exact_integration_1(self):
        current_model = KratosMultiphysics.Model()

        model_part = current_model.CreateModelPart("Main")
        model_part.SetBufferSize(3)
        model_part.AddProperties(KratosMultiphysics.Properties(1))
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)

        normal = KratosMultiphysics.Vector(3)

        # Triangle 1
        model_part.CreateNewNode(1, 0.00, 0.00, 0.00)
        model_part.CreateNewNode(2, 1.00, 0.00, 0.00)
        model_part.CreateNewNode(3, 0.00, 1.00, 0.00)

        cond1 = model_part.CreateNewCondition("SurfaceCondition3D3N", 1, [1, 2, 3], model_part.GetProperties()[1])
        normal = cond1.GetGeometry().UnitNormal()
        cond1.SetValue(KratosMultiphysics.NORMAL, normal)
        model_part.GetNode(1).SetSolutionStepValue(KratosMultiphysics.NORMAL, normal)
        model_part.GetNode(2).SetSolutionStepValue(KratosMultiphysics.NORMAL, normal)
        model_part.GetNode(3).SetSolutionStepValue(KratosMultiphysics.NORMAL, normal)

        # Creating the utility:
        exact_integration = KratosMultiphysics.ExactMortarIntegrationUtility3D3N(2)

        # Triangle 2
        normal[2] = -1.0
        model_part.CreateNewNode(4, 0.00, 0.00, 0.01)
        model_part.CreateNewNode(5, 1.00, 0.00, 0.01)
        model_part.CreateNewNode(6, 0.00, 1.00, 0.01)

        cond2 = model_part.CreateNewCondition("SurfaceCondition3D3N", 2, [4, 5, 6], model_part.GetProperties()[1])
        normal = cond2.GetGeometry().UnitNormal()
        cond2.SetValue(KratosMultiphysics.NORMAL, normal)
        model_part.GetNode(4).SetSolutionStepValue(KratosMultiphysics.NORMAL, normal)
        model_part.GetNode(5).SetSolutionStepValue(KratosMultiphysics.NORMAL, normal)
        model_part.GetNode(6).SetSolutionStepValue(KratosMultiphysics.NORMAL, normal)

        matrix_solution = KratosMultiphysics.Matrix()
        solution = exact_integration.TestGetExactIntegration(cond1, cond2, matrix_solution)

        # Debug
        #if (solution == True):
            #print("Integration accomplished", matrix_solution)

        self.assertTrue(solution)
        self.assertAlmostEqual(matrix_solution[0, 0], 1.0 / 6.0)
        self.assertAlmostEqual(matrix_solution[0, 1], 1.0 / 6.0)
        self.assertAlmostEqual(matrix_solution[0, 2], 1.0 / 6.0)
        self.assertAlmostEqual(matrix_solution[1, 0], 4.0 / 6.0)
        self.assertAlmostEqual(matrix_solution[1, 1], 1.0 / 6.0)
        self.assertAlmostEqual(matrix_solution[1, 2], 1.0 / 6.0)
        self.assertAlmostEqual(matrix_solution[2, 0], 1.0 / 6.0)
        self.assertAlmostEqual(matrix_solution[2, 1], 4.0 / 6.0)
        self.assertAlmostEqual(matrix_solution[2, 2], 1.0 / 6.0)

    def test_triangle_exact_integration_2(self):
        current_model = KratosMultiphysics.Model()

        model_part = current_model.CreateModelPart("Main")
        model_part.SetBufferSize(3)
        model_part.AddProperties(KratosMultiphysics.Properties(1))
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)

        normal = KratosMultiphysics.Vector(3)

        # Triangle 1
        model_part.CreateNewNode(1, 0.00, 0.00, 0.00)
        model_part.CreateNewNode(2, 1.00, 0.00, 0.00)
        model_part.CreateNewNode(3, 0.00, 1.00, 0.00)

        cond1 = model_part.CreateNewCondition("SurfaceCondition3D3N", 1, [1, 2, 3], model_part.GetProperties()[1])
        normal = cond1.GetGeometry().UnitNormal()
        cond1.SetValue(KratosMultiphysics.NORMAL, normal)
        model_part.GetNode(1).SetSolutionStepValue(KratosMultiphysics.NORMAL, normal)
        model_part.GetNode(2).SetSolutionStepValue(KratosMultiphysics.NORMAL, normal)
        model_part.GetNode(3).SetSolutionStepValue(KratosMultiphysics.NORMAL, normal)

        # Creating the utility:
        exact_integration = KratosMultiphysics.ExactMortarIntegrationUtility3D3N(2)

        # Triangle 2
        model_part.CreateNewNode(4, 0.00, 0.00, 0.01)
        model_part.CreateNewNode(5, 1.00, 0.00, 0.01)
        model_part.CreateNewNode(6, 1.00, 1.00, 0.01)

        cond2 = model_part.CreateNewCondition("SurfaceCondition3D3N", 2, [4, 5, 6], model_part.GetProperties()[1])
        normal = cond2.GetGeometry().UnitNormal()
        cond2.SetValue(KratosMultiphysics.NORMAL, normal)
        model_part.GetNode(4).SetSolutionStepValue(KratosMultiphysics.NORMAL, normal)
        model_part.GetNode(5).SetSolutionStepValue(KratosMultiphysics.NORMAL, normal)
        model_part.GetNode(6).SetSolutionStepValue(KratosMultiphysics.NORMAL, normal)

        matrix_solution = KratosMultiphysics.Matrix()
        solution = exact_integration.TestGetExactIntegration(cond1, cond2, matrix_solution)

        # Debug
        #if (solution == True):
            #print("Integration accomplished", matrix_solution)

        self.assertTrue(solution)
        self.assertAlmostEqual(matrix_solution[0, 0], 0.25)
        self.assertAlmostEqual(matrix_solution[0, 1], 1.0 / 12.0)
        self.assertAlmostEqual(matrix_solution[0, 2], 1.0 / 12.0)
        self.assertAlmostEqual(matrix_solution[1, 0], 0.75)
        self.assertAlmostEqual(matrix_solution[1, 1], 1.0 / 12.0)
        self.assertAlmostEqual(matrix_solution[1, 2], 1.0 / 12.0)
        self.assertAlmostEqual(matrix_solution[2, 0], 0.5)
        self.assertAlmostEqual(matrix_solution[2, 1], 1.0 / 3.0)
        self.assertAlmostEqual(matrix_solution[2, 2], 1.0 / 12.0)

    def test_triangle_exact_integration_3(self):
        current_model = KratosMultiphysics.Model()

        model_part = current_model.CreateModelPart("Main")
        model_part.SetBufferSize(3)
        model_part.AddProperties(KratosMultiphysics.Properties(1))
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)

        normal = KratosMultiphysics.Vector(3)

        # Triangle 1 and 2
        model_part.CreateNewNode(1, 0.00, 0.00, 0.00)
        model_part.CreateNewNode(2, 1.00, 0.00, 0.00)
        model_part.CreateNewNode(3, 0.00, 1.00, 0.00)
        model_part.CreateNewNode(4, 1.00, 1.00, 0.00)

        cond1 = model_part.CreateNewCondition("SurfaceCondition3D3N", 1, [1, 2, 3], model_part.GetProperties()[1])
        cond2 = model_part.CreateNewCondition("SurfaceCondition3D3N", 2, [2, 4, 3], model_part.GetProperties()[1])
        normal = cond1.GetGeometry().UnitNormal()
        cond1.SetValue(KratosMultiphysics.NORMAL, normal)
        cond2.SetValue(KratosMultiphysics.NORMAL, normal)
        model_part.GetNode(1).SetSolutionStepValue(KratosMultiphysics.NORMAL, normal)
        model_part.GetNode(2).SetSolutionStepValue(KratosMultiphysics.NORMAL, normal)
        model_part.GetNode(3).SetSolutionStepValue(KratosMultiphysics.NORMAL, normal)
        model_part.GetNode(4).SetSolutionStepValue(KratosMultiphysics.NORMAL, normal)

        # Creating the utility:
        exact_integration = KratosMultiphysics.ExactMortarIntegrationUtility3D3N(2)

        # Triangle 3 and 4
        model_part.CreateNewNode(5, 0.00, 0.00, 0.01)
        model_part.CreateNewNode(6, 1.00, 0.00, 0.01)
        model_part.CreateNewNode(7, 0.00, 1.00, 0.01)
        model_part.CreateNewNode(8, 1.00, 1.00, 0.01)

        cond3 = model_part.CreateNewCondition("SurfaceCondition3D3N", 3, [5, 6, 8], model_part.GetProperties()[1])
        cond4 = model_part.CreateNewCondition("SurfaceCondition3D3N", 4, [5, 8, 7], model_part.GetProperties()[1])
        normal = cond3.GetGeometry().UnitNormal()
        cond3.SetValue(KratosMultiphysics.NORMAL, normal)
        cond4.SetValue(KratosMultiphysics.NORMAL, normal)
        model_part.GetNode(5).SetSolutionStepValue(KratosMultiphysics.NORMAL, normal)
        model_part.GetNode(6).SetSolutionStepValue(KratosMultiphysics.NORMAL, normal)
        model_part.GetNode(7).SetSolutionStepValue(KratosMultiphysics.NORMAL, normal)
        model_part.GetNode(8).SetSolutionStepValue(KratosMultiphysics.NORMAL, normal)

        matrix_solution = KratosMultiphysics.Matrix()

        solution = exact_integration.TestGetExactIntegration(cond1, cond3, matrix_solution)

        # Debug
        #if (solution == True):
            #print("First Integration accomplished", matrix_solution)

        self.assertTrue(solution)
        self.assertAlmostEqual(matrix_solution[0, 0], 0.25)
        self.assertAlmostEqual(matrix_solution[0, 1], 1.0 / 12.0)
        self.assertAlmostEqual(matrix_solution[0, 2], 1.0 / 12.0)
        self.assertAlmostEqual(matrix_solution[1, 0], 0.75)
        self.assertAlmostEqual(matrix_solution[1, 1], 1.0 / 12.0)
        self.assertAlmostEqual(matrix_solution[1, 2], 1.0 / 12.0)
        self.assertAlmostEqual(matrix_solution[2, 0], 0.5)
        self.assertAlmostEqual(matrix_solution[2, 1], 1.0 / 3.0)
        self.assertAlmostEqual(matrix_solution[2, 2], 1.0 / 12.0)

        solution = exact_integration.TestGetExactIntegration(cond1, cond4, matrix_solution)

        # Debug
        #if (solution == True):
            #print("Second Integration accomplished", matrix_solution)

        self.assertTrue(solution)
        self.assertAlmostEqual(matrix_solution[0, 0], 1.0 / 12.0)
        self.assertAlmostEqual(matrix_solution[0, 1], 0.25)
        self.assertAlmostEqual(matrix_solution[0, 2], 1.0 / 12.0)
        self.assertAlmostEqual(matrix_solution[1, 0], 1.0 / 3.0)
        self.assertAlmostEqual(matrix_solution[1, 1], 0.5)
        self.assertAlmostEqual(matrix_solution[1, 2], 1.0 / 12.0)
        self.assertAlmostEqual(matrix_solution[2, 0], 1.0 / 12.0)
        self.assertAlmostEqual(matrix_solution[2, 1], 0.75)
        self.assertAlmostEqual(matrix_solution[2, 2], 1.0 / 12.0)

        solution = exact_integration.TestGetExactIntegration(cond2, cond3, matrix_solution)

        # Debug
        #if (solution == True):
            #print("Third Integration accomplished", matrix_solution)

        self.assertTrue(solution)
        self.assertAlmostEqual(matrix_solution[0, 0],  1.0 / 6.0)
        self.assertAlmostEqual(matrix_solution[0, 1], 1.0 / 12.0)
        self.assertAlmostEqual(matrix_solution[0, 2],  1.0 / 12.0)
        self.assertAlmostEqual(matrix_solution[1, 0], 2.0 / 3.0)
        self.assertAlmostEqual(matrix_solution[1, 1], 1.0 / 12.0)
        self.assertAlmostEqual(matrix_solution[1, 2], 1.0 / 12.0)
        self.assertAlmostEqual(matrix_solution[2, 0], 1.0 / 6.0)
        self.assertAlmostEqual(matrix_solution[2, 1], 1.0 / 3.0)
        self.assertAlmostEqual(matrix_solution[2, 2], 1.0 / 12.0)

        solution = exact_integration.TestGetExactIntegration(cond2, cond4, matrix_solution)

        # Debug
        #if (solution == True):
            #print("Fourth Integration accomplished", matrix_solution)

        self.assertTrue(solution)
        self.assertAlmostEqual(matrix_solution[0, 0],  4.0 / 6.0)
        self.assertAlmostEqual(matrix_solution[0, 1],  0.25)
        self.assertAlmostEqual(matrix_solution[0, 2],  1.0 / 12.0)
        self.assertAlmostEqual(matrix_solution[1, 0],  1.0 / 6.0)
        self.assertAlmostEqual(matrix_solution[1, 1],  0.75)
        self.assertAlmostEqual(matrix_solution[1, 2],  1.0 / 12.0)
        self.assertAlmostEqual(matrix_solution[2, 0],  1.0 / 6.0)
        self.assertAlmostEqual(matrix_solution[2, 1],  1.0 / 2.0)
        self.assertAlmostEqual(matrix_solution[2, 2],  1.0 / 12.0)

    # QUADRILATERAL
    def test_quadrilateral_exact_integration_1(self):
        current_model = KratosMultiphysics.Model()

        model_part = current_model.CreateModelPart("Main")
        model_part.SetBufferSize(3)
        model_part.AddProperties(KratosMultiphysics.Properties(1))
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)

        normal = KratosMultiphysics.Vector(3)

        # Quadrilateral 1
        model_part.CreateNewNode(1, 0.00, 0.00, 0.00)
        model_part.CreateNewNode(2, 1.00, 0.00, 0.00)
        model_part.CreateNewNode(3, 1.00, 1.00, 0.00)
        model_part.CreateNewNode(4, 0.00, 1.00, 0.00)

        cond1 = model_part.CreateNewCondition("SurfaceCondition3D4N", 1, [1, 2, 3, 4], model_part.GetProperties()[1])
        normal = cond1.GetGeometry().UnitNormal()
        cond1.SetValue(KratosMultiphysics.NORMAL, normal)
        model_part.GetNode(1).SetSolutionStepValue(KratosMultiphysics.NORMAL, normal)
        model_part.GetNode(2).SetSolutionStepValue(KratosMultiphysics.NORMAL, normal)
        model_part.GetNode(3).SetSolutionStepValue(KratosMultiphysics.NORMAL, normal)
        model_part.GetNode(4).SetSolutionStepValue(KratosMultiphysics.NORMAL, normal)

        # Creating the utility:
        exact_integration = KratosMultiphysics.ExactMortarIntegrationUtility3D4N(2)

        # Quadrilateral 2
        model_part.CreateNewNode(5, 0.00, 0.00, 0.01)
        model_part.CreateNewNode(6, 1.00, 0.00, 0.01)
        model_part.CreateNewNode(7, 1.00, 1.00, 0.01)
        model_part.CreateNewNode(8, 0.00, 1.00, 0.01)

        cond2 = model_part.CreateNewCondition("SurfaceCondition3D4N", 2, [5, 6, 7, 8], model_part.GetProperties()[1])
        normal = cond2.GetGeometry().UnitNormal()
        cond2.SetValue(KratosMultiphysics.NORMAL, normal)
        model_part.GetNode(5).SetSolutionStepValue(KratosMultiphysics.NORMAL, normal)
        model_part.GetNode(6).SetSolutionStepValue(KratosMultiphysics.NORMAL, normal)
        model_part.GetNode(7).SetSolutionStepValue(KratosMultiphysics.NORMAL, normal)
        model_part.GetNode(8).SetSolutionStepValue(KratosMultiphysics.NORMAL, normal)

        matrix_solution = KratosMultiphysics.Matrix()
        solution = exact_integration.TestGetExactIntegration(cond1, cond2, matrix_solution)

        # Debug
        #if (solution == True):
            #print("Integration accomplished", matrix_solution)

        self.assertTrue(solution)
        self.assertAlmostEqual(matrix_solution[0, 0], -1.0 / 3.0)
        self.assertAlmostEqual(matrix_solution[0, 1], -2.0 / 3.0)
        self.assertAlmostEqual(matrix_solution[0, 2],  1.0 / 6.0)
        self.assertAlmostEqual(matrix_solution[1, 0],  2.0 / 3.0)
        self.assertAlmostEqual(matrix_solution[1, 1], -2.0 / 3.0)
        self.assertAlmostEqual(matrix_solution[1, 2],  1.0 / 6.0)
        self.assertAlmostEqual(matrix_solution[2, 0],  2.0 / 3.0)
        self.assertAlmostEqual(matrix_solution[2, 1],  1.0 / 3.0)
        self.assertAlmostEqual(matrix_solution[2, 2],  1.0 / 6.0)
        self.assertAlmostEqual(matrix_solution[3, 0], -2.0 / 3.0)
        self.assertAlmostEqual(matrix_solution[3, 1], -1.0 / 3.0)
        self.assertAlmostEqual(matrix_solution[3, 2],  1.0 / 6.0)
        self.assertAlmostEqual(matrix_solution[4, 0],  1.0 / 3.0)
        self.assertAlmostEqual(matrix_solution[4, 1],  2.0 / 3.0)
        self.assertAlmostEqual(matrix_solution[4, 2],  1.0 / 6.0)
        self.assertAlmostEqual(matrix_solution[5, 0], -2.0 / 3.0)
        self.assertAlmostEqual(matrix_solution[5, 1],  2.0 / 3.0)
        self.assertAlmostEqual(matrix_solution[5, 2],  1.0 / 6.0)

    def test_quadrilateral_exact_integration_2(self):
        current_model = KratosMultiphysics.Model()
        
        model_part = current_model.CreateModelPart("Main")
        model_part.SetBufferSize(3)
        model_part.AddProperties(KratosMultiphysics.Properties(1))
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)

        normal = KratosMultiphysics.Vector(3)

        # Quadrilateral 1
        model_part.CreateNewNode(1, 0.00, 0.00, 0.00)
        model_part.CreateNewNode(2, 1.00, 0.00, 0.00)
        model_part.CreateNewNode(3, 1.00, 1.00, 0.00)
        model_part.CreateNewNode(4, 0.00, 1.00, 0.00)

        cond1 = model_part.CreateNewCondition("SurfaceCondition3D4N", 1, [1, 2, 3, 4], model_part.GetProperties()[1])
        normal = cond1.GetGeometry().UnitNormal()
        cond1.SetValue(KratosMultiphysics.NORMAL, normal)
        model_part.GetNode(1).SetSolutionStepValue(KratosMultiphysics.NORMAL, normal)
        model_part.GetNode(2).SetSolutionStepValue(KratosMultiphysics.NORMAL, normal)
        model_part.GetNode(3).SetSolutionStepValue(KratosMultiphysics.NORMAL, normal)
        model_part.GetNode(4).SetSolutionStepValue(KratosMultiphysics.NORMAL, normal)

        # Creating the utility:
        exact_integration = KratosMultiphysics.ExactMortarIntegrationUtility3D4N()

        # Quadrilateral 2
        normal[2] = -1.0
        model_part.CreateNewNode(5, 0.50, 0.50, 0.01)
        model_part.CreateNewNode(6, 1.50, 0.50, 0.01)
        model_part.CreateNewNode(7, 1.50, 1.50, 0.01)
        model_part.CreateNewNode(8, 0.50, 1.50, 0.01)

        cond2 = model_part.CreateNewCondition("SurfaceCondition3D4N", 2, [5, 6, 7, 8], model_part.GetProperties()[1])
        normal = cond2.GetGeometry().UnitNormal()
        cond2.SetValue(KratosMultiphysics.NORMAL, normal)
        model_part.GetNode(5).SetSolutionStepValue(KratosMultiphysics.NORMAL, normal)
        model_part.GetNode(6).SetSolutionStepValue(KratosMultiphysics.NORMAL, normal)
        model_part.GetNode(7).SetSolutionStepValue(KratosMultiphysics.NORMAL, normal)
        model_part.GetNode(8).SetSolutionStepValue(KratosMultiphysics.NORMAL, normal)

        matrix_solution = KratosMultiphysics.Matrix()
        solution = exact_integration.TestGetExactIntegration(cond1, cond2, matrix_solution)

        # Debug
        #if (solution == True):
            #print("Integration accomplished", matrix_solution)

        self.assertTrue(solution)
        self.assertAlmostEqual(matrix_solution[0, 0], 2.0 / 6.0)
        self.assertAlmostEqual(matrix_solution[0, 1], 1.0 / 6.0)
        self.assertAlmostEqual(matrix_solution[0, 2], 1.0 / 24.0)
        self.assertAlmostEqual(matrix_solution[1, 0], 5.0 / 6.0)
        self.assertAlmostEqual(matrix_solution[1, 1], 1.0 / 6.0)
        self.assertAlmostEqual(matrix_solution[1, 2], 1.0 / 24.0)
        self.assertAlmostEqual(matrix_solution[2, 0], 5.0 / 6.0)
        self.assertAlmostEqual(matrix_solution[2, 1], 4.0 / 6.0)
        self.assertAlmostEqual(matrix_solution[2, 2], 1.0 / 24.0)
        self.assertAlmostEqual(matrix_solution[3, 0], 1.0 / 6.0)
        self.assertAlmostEqual(matrix_solution[3, 1], 2.0 / 6.0)
        self.assertAlmostEqual(matrix_solution[3, 2], 1.0 / 24.0)
        self.assertAlmostEqual(matrix_solution[4, 0], 4.0 / 6.0)
        self.assertAlmostEqual(matrix_solution[4, 1], 5.0 / 6.0)
        self.assertAlmostEqual(matrix_solution[4, 2], 1.0 / 24.0)
        self.assertAlmostEqual(matrix_solution[5, 0], 1.0 / 6.0)
        self.assertAlmostEqual(matrix_solution[5, 1], 5.0 / 6.0)
        self.assertAlmostEqual(matrix_solution[5, 2], 1.0 / 24.0)


if __name__ == '__main__':
    KratosUnittest.main()
