import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics

class TestModelPartCombinationUtilities(KratosUnittest.TestCase):

    def test_combine_model_parts(self):
        current_model = KratosMultiphysics.Model()

        # First model part
        model_part_1= current_model.CreateModelPart("Main1")

        model_part_1.CreateNewNode(1, 0.00,0.00,0.00)
        model_part_1.CreateNewNode(2, 1.00,0.00,0.00)
        model_part_1.CreateNewNode(3, 1.00,1.00,0.00)
        model_part_1.AddProperties(KratosMultiphysics.Properties(1))
        model_part_1.CreateNewElement("Element2D3N", 1, [1,2,3], model_part_1.GetProperties()[1])
        model_part_1.CreateNewElement("Element2D3N", 2, [1,2,3], model_part_1.GetProperties()[1])

        model_part_1.CreateSubModelPart("Inlets1")
        inlets_model_part_1 = model_part_1.GetSubModelPart("Inlets1")
        inlets_model_part_1.CreateNewNode(4, 0.00,0.00,0.00)
        inlets_model_part_1.CreateNewNode(5, 1.00,0.00,0.00)
        inlets_model_part_1.CreateNewNode(6, 1.00,1.00,0.00)
        inlets_model_part_1.CreateNewElement("Element2D3N", 3, [4,5,6], model_part_1.GetProperties()[1])

        # Second model part
        model_part_2= current_model.CreateModelPart("Main2")

        model_part_2.CreateNewNode(1, 0.00,0.00,0.00)
        model_part_2.CreateNewNode(2, 1.00,0.00,0.00)
        model_part_2.CreateNewNode(3, 1.00,1.00,0.00)
        model_part_2.AddProperties(KratosMultiphysics.Properties(1))
        model_part_2.CreateNewCondition("SurfaceCondition3D3N", 1, [1,2,3], model_part_2.GetProperties()[1])
        model_part_2.CreateNewCondition("SurfaceCondition3D3N", 2, [1,2,3], model_part_2.GetProperties()[1])

        model_part_2.CreateSubModelPart("Inlets2")
        inlets_model_part_2 = model_part_2.GetSubModelPart("Inlets2")
        inlets_model_part_2.CreateNewNode(4, 0.00,0.00,0.00)
        inlets_model_part_2.CreateNewNode(5, 1.00,0.00,0.00)
        inlets_model_part_2.CreateNewNode(6, 1.00,1.00,0.00)
        inlets_model_part_2.CreateNewCondition("SurfaceCondition3D3N", 3, [4,5,6], model_part_2.GetProperties()[1])

        inlets_model_part_2.CreateSubModelPart("Inlet1")
        inlets_model_part_2.CreateSubModelPart("Inlet2")
        inlet2_model_part_2 = inlets_model_part_2.GetSubModelPart("Inlet2")
        inlet2_model_part_2.CreateNewNode(7, 0.00,0.00,0.00)
        inlet2_model_part_2.CreateNewNode(8, 1.00,0.00,0.00)
        inlet2_model_part_2.CreateNewNode(9, 1.00,1.00,0.00)
        inlet2_model_part_2.CreateNewCondition("SurfaceCondition3D3N", 4, [7,8,9], model_part_2.GetProperties()[1])

        # Combine model parts
        param = KratosMultiphysics.Parameters("""{
          "model_parts_list"         : ["Main1", "Main2"],
          "combined_model_part_name" : "CombinedModelParts",
          "echo_level"               : 0
        }""")
        combined_model_part = KratosMultiphysics.ModelPartCombinationUtilities(current_model).CombineModelParts(param)

        # Checks
        self.assertEqual(combined_model_part.Name, "CombinedModelParts")
        self.assertEqual(combined_model_part.NumberOfProperties(), 1)
        self.assertEqual(combined_model_part.NumberOfNodes(), 15)
        self.assertEqual(combined_model_part.NumberOfElements(), 3)
        self.assertEqual(combined_model_part.NumberOfConditions(), 4)
        self.assertTrue(combined_model_part.HasSubModelPart("Inlets1"))
        self.assertTrue(combined_model_part.HasSubModelPart("Inlets2"))
        self.assertTrue(combined_model_part.GetSubModelPart("Inlets2").HasSubModelPart("Inlet1"))
        self.assertTrue(combined_model_part.GetSubModelPart("Inlets2").HasSubModelPart("Inlet2"))

    @KratosUnittest.expectedFailure
    def test_combine_model_parts_fail(self):
        current_model = KratosMultiphysics.Model()

        # First model part
        model_part_1= current_model.CreateModelPart("Main1")

        model_part_1.CreateNewNode(1, 0.00,0.00,0.00)
        model_part_1.CreateNewNode(2, 1.00,0.00,0.00)
        model_part_1.CreateNewNode(3, 1.00,1.00,0.00)
        model_part_1.AddProperties(KratosMultiphysics.Properties(1))
        model_part_1.CreateNewElement("Element2D3N", 1, [1,2,3], model_part_1.GetProperties()[1])
        model_part_1.CreateNewElement("Element2D3N", 2, [1,2,3], model_part_1.GetProperties()[1])

        model_part_1.CreateSubModelPart("Inlets")
        model_part_1.CreateSubModelPart("Outlet")

        # Second model part
        model_part_2= current_model.CreateModelPart("Main2")

        model_part_2.CreateNewNode(1, 0.00,0.00,0.00)
        model_part_2.CreateNewNode(2, 1.00,0.00,0.00)
        model_part_2.CreateNewNode(3, 1.00,1.00,0.00)
        model_part_2.AddProperties(KratosMultiphysics.Properties(2))
        model_part_2.CreateNewCondition("SurfaceCondition3D3N", 1, [1,2,3], model_part_2.GetProperties()[2])
        model_part_2.CreateNewCondition("SurfaceCondition3D3N", 2, [1,2,3], model_part_2.GetProperties()[2])

        model_part_2.CreateSubModelPart("Inlets")
        model_part_2.CreateSubModelPart("Outlet")

        # Combine model parts
        param = KratosMultiphysics.Parameters("""{
          "model_parts_list"         : ["Main1", "Main2"]
        }""")
        KratosMultiphysics.ModelPartCombinationUtilities(current_model).CombineModelParts(param)

        # TODO: ADD CHECKS

if __name__ == '__main__':
    KratosUnittest.main()
