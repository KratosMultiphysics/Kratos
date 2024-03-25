import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics

class TestConnectivityPreserveModeler(KratosUnittest.TestCase):

    def test_connectivity_preserve_modeler(self):
        current_model = KratosMultiphysics.Model()
        model_part1 = current_model.CreateModelPart("Main")
        model_part1.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)

        model_part1.CreateNewNode(1,0.0,0.1,0.2)
        model_part1.CreateNewNode(2,2.0,0.1,0.2)
        model_part1.CreateNewNode(3,1.0,1.1,0.2)
        model_part1.CreateNewNode(4,2.0,3.1,10.2)

        sub1 = model_part1.CreateSubModelPart("sub1")
        sub2 = model_part1.CreateSubModelPart("sub2")
        subsub1 = sub1.CreateSubModelPart("subsub1")
        subsub1.AddNodes([1,2])
        sub2.AddNodes([3])

        model_part1.CreateNewElement("Element2D3N", 1, [1,2,3], model_part1.GetProperties()[1])
        model_part1.CreateNewElement("Element2D3N", 2, [1,2,4], model_part1.GetProperties()[1])

        model_part1.CreateNewCondition("LineCondition2D2N", 2, [2,4], model_part1.GetProperties()[1])
        sub1.AddConditions([2])

        current_model = KratosMultiphysics.Model()
        new_model_part = current_model.CreateModelPart("Other")

        modeler = KratosMultiphysics.ConnectivityPreserveModeler()
        modeler.GenerateModelPart(model_part1, new_model_part, "Element2D3N", "LineCondition2D2N")

        self.assertEqual(len(model_part1.Nodes) , len(new_model_part.Nodes))
        self.assertEqual(len(model_part1.Conditions) , len(new_model_part.Conditions))
        self.assertEqual(len(model_part1.Elements) , len(new_model_part.Elements))

        ##assign a value to an element in model_part1, the corresponding element in the other model part
        # will change as they are pointing to the same geometry
        model_part1.Elements[1].SetValue(KratosMultiphysics.DISTANCE, 1.0)
        self.assertEqual(model_part1.Elements[1].GetValue(KratosMultiphysics.DISTANCE) , 1.0)
        self.assertEqual(new_model_part.Elements[1].GetValue(KratosMultiphysics.DISTANCE) , 1.0)

        ##assign a value to an condition in model_part1, the corresponding condition in the other model part
        # will change as they are pointing to the same geometry
        model_part1.Conditions[2].SetValue(KratosMultiphysics.DISTANCE, 1.0)
        self.assertEqual(model_part1.Conditions[2].GetValue(KratosMultiphysics.DISTANCE) , 1.0)
        self.assertEqual(new_model_part.Conditions[2].GetValue(KratosMultiphysics.DISTANCE) , 1.0)

        #assign a value to a node of new_model_part --> affects both model parts
        new_model_part.Nodes[1].SetSolutionStepValue(KratosMultiphysics.DISTANCE,0,2.0)
        self.assertEqual(model_part1.Nodes[1].GetSolutionStepValue(KratosMultiphysics.DISTANCE)   , 2.0)
        self.assertEqual(new_model_part.Nodes[1].GetSolutionStepValue(KratosMultiphysics.DISTANCE), 2.0)

        #test if submodelparts are created correctly
        for part in model_part1.SubModelParts:
            new_part = new_model_part.GetSubModelPart(part.Name)
            for node in part.Nodes:
                self.assertTrue( node.Id in new_part.Nodes)
            for cond in part.Conditions:
                self.assertTrue( cond.Id in new_part.Conditions)
            for elem in part.Elements:
                self.assertTrue( elem.Id in new_part.Elements)

        model_part1.GetSubModelPart("sub1").Conditions[2].SetValue(KratosMultiphysics.TEMPERATURE, 1234.0)
        self.assertEqual(model_part1.Conditions[2].GetValue(KratosMultiphysics.TEMPERATURE), 1234.0)
        self.assertEqual(model_part1.GetSubModelPart("sub1").Conditions[2].GetValue(KratosMultiphysics.TEMPERATURE), 1234.0)
        self.assertEqual(new_model_part.Conditions[2].GetValue(KratosMultiphysics.TEMPERATURE), 1234.0)
        self.assertEqual(new_model_part.GetSubModelPart("sub1").Conditions[2].GetValue(KratosMultiphysics.TEMPERATURE), 1234.0)

    def test_repeated_call(self):
        current_model = KratosMultiphysics.Model()
        model_part1 = current_model.CreateModelPart("Main")
        model_part1.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)

        model_part1.CreateNewNode(1,0.0,0.1,0.2)
        model_part1.CreateNewNode(2,2.0,0.1,0.2)
        model_part1.CreateNewNode(3,1.0,1.1,0.2)
        model_part1.CreateNewNode(4,2.0,3.1,10.2)

        sub1 = model_part1.CreateSubModelPart("sub1")
        subsub1 = sub1.CreateSubModelPart("subsub1")
        subsub1.AddNodes([1,2])

        model_part1.CreateNewElement("Element2D3N", 1, [1,2,3], model_part1.GetProperties()[1])
        model_part1.CreateNewElement("Element2D3N", 2, [1,2,4], model_part1.GetProperties()[1])

        model_part1.CreateNewCondition("LineCondition2D2N", 2, [2,4], model_part1.GetProperties()[1])
        model_part1.CreateNewCondition("LineCondition2D2N", 1, [1,2], model_part1.GetProperties()[1])
        sub1.AddConditions([2])

        current_model = KratosMultiphysics.Model()
        new_model_part = current_model.CreateModelPart("New1")
        new_model_part2 = current_model.CreateModelPart("New2")

        modeler = KratosMultiphysics.ConnectivityPreserveModeler()
        modeler.GenerateModelPart(model_part1, new_model_part, "Element2D3N", "LineCondition2D2N")
        self.assertEqual(len(model_part1.Nodes) , len(new_model_part.Nodes))
        self.assertEqual(len(model_part1.Conditions) , len(new_model_part.Conditions))
        self.assertEqual(len(model_part1.Elements) , len(new_model_part.Elements))

        modeler.GenerateModelPart(model_part1, new_model_part2, "Element2D3N", "LineCondition2D2N")
        self.assertEqual(len(model_part1.Nodes) , len(new_model_part2.Nodes))
        self.assertEqual(len(model_part1.Conditions) , len(new_model_part2.Conditions))
        self.assertEqual(len(model_part1.Elements) , len(new_model_part2.Elements))

        modeler.GenerateModelPart(model_part1, new_model_part, "Element2D3N", "LineCondition2D2N")
        self.assertEqual(len(model_part1.Nodes) , len(new_model_part.Nodes))
        self.assertEqual(len(model_part1.Conditions) , len(new_model_part.Conditions))
        self.assertEqual(len(model_part1.Elements) , len(new_model_part.Elements))

        modeler.GenerateModelPart(model_part1, new_model_part2, "Element2D3N", "LineCondition2D2N")
        self.assertEqual(len(model_part1.Nodes) , len(new_model_part2.Nodes))
        self.assertEqual(len(model_part1.Conditions) , len(new_model_part2.Conditions))
        self.assertEqual(len(model_part1.Elements) , len(new_model_part2.Elements))

        self.assertEqual(len(model_part1.Nodes) , len(new_model_part.Nodes))
        self.assertEqual(len(model_part1.Conditions) , len(new_model_part.Conditions))
        self.assertEqual(len(model_part1.Elements) , len(new_model_part.Elements))

    def test_variable_list_merging(self):
        current_model = KratosMultiphysics.Model()
        model_part1 = current_model.CreateModelPart("mp1")
        model_part1.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        model_part1.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)

        model_part2 = current_model.CreateModelPart("mp2")
        model_part2.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)
        model_part2.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)

        #check before merging variable lists
        self.assertTrue(model_part1.HasNodalSolutionStepVariable(KratosMultiphysics.DISTANCE))
        self.assertTrue(model_part1.HasNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT))
        self.assertFalse(model_part1.HasNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE))
        self.assertFalse(model_part1.HasNodalSolutionStepVariable(KratosMultiphysics.VELOCITY))
        self.assertFalse(model_part2.HasNodalSolutionStepVariable(KratosMultiphysics.DISTANCE))
        self.assertFalse(model_part2.HasNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT))
        self.assertTrue(model_part2.HasNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE))
        self.assertTrue(model_part2.HasNodalSolutionStepVariable(KratosMultiphysics.VELOCITY))

        KratosMultiphysics.MergeVariableListsUtility().Merge(model_part1, model_part2)

        #check after merging variable lists
        self.assertTrue(model_part1.HasNodalSolutionStepVariable(KratosMultiphysics.DISTANCE))
        self.assertTrue(model_part1.HasNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT))
        self.assertTrue(model_part1.HasNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE))
        self.assertTrue(model_part1.HasNodalSolutionStepVariable(KratosMultiphysics.VELOCITY))
        self.assertTrue(model_part2.HasNodalSolutionStepVariable(KratosMultiphysics.DISTANCE))
        self.assertTrue(model_part2.HasNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT))
        self.assertTrue(model_part2.HasNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE))
        self.assertTrue(model_part2.HasNodalSolutionStepVariable(KratosMultiphysics.VELOCITY))

    def test_element_only_copy(self):
        current_model = KratosMultiphysics.Model()
        model_part1 = current_model.CreateModelPart("Main")
        model_part1.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)

        model_part1.CreateNewNode(1,0.0,0.1,0.2)
        model_part1.CreateNewNode(2,2.0,0.1,0.2)
        model_part1.CreateNewNode(3,1.0,1.1,0.2)
        model_part1.CreateNewNode(4,2.0,3.1,10.2)

        sub1 = model_part1.CreateSubModelPart("sub1")

        model_part1.CreateNewElement("Element2D3N", 1, [1,2,3], model_part1.GetProperties()[1])
        model_part1.CreateNewElement("Element2D3N", 2, [1,2,4], model_part1.GetProperties()[1])

        model_part1.CreateNewCondition("LineCondition2D2N", 2, [2,4], model_part1.GetProperties()[1])
        sub1.AddConditions([2])

        element_model_part = current_model.CreateModelPart("ElementCopy")

        modeler = KratosMultiphysics.ConnectivityPreserveModeler()

        modeler.GenerateModelPart(model_part1, element_model_part, "Element2D3N")

        self.assertEqual(len(element_model_part.Nodes) , len(model_part1.Nodes))
        self.assertEqual(len(element_model_part.Elements) , len(model_part1.Elements))
        self.assertEqual(len(element_model_part.Conditions) , 0)

        element_sub1_copy = element_model_part.GetSubModelPart("sub1")

        self.assertEqual(len(element_sub1_copy.Nodes) , len(sub1.Nodes))
        self.assertEqual(len(element_sub1_copy.Elements) , len(sub1.Elements))
        self.assertEqual(len(element_sub1_copy.Conditions) , 0)

        condition_model_part = current_model.CreateModelPart("ConditionCopy")

        modeler.GenerateModelPart(model_part1, condition_model_part, "LineCondition2D2N")

        self.assertEqual(len(condition_model_part.Nodes) , len(model_part1.Nodes))
        self.assertEqual(len(condition_model_part.Elements) , 0)
        self.assertEqual(len(condition_model_part.Conditions) , len(model_part1.Conditions))

        condition_sub1_copy = condition_model_part.GetSubModelPart("sub1")

        self.assertEqual(len(condition_sub1_copy.Nodes) , len(sub1.Nodes))
        self.assertEqual(len(condition_sub1_copy.Elements) , 0)
        self.assertEqual(len(condition_sub1_copy.Conditions) , len(sub1.Conditions))


    def test_properties_shared(self):
        current_model = KratosMultiphysics.Model()
        model_part_1 = current_model.CreateModelPart("Main")
        model_part_2 = current_model.CreateModelPart("Destination")

        model_part_1.CreateNewNode(1, 0.0, 0.1, 0.0)
        model_part_1.CreateNewNode(2, 2.0, 0.1, 0.0)
        model_part_1.CreateNewNode(3, 1.0, 1.1, 0.0)

        properties_0 = model_part_1.CreateNewProperties(0)
        properties_1 = model_part_1.CreateNewProperties(1)
        properties_0.SetValue(KratosMultiphysics.DENSITY, 0.1)
        properties_1.SetValue(KratosMultiphysics.DENSITY, 1.0)

        model_part_1.CreateNewElement("Element2D3N", 1, [1,2,3], model_part_1.GetProperties()[1])

        KratosMultiphysics.ConnectivityPreserveModeler().GenerateModelPart(model_part_1, model_part_2, "Element2D3N")

        self.assertEqual(model_part_2.NumberOfProperties(), 2)
        self.assertEqual(model_part_2.GetProperties()[0].GetValue(KratosMultiphysics.DENSITY), 0.1)
        self.assertEqual(model_part_2.GetProperties()[1].GetValue(KratosMultiphysics.DENSITY), 1.0)

        model_part_2.RemoveProperties(0)
        model_part_2.GetProperties()[1].SetValue(KratosMultiphysics.DENSITY, 2.0)

        self.assertEqual(model_part_1.NumberOfProperties(), 2)
        self.assertEqual(model_part_1.GetProperties()[0].GetValue(KratosMultiphysics.DENSITY), 0.1)
        self.assertEqual(model_part_1.GetProperties()[1].GetValue(KratosMultiphysics.DENSITY), 2.0)


    def test_setup_with_json_parameters(self):
        current_model = KratosMultiphysics.Model()
        model_part1 = current_model.CreateModelPart("Main")
        model_part1.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)

        model_part1.CreateNewNode(1,0.0,0.1,0.2)
        model_part1.CreateNewNode(2,2.0,0.1,0.2)
        model_part1.CreateNewNode(3,1.0,1.1,0.2)
        model_part1.CreateNewNode(4,2.0,3.1,10.2)

        sub1 = model_part1.CreateSubModelPart("sub1")
        sub2 = model_part1.CreateSubModelPart("sub2")
        subsub1 = sub1.CreateSubModelPart("subsub1")
        subsub1.AddNodes([1,2])
        sub2.AddNodes([3])

        model_part1.CreateNewElement("Element2D3N", 1, [1,2,3], model_part1.GetProperties()[1])
        model_part1.CreateNewElement("Element2D3N", 2, [1,2,4], model_part1.GetProperties()[1])

        model_part1.CreateNewCondition("LineCondition2D2N", 2, [2,4], model_part1.GetProperties()[1])
        sub1.AddConditions([2])

        # Note: unlike other modes of operation, the json configuration requires both parts to be in the same Model
        new_model_part = current_model.CreateModelPart("Other")

        modeler = KratosMultiphysics.ConnectivityPreserveModeler().Create(
            current_model,
            KratosMultiphysics.Parameters('''{
                "origin_model_part_name": "Main",
                "destination_model_part_name": "Other",
                "reference_element": "Element2D3N",
                "reference_condition": "LineCondition2D2N"
            }''')
        )
        modeler.SetupModelPart()

        self.assertEqual(len(model_part1.Nodes) , len(new_model_part.Nodes))
        self.assertEqual(len(model_part1.Conditions) , len(new_model_part.Conditions))
        self.assertEqual(len(model_part1.Elements) , len(new_model_part.Elements))

        # In SetupModelPart, the destination model part can be created if it does not exist
        self.assertFalse(current_model.HasModelPart("CreatedModelPart"))
        modeler = KratosMultiphysics.ConnectivityPreserveModeler().Create(
            current_model,
            KratosMultiphysics.Parameters('''{
                "origin_model_part_name": "Main",
                "destination_model_part_name": "CreatedModelPart",
                "reference_element": "Element2D3N",
                "reference_condition": "LineCondition2D2N"
            }''')
        )
        modeler.SetupModelPart()

        created_model_part = current_model.GetModelPart("CreatedModelPart")
        self.assertEqual(len(model_part1.Nodes) , len(created_model_part.Nodes))
        self.assertEqual(len(model_part1.Conditions) , len(created_model_part.Conditions))
        self.assertEqual(len(model_part1.Elements) , len(created_model_part.Elements))

if __name__ == '__main__':
    KratosUnittest.main()
