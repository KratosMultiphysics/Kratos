from __future__ import print_function, absolute_import, division

import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics import *

class TestModelPart(KratosUnittest.TestCase):

    def test_model_part_sub_model_parts(self):
        model_part = ModelPart("Main")

        self.assertEqual(model_part.NumberOfSubModelParts(), 0)

        model_part.CreateSubModelPart("Inlets")

        self.assertTrue(model_part.HasSubModelPart("Inlets"))
        self.assertEqual(model_part.NumberOfSubModelParts(), 1)
        self.assertEqual(model_part.GetSubModelPart("Inlets").Name, "Inlets")

        model_part.CreateSubModelPart("Temp")
        model_part.CreateSubModelPart("Outlet")

        self.assertTrue(model_part.HasSubModelPart("Temp"))
        self.assertTrue(model_part.HasSubModelPart("Outlet"))
        self.assertEqual(model_part.NumberOfSubModelParts(), 3)
        self.assertEqual(model_part.GetSubModelPart("Inlets").Name, "Inlets")
        self.assertEqual(model_part.GetSubModelPart("Outlet").Name, "Outlet")

        sub_model_part_1 = model_part.GetSubModelPart("Inlets")
        sub_model_part_1.CreateSubModelPart("Inlet1")
        sub_model_part_1.CreateSubModelPart("Inlet2")

        self.assertEqual(model_part.NumberOfSubModelParts(), 3)
        self.assertEqual(model_part.GetSubModelPart("Inlets").Name, "Inlets")
        self.assertEqual(model_part.GetSubModelPart("Outlet").Name, "Outlet")

        #print ("Removing Temp....")
        model_part.RemoveSubModelPart("Temp")
        #print ("Temp removed!")

        self.assertFalse(model_part.HasSubModelPart("Temp"))
        self.assertEqual(model_part.NumberOfSubModelParts(), 2)
        self.assertEqual(model_part.GetSubModelPart("Inlets").Name, "Inlets")
        self.assertEqual(model_part.GetSubModelPart("Outlet").Name, "Outlet")

        #print ("Removing Inlets....")
        model_part.RemoveSubModelPart(sub_model_part_1)
        #print ("Inlets removed!")

        self.assertFalse(model_part.HasSubModelPart("Inlets"))
        self.assertEqual(model_part.NumberOfSubModelParts(), 1)
        self.assertEqual(model_part.GetSubModelPart("Outlet").Name, "Outlet")

       #print ("Removing Outlet....")
        model_part.RemoveSubModelPart("Outlet")
        #print ("Outlet removed!")

        self.assertFalse(model_part.HasSubModelPart("Inlets"))
        self.assertEqual(model_part.NumberOfSubModelParts(), 0)
        #print (model_part)

    def test_model_part_nodes(self):
        model_part = ModelPart("Main")

        self.assertEqual(model_part.NumberOfNodes(), 0)
        self.assertEqual(model_part.NumberOfNodes(0), 0)

        model_part.CreateNewNode(1, 1.00,0.00,0.00)

        self.assertEqual(model_part.NumberOfNodes(), 1)
        self.assertEqual(model_part.NumberOfNodes(0), 1)

        model_part.CreateNewNode(1, 0.00,0.00,0.00) # This overwrites the previous one

        self.assertEqual(model_part.NumberOfNodes(), 1)
        self.assertEqual(model_part.GetNode(1).Id, 1)
        self.assertEqual(model_part.GetNode(1,0).X, 0.00)
        self.assertEqual(len(model_part.Nodes), 1)

        model_part.CreateNewNode(2000, 2.00,0.00,0.00)

        self.assertEqual(model_part.NumberOfNodes(), 2)
        self.assertEqual(model_part.GetNode(1).Id, 1)
        self.assertEqual(model_part.GetNode(2000).Id, 2000)
        self.assertEqual(model_part.GetNode(2000).X, 2.00)

        model_part.CreateNewNode(2, 2.00,0.00,0.00)

        self.assertEqual(model_part.NumberOfNodes(), 3)
        self.assertEqual(model_part.GetNode(1).Id, 1)
        self.assertEqual(model_part.GetNode(2).Id, 2)
        self.assertEqual(model_part.GetNode(1).X, 0.00)
        self.assertEqual(model_part.GetNode(2).X, 2.00)

        model_part.RemoveNode(2000)

        self.assertEqual(model_part.NumberOfNodes(), 2)

        model_part.CreateSubModelPart("Inlets")
        model_part.CreateSubModelPart("Temp")
        model_part.CreateSubModelPart("Outlet")
        inlets_model_part = model_part.GetSubModelPart("Inlets")
        inlets_model_part.CreateNewNode(3, 3.00,0.00,0.00)

        self.assertEqual(inlets_model_part.NumberOfNodes(), 1)
        self.assertEqual(inlets_model_part.GetNode(3).Id, 3)
        self.assertEqual(inlets_model_part.GetNode(3).X, 3.00)
        self.assertEqual(model_part.NumberOfNodes(), 3)
        self.assertEqual(model_part.GetNode(3).Id, 3)
        self.assertEqual(model_part.GetNode(3).X, 3.00)

        inlets_model_part.CreateSubModelPart("Inlet1")
        inlets_model_part.CreateSubModelPart("Inlet2")
        inlet2_model_part = inlets_model_part.GetSubModelPart("Inlet2")
        inlet2_model_part.CreateNewNode(4, 4.00,0.00,0.00)

        self.assertEqual(inlet2_model_part.NumberOfNodes(), 1)
        self.assertEqual(inlet2_model_part.GetNode(4).Id, 4)
        self.assertEqual(inlet2_model_part.GetNode(4).X, 4.00)
        self.assertEqual(inlets_model_part.NumberOfNodes(), 2)
        self.assertEqual(inlets_model_part.GetNode(4).Id, 4)
        self.assertEqual(inlets_model_part.GetNode(4).X, 4.00)
        self.assertEqual(model_part.NumberOfNodes(), 4)
        self.assertEqual(model_part.GetNode(4).Id, 4)

        inlets_model_part.CreateNewNode(5, 5.00,0.00,0.00)
        inlets_model_part.CreateNewNode(6, 6.00,0.00,0.00)
        inlet2_model_part.CreateNewNode(7, 7.00,0.00,0.00)
        inlet2_model_part.CreateNewNode(8, 8.00,0.00,0.00)

        self.assertEqual(inlet2_model_part.NumberOfNodes(), 3)
        self.assertEqual(inlets_model_part.NumberOfNodes(), 6)
        self.assertEqual(model_part.NumberOfNodes(), 8)
        self.assertEqual(model_part.GetNode(4).Id, 4)

        inlets_model_part.RemoveNode(4)

        self.assertEqual(inlet2_model_part.NumberOfNodes(), 2)
        self.assertEqual(inlets_model_part.NumberOfNodes(), 5)
        self.assertEqual(model_part.NumberOfNodes(), 8) # the parent model part remains intact
        self.assertEqual(model_part.GetNode(4).Id, 4)

        inlets_model_part.RemoveNodeFromAllLevels(4) # Remove from all levels will delete it from

        self.assertEqual(inlet2_model_part.NumberOfNodes(), 2)
        self.assertEqual(inlets_model_part.NumberOfNodes(), 5)
        self.assertEqual(model_part.NumberOfNodes(), 7)

    def test_model_part_tables(self):
        model_part = ModelPart("Main")

        self.assertEqual(model_part.NumberOfTables(), 0)

        table = PiecewiseLinearTable()
        table.AddRow(0.00,1.00)
        table.AddRow(1.00,2.00)
        table.AddRow(2.00,2.00)
        model_part.AddTable(1, table)

        self.assertEqual(model_part.NumberOfTables(), 1)
        self.assertEqual(model_part.GetTable(1).GetValue(4.00), 2.00)

        table.AddRow(3.00,3.00)

        self.assertEqual(model_part.GetTable(1).GetValue(4.00), 4.00)

        #model_part.RemoveTable(1)

        #self.assertEqual(model_part.NumberOfTables(), 0)

    def test_model_part_properties(self):
        model_part = ModelPart("Main")

        self.assertEqual(model_part.NumberOfProperties(), 0)
        self.assertEqual(model_part.NumberOfProperties(0), 0)

        model_part.AddProperties(Properties(1))

        self.assertEqual(model_part.NumberOfProperties(), 1)
        self.assertEqual(model_part.GetProperties()[1].Id, 1)
        self.assertEqual(model_part.GetProperties(0)[1].Id, 1)
        self.assertEqual(len(model_part.Properties), 1)

        model_part.AddProperties(Properties(2000))

        self.assertEqual(model_part.NumberOfProperties(), 2)
        self.assertEqual(model_part.GetProperties()[1].Id, 1)
        self.assertEqual(model_part.GetProperties()[2000].Id, 2000)

        model_part.AddProperties(Properties(2))

        self.assertEqual(model_part.NumberOfProperties(), 3)
        self.assertEqual(model_part.GetProperties()[1].Id, 1)
        self.assertEqual(model_part.GetProperties()[2].Id, 2)

        model_part.RemoveProperties(2000)

        self.assertEqual(model_part.NumberOfProperties(), 2)

        model_part.CreateSubModelPart("Inlets")
        model_part.CreateSubModelPart("Temp")
        model_part.CreateSubModelPart("Outlet")
        inlets_model_part = model_part.GetSubModelPart("Inlets")
        inlets_model_part.AddProperties(Properties(3))

        self.assertEqual(inlets_model_part.NumberOfProperties(), 1)
        self.assertEqual(inlets_model_part.GetProperties()[3].Id, 3)
        self.assertEqual(model_part.NumberOfProperties(), 3)
        self.assertEqual(model_part.GetProperties()[3].Id, 3)

        inlets_model_part.CreateSubModelPart("Inlet1")
        inlets_model_part.CreateSubModelPart("Inlet2")
        inlet2_model_part = inlets_model_part.GetSubModelPart("Inlet2")
        inlet2_model_part.AddProperties(Properties(4))

        self.assertEqual(inlet2_model_part.NumberOfProperties(), 1)
        self.assertEqual(inlet2_model_part.GetProperties()[4].Id, 4)
        self.assertEqual(inlets_model_part.NumberOfProperties(), 2)
        self.assertEqual(inlets_model_part.GetProperties()[4].Id, 4)
        self.assertEqual(model_part.NumberOfProperties(), 4)
        self.assertEqual(model_part.GetProperties()[4].Id, 4)

        inlets_model_part.AddProperties(Properties(5))
        inlets_model_part.AddProperties(Properties(6))
        inlet2_model_part.AddProperties(Properties(7))
        inlet2_model_part.AddProperties(Properties(8))

        self.assertEqual(inlet2_model_part.NumberOfProperties(), 3)
        self.assertEqual(inlets_model_part.NumberOfProperties(), 6)
        self.assertEqual(model_part.NumberOfProperties(), 8)
        self.assertEqual(model_part.GetProperties()[4].Id, 4)

        inlets_model_part.RemoveProperties(4)

        self.assertEqual(inlet2_model_part.NumberOfProperties(), 2)
        self.assertEqual(inlets_model_part.NumberOfProperties(), 5)
        self.assertEqual(model_part.NumberOfProperties(), 8) # the parent model part remains intact
        self.assertEqual(model_part.GetProperties()[4].Id, 4)

        inlets_model_part.RemovePropertiesFromAllLevels(4) # Remove from all levels will delete it from

        self.assertEqual(inlet2_model_part.NumberOfProperties(), 2)
        self.assertEqual(inlets_model_part.NumberOfProperties(), 5)
        self.assertEqual(model_part.NumberOfProperties(), 7)

    def test_model_part_elements(self):
        model_part = ModelPart("Main")

        self.assertEqual(model_part.NumberOfElements(), 0)
        self.assertEqual(model_part.NumberOfElements(0), 0)

        model_part.CreateNewNode(1, 0.00,0.00,0.00)
        model_part.CreateNewNode(2, 1.00,0.00,0.00)
        model_part.CreateNewNode(3, 1.00,1.00,0.00)
        model_part.AddProperties(Properties(1))
        model_part.CreateNewElement("Element2D3N", 1, [1,2,3], model_part.GetProperties()[1])

        self.assertEqual(model_part.NumberOfElements(), 1)
        self.assertEqual(model_part.NumberOfElements(0), 1)

        model_part.CreateNewElement("Element2D3N", 1, [1,2,3], model_part.GetProperties()[1])

        self.assertEqual(model_part.NumberOfElements(), 1)
        self.assertEqual(model_part.GetElement(1).Id, 1)
        self.assertEqual(model_part.GetElement(1,0).Id, 1)
        self.assertEqual(model_part.Elements[1].Id, 1)
        self.assertEqual(len(model_part.Elements), 1)

        model_part.CreateNewElement("Element2D3N", 2000, [1,2,3], model_part.GetProperties()[1])

        self.assertEqual(model_part.NumberOfElements(), 2)
        self.assertEqual(model_part.GetElement(1).Id, 1)
        self.assertEqual(model_part.GetElement(2000).Id, 2000)

        model_part.CreateNewElement("Element2D3N", 2, [1,2,3], model_part.GetProperties()[1])

        self.assertEqual(model_part.NumberOfElements(), 3)
        self.assertEqual(model_part.GetElement(1).Id, 1)
        self.assertEqual(model_part.GetElement(2).Id, 2)

        model_part.RemoveElement(2000)

        self.assertEqual(model_part.NumberOfElements(), 2)

        model_part.CreateSubModelPart("Inlets")
        model_part.CreateSubModelPart("Temp")
        model_part.CreateSubModelPart("Outlet")
        inlets_model_part = model_part.GetSubModelPart("Inlets")
        inlets_model_part.CreateNewNode(4, 0.00,0.00,0.00)
        inlets_model_part.CreateNewNode(5, 1.00,0.00,0.00)
        inlets_model_part.CreateNewNode(6, 1.00,1.00,0.00)
        inlets_model_part.CreateNewElement("Element2D3N", 3, [4,5,6], model_part.GetProperties()[1])

        self.assertEqual(inlets_model_part.NumberOfElements(), 1)
        self.assertEqual(inlets_model_part.GetElement(3).Id, 3)
        self.assertEqual(model_part.NumberOfElements(), 3)
        self.assertEqual(model_part.GetElement(3).Id, 3)

        inlets_model_part.CreateSubModelPart("Inlet1")
        inlets_model_part.CreateSubModelPart("Inlet2")
        inlet2_model_part = inlets_model_part.GetSubModelPart("Inlet2")
        inlet2_model_part.CreateNewNode(7, 0.00,0.00,0.00)
        inlet2_model_part.CreateNewNode(8, 1.00,0.00,0.00)
        inlet2_model_part.CreateNewNode(9, 1.00,1.00,0.00)
        inlet2_model_part.CreateNewElement("Element2D3N", 4, [7,8,9], model_part.GetProperties()[1])

        self.assertEqual(inlet2_model_part.NumberOfElements(), 1)
        self.assertEqual(inlet2_model_part.GetElement(4).Id, 4)
        self.assertEqual(inlets_model_part.NumberOfElements(), 2)
        self.assertEqual(inlets_model_part.GetElement(4).Id, 4)
        self.assertEqual(model_part.NumberOfElements(), 4)
        self.assertEqual(model_part.GetElement(4).Id, 4)

        inlets_model_part.CreateNewElement("Element2D3N", 5, [7,8,9], model_part.GetProperties()[1])
        inlets_model_part.CreateNewElement("Element2D3N", 6, [7,8,9], model_part.GetProperties()[1])
        inlet2_model_part.CreateNewElement("Element2D3N", 7, [7,8,9], model_part.GetProperties()[1])
        inlet2_model_part.CreateNewElement("Element2D3N", 8, [7,8,9], model_part.GetProperties()[1])

        self.assertEqual(inlet2_model_part.NumberOfElements(), 3)
        self.assertEqual(inlets_model_part.NumberOfElements(), 6)
        self.assertEqual(model_part.NumberOfElements(), 8)
        self.assertEqual(model_part.GetElement(4).Id, 4)

        inlets_model_part.RemoveElement(4)

        self.assertEqual(inlet2_model_part.NumberOfElements(), 2)
        self.assertEqual(inlets_model_part.NumberOfElements(), 5)
        self.assertEqual(model_part.NumberOfElements(), 8) # the parent model part remains intact
        self.assertEqual(model_part.GetElement(4).Id, 4)

        inlets_model_part.RemoveElementFromAllLevels(4) # Remove from all levels will delete it from

        self.assertEqual(inlet2_model_part.NumberOfElements(), 2)
        self.assertEqual(inlets_model_part.NumberOfElements(), 5)
        self.assertEqual(model_part.NumberOfElements(), 7)

    def test_model_part_conditions(self):
        model_part = ModelPart("Main")

        self.assertEqual(model_part.NumberOfConditions(), 0)
        self.assertEqual(model_part.NumberOfConditions(0), 0)

        model_part.CreateNewNode(1, 0.00,0.00,0.00)
        model_part.CreateNewNode(2, 1.00,0.00,0.00)
        model_part.CreateNewNode(3, 1.00,1.00,0.00)
        model_part.AddProperties(Properties(1))
        model_part.CreateNewCondition("Condition3D", 1, [1,2,3], model_part.GetProperties()[1])

        self.assertEqual(model_part.NumberOfConditions(), 1)
        self.assertEqual(model_part.NumberOfConditions(0), 1)

        model_part.CreateNewCondition("Condition3D", 1, [1,2,3], model_part.GetProperties()[1])

        self.assertEqual(model_part.NumberOfConditions(), 1)
        self.assertEqual(model_part.GetCondition(1).Id, 1)
        self.assertEqual(model_part.GetCondition(1,0).Id, 1)
        self.assertEqual(model_part.Conditions[1].Id, 1)
        self.assertEqual(len(model_part.Conditions), 1)

        model_part.CreateNewCondition("Condition2D", 2000, [2,3], model_part.GetProperties()[1])

        self.assertEqual(model_part.NumberOfConditions(), 2)
        self.assertEqual(model_part.GetCondition(1).Id, 1)
        self.assertEqual(model_part.GetCondition(2000).Id, 2000)

        model_part.CreateNewCondition("Condition3D", 2, [1,2,3], model_part.GetProperties()[1])

        self.assertEqual(model_part.NumberOfConditions(), 3)
        self.assertEqual(model_part.GetCondition(1).Id, 1)
        self.assertEqual(model_part.GetCondition(2).Id, 2)

        model_part.RemoveCondition(2000)

        self.assertEqual(model_part.NumberOfConditions(), 2)

        model_part.CreateSubModelPart("Inlets")
        model_part.CreateSubModelPart("Temp")
        model_part.CreateSubModelPart("Outlet")
        inlets_model_part = model_part.GetSubModelPart("Inlets")
        inlets_model_part.CreateNewNode(4, 0.00,0.00,0.00)
        inlets_model_part.CreateNewNode(5, 1.00,0.00,0.00)
        inlets_model_part.CreateNewNode(6, 1.00,1.00,0.00)
        inlets_model_part.CreateNewCondition("Condition3D", 3, [4,5,6], model_part.GetProperties()[1])

        self.assertEqual(inlets_model_part.NumberOfConditions(), 1)
        self.assertEqual(inlets_model_part.GetCondition(3).Id, 3)
        self.assertEqual(model_part.NumberOfConditions(), 3)
        self.assertEqual(model_part.GetCondition(3).Id, 3)

        inlets_model_part.CreateSubModelPart("Inlet1")
        inlets_model_part.CreateSubModelPart("Inlet2")
        inlet2_model_part = inlets_model_part.GetSubModelPart("Inlet2")
        inlet2_model_part.CreateNewNode(7, 0.00,0.00,0.00)
        inlet2_model_part.CreateNewNode(8, 1.00,0.00,0.00)
        inlet2_model_part.CreateNewNode(9, 1.00,1.00,0.00)
        inlet2_model_part.CreateNewCondition("Condition3D", 4, [7,8,9], model_part.GetProperties()[1])

        self.assertEqual(inlet2_model_part.NumberOfConditions(), 1)
        self.assertEqual(inlet2_model_part.GetCondition(4).Id, 4)
        self.assertEqual(inlets_model_part.NumberOfConditions(), 2)
        self.assertEqual(inlets_model_part.GetCondition(4).Id, 4)
        self.assertEqual(model_part.NumberOfConditions(), 4)
        self.assertEqual(model_part.GetCondition(4).Id, 4)

        inlets_model_part.CreateNewCondition("Condition3D", 5, [7,8,9], model_part.GetProperties()[1])
        inlets_model_part.CreateNewCondition("Condition3D", 6, [7,8,9], model_part.GetProperties()[1])
        inlet2_model_part.CreateNewCondition("Condition3D", 7, [7,8,9], model_part.GetProperties()[1])
        inlet2_model_part.CreateNewCondition("Condition3D", 8, [7,8,9], model_part.GetProperties()[1])

        self.assertEqual(inlet2_model_part.NumberOfConditions(), 3)
        self.assertEqual(inlets_model_part.NumberOfConditions(), 6)
        self.assertEqual(model_part.NumberOfConditions(), 8)
        self.assertEqual(model_part.GetCondition(4).Id, 4)

        inlets_model_part.RemoveCondition(4)

        self.assertEqual(inlet2_model_part.NumberOfConditions(), 2)
        self.assertEqual(inlets_model_part.NumberOfConditions(), 5)
        self.assertEqual(model_part.NumberOfConditions(), 8) # the parent model part remains intact
        self.assertEqual(model_part.GetCondition(4).Id, 4)

        inlets_model_part.RemoveConditionFromAllLevels(4) # Remove from all levels will delete it from

        self.assertEqual(inlet2_model_part.NumberOfConditions(), 2)
        self.assertEqual(inlets_model_part.NumberOfConditions(), 5)
        self.assertEqual(model_part.NumberOfConditions(), 7)

    def test_modelpart_variables_list(self):
        model_part = ModelPart("Main")
        model_part.AddNodalSolutionStepVariable(VELOCITY)
        model_part.AddNodalSolutionStepVariable(VELOCITIES)

        model_part.CreateNewNode(1, 0.00,0.00,0.00)
        model_part.CreateNewNode(2, 1.00,0.00,0.00)
        model_part.CreateNewNode(3, 1.00,1.00,0.00)

        self.assertTrue(model_part.Nodes[1].SolutionStepsDataHas(VELOCITY))
        self.assertTrue(model_part.Nodes[1].SolutionStepsDataHas(VELOCITIES))

    def test_modelpart_buffersize(self):
        model_part = ModelPart("Main")
        model_part.SetBufferSize(3)
        
        model_part.CreateSubModelPart("submodel")
        submodel = model_part.GetSubModelPart("submodel")
        self.assertEqual(model_part.GetBufferSize(), submodel.GetBufferSize() )

if __name__ == '__main__':
    KratosUnittest.main()
