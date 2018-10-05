from __future__ import print_function, absolute_import, division

import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics import *

import sys

class TestModelPart(KratosUnittest.TestCase):

    def setUp(self):
        if (sys.version_info < (3, 2)):
            self.assertRaisesRegex = self.assertRaisesRegexp

    def test_model_part_sub_model_parts(self):
        current_model = Model()

        model_part= current_model.CreateModelPart("Main")

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

    def test_variables_list(self):
        current_model = Model()

        model_part= current_model.CreateModelPart("Main")

        self.assertEqual(model_part.GetNodalSolutionStepDataSize(), 0)

        model_part.AddNodalSolutionStepVariable(TEMPERATURE)

        self.assertEqual(model_part.GetNodalSolutionStepDataSize(), 1)

        model_part.AddNodalSolutionStepVariable(DISPLACEMENT)

        self.assertEqual(model_part.GetNodalSolutionStepDataSize(), 4)

        model_part.CreateSubModelPart("Inlets")
        sub_model_part_1 = model_part.GetSubModelPart("Inlets")

        self.assertEqual(sub_model_part_1.GetNodalSolutionStepDataSize(), 4)

        model_part.AddNodalSolutionStepVariable(VELOCITY)

        self.assertEqual(model_part.GetNodalSolutionStepDataSize(), 7)
        self.assertEqual(sub_model_part_1.GetNodalSolutionStepDataSize(), 7)

        sub_model_part_1.AddNodalSolutionStepVariable(PRESSURE)

        self.assertEqual(model_part.GetNodalSolutionStepDataSize(), 8)
        self.assertEqual(sub_model_part_1.GetNodalSolutionStepDataSize(), 8)






    def test_model_part_nodes(self):

        current_model = Model()

        model_part= current_model.CreateModelPart("Main")

        self.assertEqual(model_part.NumberOfNodes(), 0)
        self.assertEqual(model_part.NumberOfNodes(0), 0)

        model_part.CreateNewNode(1, 1.00,0.00,0.00)

        self.assertEqual(model_part.NumberOfNodes(), 1)
        self.assertEqual(model_part.NumberOfNodes(0), 1)

        #trying to create a node with Id 1 and coordinates which are different from the ones of the existing node 1. Error
        with self.assertRaises(RuntimeError):
            model_part.CreateNewNode(1, 0.00,0.00,0.00)

        #here i try to create a node with Id 1 but the coordinates coincide with the ones of the existing node. EXISTING NODE is returned and no error is thrown
        model_part.CreateNewNode(1, 1.00,0.00,0.00)
        self.assertEqual(model_part.NumberOfNodes(), 1)
        self.assertEqual(model_part.GetNode(1).Id, 1)
        self.assertEqual(model_part.GetNode(1,0).X, 1.00)

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
        self.assertEqual(model_part.GetNode(1).X, 1.00) #here the coordinates are still  the same as the original node
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
        current_model = Model()

        model_part= current_model.CreateModelPart("Main")

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
        current_model = Model()

        model_part= current_model.CreateModelPart("Main")

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
        current_model = Model()

        model_part= current_model.CreateModelPart("Main")

        self.assertEqual(model_part.NumberOfElements(), 0)
        self.assertEqual(model_part.NumberOfElements(0), 0)

        model_part.CreateNewNode(1, 0.00,0.00,0.00)
        model_part.CreateNewNode(2, 1.00,0.00,0.00)
        model_part.CreateNewNode(3, 1.00,1.00,0.00)
        model_part.AddProperties(Properties(1))
        model_part.CreateNewElement("Element2D3N", 1, [1,2,3], model_part.GetProperties()[1])

        self.assertEqual(model_part.NumberOfElements(), 1)
        self.assertEqual(model_part.NumberOfElements(0), 1)

        #an error is thrown if i try to create an element with the same Id
        with self.assertRaises(RuntimeError):
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
        current_model = Model()

        model_part= current_model.CreateModelPart("Main")

        self.assertEqual(model_part.NumberOfConditions(), 0)
        self.assertEqual(model_part.NumberOfConditions(0), 0)

        model_part.CreateNewNode(1, 0.00,0.00,0.00)
        model_part.CreateNewNode(2, 1.00,0.00,0.00)
        model_part.CreateNewNode(3, 1.00,1.00,0.00)
        model_part.AddProperties(Properties(1))
        model_part.CreateNewCondition("Condition3D", 1, [1,2,3], model_part.GetProperties()[1])

        self.assertEqual(model_part.NumberOfConditions(), 1)
        self.assertEqual(model_part.NumberOfConditions(0), 1)

        with self.assertRaises(RuntimeError):
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
        current_model = Model()

        model_part= current_model.CreateModelPart("Main")
        model_part.AddNodalSolutionStepVariable(VELOCITY)

        model_part.CreateNewNode(1, 0.00,0.00,0.00)
        model_part.CreateNewNode(2, 1.00,0.00,0.00)
        model_part.CreateNewNode(3, 1.00,1.00,0.00)

        self.assertTrue(model_part.Nodes[1].SolutionStepsDataHas(VELOCITY))

    def test_modelpart_buffersize(self):
        current_model = Model()

        model_part= current_model.CreateModelPart("Main")
        model_part.SetBufferSize(3)

        model_part.CreateSubModelPart("submodel")
        submodel = model_part.GetSubModelPart("submodel")
        self.assertEqual(model_part.GetBufferSize(), submodel.GetBufferSize() )

    def test_add_node(self):
        current_model = Model()

        model_part1= current_model.CreateModelPart("Main")
        sub1 = model_part1.CreateSubModelPart("sub1")
        sub2 = model_part1.CreateSubModelPart("sub2")

        model_part2= current_model.CreateModelPart("Other")

        model_part1.CreateNewNode(1,0.0,0.1,0.2)
        model_part1.CreateNewNode(2,2.0,0.1,0.2)

        n1 = model_part2.CreateNewNode(1,1.0,1.1,0.2)
        n3 = model_part2.CreateNewNode(3,2.0,3.1,0.2)
        n4 = model_part2.CreateNewNode(4,2.0,3.1,10.2)

        #this should add node 3 to both sub1 and model_part1, but not to sub2
        sub1.AddNode( model_part2.Nodes[3], 0 )
        #self.assertTrue( n3.Id in sub1.Nodes )
        #self.assertTrue( n3.Id in model_part1.Nodes )
        #self.assertFalse( n3.Id in sub2.Nodes )
        self.assertTrue( n3 in sub1.Nodes )
        self.assertTrue( n3 in model_part1.Nodes )
        self.assertFalse( n3 in sub2.Nodes )


        ##next should throw an exception, since we try to add a node with Id1 which already exists
        with self.assertRaisesRegex(RuntimeError, "Error\: attempting to add pNewNode with Id \:1, unfortunately a \(different\) node with the same Id already exists\n"):
            sub2.AddNode( n1, 0 )

        #create two extra nodes in the model model_part2
        n5 = model_part2.CreateNewNode(5,2.0,3.1,0.2)
        n6 = model_part2.CreateNewNode(6,2.0,3.1,10.2)

        ### here we test adding a list of nodes at once
        #now add node 4 and 5 to the model_part1 by Id - here it fails since we did not yet add node 4
        with self.assertRaisesRegex(RuntimeError, "Error: while adding nodes to submodelpart, the node with Id 4 does not exist in the root model part"):
            sub1.AddNodes([4,5])

        model_part1.AddNode( n4, 0 )
        model_part1.AddNode( n5, 0 )


        sub1.AddNodes([4,5]) #now it works, since we already added the nodes
        self.assertTrue( n4.Id in sub1.Nodes )
        self.assertTrue( n5.Id in sub1.Nodes )
        self.assertFalse( n5.Id in sub2.Nodes )

    def test_add_condition(self):
        current_model = Model()

        model_part1= current_model.CreateModelPart("Main")
        sub1 = model_part1.CreateSubModelPart("sub1")
        sub2 = model_part1.CreateSubModelPart("sub2")

        model_part2= current_model.CreateModelPart("Other")

        model_part1.CreateNewNode(1,0.0,0.1,0.2)
        model_part1.CreateNewNode(2,2.0,0.1,0.2)

        n1 = model_part2.CreateNewNode(1,1.0,1.1,0.2)
        n3 = model_part2.CreateNewNode(3,2.0,3.1,0.2)
        n4 = model_part2.CreateNewNode(4,2.0,3.1,10.2)

        model_part1.CreateNewCondition("Condition2D", 1, [1,2], sub1.GetProperties()[1])
        model_part1.CreateNewCondition("Condition2D", 2, [1,2], sub1.GetProperties()[1])

        c1 = model_part2.CreateNewCondition("Condition3D", 1, [1,3,4], model_part2.GetProperties()[1])
        c3 = model_part2.CreateNewCondition("Condition3D", 3, [1,3,4], model_part2.GetProperties()[1])

        #this should add condition 3 to both sub1 and model_part1, but not to sub2
        sub1.AddCondition( model_part2.Conditions[3], 0 )
        self.assertTrue( c3.Id in sub1.Conditions )
        self.assertTrue( c3.Id in model_part1.Conditions )
        self.assertFalse( c3.Id in sub2.Conditions )

        ##next should throw an exception, since we try to add a condition with Id1 which already exists
        with self.assertRaisesRegex(RuntimeError, "Error\: attempting to add pNewCondition with Id \:1, unfortunately a \(different\) condition with the same Id already exists\n"):
            sub2.AddCondition( c1, 0 )

        ##now we add two conditions at once
        c4 = model_part2.CreateNewCondition("Condition3D", 4, [1,3,4], model_part2.GetProperties()[1])
        c5 = model_part2.CreateNewCondition("Condition3D", 5, [1,3,4], model_part2.GetProperties()[1])

        ### here we test adding a list of conditions at once
        #now add node 4 and 5 to the model_part1 by Id - here it fails since we did not yet add node 4
        with self.assertRaisesRegex(RuntimeError, "Error: the condition with Id 4 does not exist in the root model part"):
            sub1.AddConditions([4,5])

        model_part1.AddCondition( c4, 0 )
        model_part1.AddCondition( c5, 0 )


        sub1.AddConditions([4,5]) #now it works, since we already added the nodes
        self.assertTrue( c4.Id in sub1.Conditions )
        self.assertTrue( c5.Id in sub1.Conditions )
        self.assertFalse( c5.Id in sub2.Conditions )


    def test_add_element(self):
        current_model = Model()

        model_part1= current_model.CreateModelPart("Main")
        sub1 = model_part1.CreateSubModelPart("sub1")
        sub2 = model_part1.CreateSubModelPart("sub2")

        model_part2= current_model.CreateModelPart("Other")

        model_part1.CreateNewNode(1,0.0,0.1,0.2)
        model_part1.CreateNewNode(2,2.0,0.1,0.2)

        n1 = model_part2.CreateNewNode(1,1.0,1.1,0.2)
        n3 = model_part2.CreateNewNode(3,2.0,3.1,0.2)
        n4 = model_part2.CreateNewNode(4,2.0,3.1,10.2)

        model_part1.CreateNewElement("Element2D2N", 1, [1,2], sub1.GetProperties()[1])
        model_part1.CreateNewElement("Element2D2N", 2, [1,2], sub1.GetProperties()[1])

        c1 = model_part2.CreateNewElement("Element2D2N", 1, [3,4], model_part2.GetProperties()[1])
        c3 = model_part2.CreateNewElement("Element2D2N", 3, [3,4], model_part2.GetProperties()[1])

        #this should add condition 3 to both sub1 and model_part1, but not to sub2
        sub1.AddElement( model_part2.Elements[3], 0 )
        self.assertTrue( c3.Id in sub1.Elements )
        self.assertTrue( c3.Id in model_part1.Elements )
        self.assertFalse( c3.Id in sub2.Elements )

        ##next should throw an exception, since we try to add a node with Id1 which already exists
        with self.assertRaisesRegex(RuntimeError, "Error\: attempting to add pNewElement with Id \:1, unfortunately a \(different\) element with the same Id already exists\n"):
            sub2.AddElement( c1, 0 )

        e4 = model_part2.CreateNewElement("Element2D2N", 4, [1,3], model_part2.GetProperties()[1])
        e5 = model_part2.CreateNewElement("Element2D2N", 5, [1,3], model_part2.GetProperties()[1])


       ### here we test adding a list of elements at once
        #now add node 4 and 5 to the model_part1 by Id - here it fails since we did not yet add node 4
        with self.assertRaisesRegex(RuntimeError, "Error: the element with Id 4 does not exist in the root model part"):
            sub1.AddElements([4,5])

        model_part1.AddElement( e4, 0 )
        model_part1.AddElement( e5, 0 )


        sub1.AddElements([4,5]) #now it works, since we already added the nodes
        self.assertTrue( e4.Id in sub1.Elements )
        self.assertTrue( e5.Id in sub1.Elements )
        self.assertFalse( e5.Id in sub2.Elements )

    def test_model_part_iterators(self):
        current_model = Model()

        model_part1= current_model.CreateModelPart("Main")
        sub1 = model_part1.CreateSubModelPart("sub1")
        sub2 = model_part1.CreateSubModelPart("sub2")

        subsub1 = sub1.CreateSubModelPart("subsub1")

        names = set(["sub1","sub2"])

        counter = 0

        for subpart in model_part1.SubModelParts:
            part_name = subpart.Name
            if part_name in names:
                counter+=1

            if(subpart.Name == "sub1"):
                for subsubpart in subpart.SubModelParts:
                    self.assertEqual(subsubpart.Name,"subsub1")
        self.assertEqual(counter, 2)

    def test_model_part_has_solution_step_variable(self):
        current_model = Model()
        model_part = current_model.CreateModelPart("Main")
        model_part.AddNodalSolutionStepVariable(VELOCITY)

        self.assertTrue(model_part.HasNodalSolutionStepVariable(VELOCITY))
        self.assertFalse(model_part.HasNodalSolutionStepVariable(PRESSURE))

    def test_model_part_master_slave_constraint(self):
        current_model = Model()
        model_part = current_model.CreateModelPart("Main")
        model_part.AddNodalSolutionStepVariable(PRESSURE)
        n1 = model_part.CreateNewNode(1, 1.0,1.1,0.2)
        n2 = model_part.CreateNewNode(2, 2.0,3.1,0.2)
        VariableUtils().AddDof(PRESSURE,model_part)

        c1 = MasterSlaveConstraint(10)
        model_part.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", 1, n1, PRESSURE, n2, PRESSURE, 0.5, 0.0)

        model_part.AddMasterSlaveConstraint(c1)

        consts = model_part.GetMasterSlaveConstraints()

        self.assertTrue(len(consts) == 2)

        self.assertTrue( c1.Id in model_part.MasterSlaveConstraints )

        #now try to add to submodelparts
        sub1 = model_part.CreateSubModelPart("sub1")
        sub2 = model_part.CreateSubModelPart("sub2")
        subsub1 = sub1.CreateSubModelPart("subsub1")

        ss1 = subsub1.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", 2, n1, PRESSURE, n2, PRESSURE, 0.5, 0.0)

        self.assertTrue(ss1 in subsub1.MasterSlaveConstraints)
        self.assertTrue(ss1 in sub1.MasterSlaveConstraints)
        self.assertTrue(ss1 in model_part.MasterSlaveConstraints)
        self.assertFalse(ss1 in sub2.MasterSlaveConstraints)

        sub1.RemoveMasterSlaveConstraint(ss1)
        self.assertFalse(ss1 in subsub1.MasterSlaveConstraints)
        self.assertFalse(ss1 in sub1.MasterSlaveConstraints)
        self.assertTrue(ss1 in model_part.MasterSlaveConstraints)

        subsub1.RemoveMasterSlaveConstraintFromAllLevels(ss1)

        self.assertFalse(ss1 in model_part.MasterSlaveConstraints)

    def test_no_constructor(self):
        with self.assertRaisesRegex(TypeError, "Kratos.ModelPart: No constructor defined!"):
            ModelPart()


if __name__ == '__main__':
    KratosUnittest.main()
