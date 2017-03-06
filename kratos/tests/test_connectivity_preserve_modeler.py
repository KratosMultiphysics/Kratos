from __future__ import print_function, absolute_import, division

import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics 

class TestConnectivityPreserveModeler(KratosUnittest.TestCase):

    def test_connectivity_preserve_modeler(self):
        model_part1 = KratosMultiphysics.ModelPart("Main")
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
        
        model_part1.CreateNewCondition("Condition2D2N", 2, [2,4], model_part1.GetProperties()[1])
        
        new_model_part = KratosMultiphysics.ModelPart("Other")
        modeler = KratosMultiphysics.ConnectivityPreserveModeler()
        modeler.GenerateModelPart(model_part1, new_model_part, "Element2D3N", "Condition2D2N")
        
        self.assertEqual(len(model_part1.Nodes) , len(new_model_part.Nodes))
        self.assertEqual(len(model_part1.Conditions) , len(new_model_part.Conditions))
        self.assertEqual(len(model_part1.Elements) , len(new_model_part.Elements))
        
        ##assign a value to an element in model_part1, the corresponding element in the other model part will not be changed
        model_part1.Elements[1].SetValue(KratosMultiphysics.DISTANCE, 1.0)
        self.assertEqual(model_part1.Elements[1].GetValue(KratosMultiphysics.DISTANCE) , 1.0)
        self.assertEqual(new_model_part.Elements[1].GetValue(KratosMultiphysics.DISTANCE) , 0.0)
    
        ##assign a value to an element in model_part1, the corresponding element in the other model part will not be changed
        model_part1.Conditions[1].SetValue(KratosMultiphysics.DISTANCE, 1.0)
        self.assertEqual(model_part1.Conditions[1].GetValue(KratosMultiphysics.DISTANCE) , 1.0)
        self.assertEqual(new_model_part.Conditions[1].GetValue(KratosMultiphysics.DISTANCE) , 0.0)

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
                
if __name__ == '__main__':
    KratosUnittest.main()
