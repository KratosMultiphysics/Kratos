import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics as KM
import numpy as np

class TestGraphUtilities(KratosUnittest.TestCase):

    def test_graph_and_connected_components_non_consecutive(self):
        current_model = KM.Model()
        model_part = current_model.CreateModelPart("Main")
        model_part.AddNodalSolutionStepVariable(KM.TEMPERATURE)
        model_part.AddNodalSolutionStepVariable(KM.DISTANCE)

        #modelpart being created: note that it contains two connected components
        # 3 - 4 - 5 - 6       7 - 8
        # | / | / | / |       | / |
        # 100-2 - 9 - 101      110- 120
        model_part.CreateNewNode(100,0.0,0.0,0.0)
        model_part.CreateNewNode(2,1.0,0.0,0.0)
        model_part.CreateNewNode(3,0.0,1.0,0.0)
        model_part.CreateNewNode(4,1.0,1.0,0.0)
        model_part.CreateNewNode(5,2.0,1.0,0.0)
        model_part.CreateNewNode(6,3.0,1.0,0.0)
        model_part.CreateNewNode(7,4.0,1.0,0.0)
        model_part.CreateNewNode(8,5.0,1.0,0.0)
        model_part.CreateNewNode(9,2.0,0.0,0.0)
        model_part.CreateNewNode(101,3.0,0.0,0.0)
        model_part.CreateNewNode(110,4.0,0.0,0.0)
        model_part.CreateNewNode(120,5.0,0.0,0.0)

        for node in model_part.Nodes:
            node.AddDof(KM.TEMPERATURE)
            node.Free(KM.TEMPERATURE)

        model_part.CreateNewElement("Element2D3N", 1, [100,4,3], model_part.GetProperties()[1])
        model_part.CreateNewElement("Element2D3N", 2, [100,2,4], model_part.GetProperties()[1])
        model_part.CreateNewElement("Element2D3N", 3, [2,9,5], model_part.GetProperties()[1])
        model_part.CreateNewElement("Element2D3N", 4, [5,4,2], model_part.GetProperties()[1])
        model_part.CreateNewElement("Element2D3N", 5, [6,9,101], model_part.GetProperties()[1])
        model_part.CreateNewElement("Element2D3N", 6, [9,6,5], model_part.GetProperties()[1])
        model_part.CreateNewElement("Element2D3N", 7, [110,8,7], model_part.GetProperties()[1])
        model_part.CreateNewElement("Element2D3N", 8, [110,120,8], model_part.GetProperties()[1])

        #note that we first compute the graph
        row_indices,col_indices = KM.ModelPartGraphUtilities.ComputeCSRGraph(model_part)

        #and then the connected components. This is done this way so that the graph can be reused.
        #we also need to pass the Nodes list to this second function as nodes can be numbered non consecutively,
        #however the output is to a vector, which is then designed for use together with VariableUtils().SetSolutionStepValuesVector
        #on that list of nodes
        ncolors,colors = KM.ModelPartGraphUtilities.ComputeConnectedComponents(model_part.Nodes,row_indices,col_indices)

        self.assertEqual(ncolors,2)

        KM.VariableUtils().SetSolutionStepValuesVector(model_part.Nodes, KM.TEMPERATURE, colors, 0)

        #with color 0
        block1_ids = [100,2,3,4,5,6,9,101]
        for n_id in block1_ids:
            self.assertAlmostEqual(model_part.Nodes[n_id].GetSolutionStepValue(KM.TEMPERATURE,0), 0.0)

        #with color 1
        block1_ids = [7,8,110,120]
        for n_id in block1_ids:
            self.assertAlmostEqual(model_part.Nodes[n_id].GetSolutionStepValue(KM.TEMPERATURE,0), 1.0)


        ################## now let's define a distance fiel, so that it is positive everywhere except in nodes 2,5
        for node in model_part.Nodes:
            node.SetSolutionStepValue(KM.DISTANCE,0,1.0)
        model_part.Nodes[2].SetSolutionStepValue(KM.DISTANCE,0,-1.0)
        model_part.Nodes[5].SetSolutionStepValue(KM.DISTANCE,0,-1.0)

        #here we obtain the list of distances
        distances = KM.VariableUtils().GetSolutionStepValuesVector(model_part.Nodes, KM.DISTANCE, 0)

        #now set active=false when distance<0
        active_nodes = np.where(np.array(distances)>0, True, False)

        #compute the connected components taking into account active_nodes - note that the graph is the same as before
        ncolors,colors = KM.ModelPartGraphUtilities.ComputeConnectedComponentsWithActiveNodesCheck(model_part.Nodes,row_indices,col_indices,active_nodes)

        self.assertEqual(ncolors,3)

        KM.VariableUtils().SetSolutionStepValuesVector(model_part.Nodes, KM.TEMPERATURE, colors, 0)

        #with color 0
        block1_ids = [3,4,100]
        for n_id in block1_ids:
            self.assertAlmostEqual(model_part.Nodes[n_id].GetSolutionStepValue(KM.TEMPERATURE,0), 0.0)

        #with color 1
        block1_ids = [6,9,101]
        for n_id in block1_ids:
            self.assertAlmostEqual(model_part.Nodes[n_id].GetSolutionStepValue(KM.TEMPERATURE,0), 1.0)

        #with color 2
        block1_ids = [7,8,110,120]
        for n_id in block1_ids:
            self.assertAlmostEqual(model_part.Nodes[n_id].GetSolutionStepValue(KM.TEMPERATURE,0), 2.0)

        ############# check automatic application of fixity
        model_part.Nodes[2].Fix(KM.TEMPERATURE) #this one will be ignored as inactive
        model_part.Nodes[100].Fix(KM.TEMPERATURE)
        KM.ModelPartGraphUtilities.ApplyMinimalScalarFixity(model_part.Nodes,KM.TEMPERATURE, colors, ncolors)

        fixed_list = [node.IsFixed(KM.TEMPERATURE) for node in model_part.Nodes]

        expected_fixed_list = [True, False, False, False, True, True, False, False, True, False, False, False]
        self.assertListEqual(fixed_list, expected_fixed_list)


if __name__ == '__main__':
    KratosUnittest.main()