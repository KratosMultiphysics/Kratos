import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as kratos_unittest

class TestVertexMoprhingContainerVariableDataMapper(kratos_unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = Kratos.Model()
        cls.origin_model_part = cls.model.CreateModelPart("origin")
        cls.origin_model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
        cls.origin_model_part.AddNodalSolutionStepVariable(Kratos.DENSITY)
        cls.origin_model_part.AddNodalSolutionStepVariable(Kratos.VELOCITY)
        cls.origin_model_part.AddNodalSolutionStepVariable(Kratos.ACCELERATION)
        cls.origin_model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 3

        cls.destination_model_part = cls.model.CreateModelPart("destination")
        cls.destination_model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
        cls.destination_model_part.AddNodalSolutionStepVariable(Kratos.DENSITY)
        cls.destination_model_part.AddNodalSolutionStepVariable(Kratos.VELOCITY)
        cls.destination_model_part.AddNodalSolutionStepVariable(Kratos.ACCELERATION)
        cls.destination_model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 3

        origin_grid_size = 3
        for i in range(9):
            cls.origin_model_part.CreateNewNode(i+1, i // origin_grid_size, i % origin_grid_size, 0)

        destination_grid_size = 2
        for i in range(4):
            cls.destination_model_part.CreateNewNode(i+1, 1.0 / destination_grid_size + (i // destination_grid_size), 1.0 / destination_grid_size + (i % destination_grid_size), 0)

        node: Kratos.Node
        for node in cls.origin_model_part.Nodes:
            node.SetSolutionStepValue(Kratos.PRESSURE, node.Id)
            node.SetSolutionStepValue(Kratos.VELOCITY, Kratos.Array3([node.Id + 1, node.Id + 2, node.Id + 3]))

        cls.mapper_parameters = Kratos.Parameters("""{
            "filter_function_type"      : "linear",
            "filter_radius"             : 1.0,
            "max_nodes_in_filter_radius": 10000
        }""")

    def test_MapScalar(self):
        mapper = KratosOA.VertexMorphingNodalContainerVariableDataMapper(self.origin_model_part, self.destination_model_part, self.mapper_parameters.Clone())
        mapper.Update()

        origin_data = KratosOA.HistoricalContainerVariableDataHolder(self.origin_model_part)
        origin_data.ReadDataFromContainerVariable(Kratos.PRESSURE)

        destination_data = KratosOA.HistoricalContainerVariableDataHolder(self.destination_model_part)
        mapper.Map(origin_data, destination_data)
        destination_data.AssignDataToContainerVariable(Kratos.PRESSURE)

        neighbour_nodes_map = {
            1: [1, 2, 4, 5],
            2: [2, 3, 5, 6],
            3: [4, 5, 7, 8],
            4: [5, 6, 8, 9]
        }
        for dest_id, orig_ids in neighbour_nodes_map.items():
            dest_node: Kratos.Node = self.destination_model_part.GetNode(dest_id)
            v = 0.0
            for orig_id in orig_ids:
                orig_node: Kratos.Node = self.origin_model_part.GetNode(orig_id)
                v += orig_node.GetSolutionStepValue(Kratos.PRESSURE) * 0.25
            self.assertAlmostEqual(dest_node.GetSolutionStepValue(Kratos.PRESSURE), v, 12)

        inverse_neighbour_map = {
            1: [1],
            2: [1, 2],
            3: [2],
            4: [1, 3],
            5: [1, 2, 3, 4],
            6: [2, 4],
            7: [3],
            8: [3, 4],
            9: [4]
        }

        inverse_origin_data = KratosOA.HistoricalContainerVariableDataHolder(self.origin_model_part)
        mapper.InverseMap(inverse_origin_data, destination_data)
        inverse_origin_data.AssignDataToContainerVariable(Kratos.DENSITY)

        for orig_id, dest_ids in inverse_neighbour_map.items():
            orig_node: Kratos.Node = self.origin_model_part.GetNode(orig_id)
            v = 0.0
            for dest_id in dest_ids:
                dest_node: Kratos.Node = self.destination_model_part.GetNode(dest_id)
                v += dest_node.GetSolutionStepValue(Kratos.PRESSURE) * 0.25
            self.assertAlmostEqual(orig_node.GetSolutionStepValue(Kratos.DENSITY), v, 12)

    def test_MapArray(self):
        mapper = KratosOA.VertexMorphingNodalContainerVariableDataMapper(self.origin_model_part, self.destination_model_part, self.mapper_parameters.Clone())
        mapper.Update()

        origin_data = KratosOA.HistoricalContainerVariableDataHolder(self.origin_model_part)
        origin_data.ReadDataFromContainerVariable(Kratos.VELOCITY)

        destination_data = KratosOA.HistoricalContainerVariableDataHolder(self.destination_model_part)
        mapper.Map(origin_data, destination_data)
        destination_data.AssignDataToContainerVariable(Kratos.VELOCITY)

        neighbour_nodes_map = {
            1: [1, 2, 4, 5],
            2: [2, 3, 5, 6],
            3: [4, 5, 7, 8],
            4: [5, 6, 8, 9]
        }
        for dest_id, orig_ids in neighbour_nodes_map.items():
            dest_node: Kratos.Node = self.destination_model_part.GetNode(dest_id)
            v = Kratos.Array3([0, 0, 0])
            for orig_id in orig_ids:
                orig_node: Kratos.Node = self.origin_model_part.GetNode(orig_id)
                v += orig_node.GetSolutionStepValue(Kratos.VELOCITY) * 0.25
            self.assertVectorAlmostEqual(dest_node.GetSolutionStepValue(Kratos.VELOCITY), v, 12)

        inverse_neighbour_map = {
            1: [1],
            2: [1, 2],
            3: [2],
            4: [1, 3],
            5: [1, 2, 3, 4],
            6: [2, 4],
            7: [3],
            8: [3, 4],
            9: [4]
        }

        inverse_origin_data = KratosOA.HistoricalContainerVariableDataHolder(self.origin_model_part)
        mapper.InverseMap(inverse_origin_data, destination_data)
        inverse_origin_data.AssignDataToContainerVariable(Kratos.ACCELERATION)

        for orig_id, dest_ids in inverse_neighbour_map.items():
            orig_node: Kratos.Node = self.origin_model_part.GetNode(orig_id)
            v = Kratos.Array3([0, 0, 0])
            for dest_id in dest_ids:
                dest_node: Kratos.Node = self.destination_model_part.GetNode(dest_id)
                v += dest_node.GetSolutionStepValue(Kratos.VELOCITY) * 0.25
            self.assertVectorAlmostEqual(orig_node.GetSolutionStepValue(Kratos.ACCELERATION), v, 12)

if __name__ == "__main__":
    Kratos.Tester.SetVerbosity(Kratos.Tester.Verbosity.PROGRESS)  # TESTS_OUTPUTS
    kratos_unittest.main()