import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as KratosUtils
from KratosMultiphysics.testing.utilities import ReadModelPart

class TestModelPartOperationUtilities(KratosUnittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")
        cls.model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 2
        with KratosUnittest.WorkFolderScope("test_files/mdpa_files", __file__):
            ReadModelPart("model_part_operations", cls.model_part)

        cls.list_of_model_parts = [
            "test.Parts_Solid_sub_domain_1",
            "test.Parts_Solid_sub_domain_2",
            "test.Parts_Solid_sub_domain_3",
            "test.Parts_Solid_sub_domain_4"]

        cls.model_parts_list = [cls.model[model_part_name] for model_part_name in cls.list_of_model_parts]

    @classmethod
    def tearDownClass(cls) -> None:
        with KratosUnittest.WorkFolderScope("test_files/mdpa_files", __file__):
            # Clean up temporary files
            KratosUtils.DeleteFileIfExisting("model_part_operations.time")

    @staticmethod
    def __GenerateUniqueIdsList(container) -> 'list[int]':
        entity_ids = []
        for entity in container:
            entity_ids.append(entity.Id)
        return entity_ids

    def __CheckNeighbours(self, node_ids: 'list[int]', element_ids: 'list[int]', merged_model_part: Kratos.ModelPart, merged_model_part_with_neighbours: Kratos.ModelPart):
        merged_with_neighbours_node_ids = list(set(TestModelPartOperationUtilities.__GenerateUniqueIdsList(merged_model_part_with_neighbours.Nodes)))
        merged_with_neighbours_element_ids = list(set(TestModelPartOperationUtilities.__GenerateUniqueIdsList(merged_model_part_with_neighbours.Elements)))

        # check for neighbour nodes
        process = Kratos.FindGlobalNodalNeighboursProcess(self.model_part)
        process.Execute()
        node_node_ids_map = process.GetNeighbourIds(self.model_part.Nodes)

        for node in merged_model_part.Nodes:
            node_ids.extend(node_node_ids_map[node.Id])
        node_ids = list(set(node_ids))
        self.assertEqual(sorted(node_ids), sorted(merged_with_neighbours_node_ids))

        # check for neigghbour elements
        process = Kratos.FindGlobalNodalElementalNeighboursProcess(self.model_part)
        process.Execute()
        node_element_ids_map = process.GetNeighbourIds(self.model_part.Nodes)
        for node in merged_model_part.Nodes:
            element_ids.extend(node_element_ids_map[node.Id])
        element_ids = list(set(element_ids))
        self.assertEqual(sorted(element_ids), sorted(merged_with_neighbours_element_ids))

    @KratosUnittest.skipIf(Kratos.IsDistributedRun(), "only the test does not support MPI")
    def test_CheckValidityOfModelPartsForOperations(self):
        # check against empty model part
        empty_model_part = self.model.CreateModelPart("empty_test")
        with self.assertRaises(RuntimeError):
            Kratos.ModelPartOperationUtilities.CheckValidityOfModelPartsForOperations(empty_model_part, self.model_parts_list, True)
        self.assertFalse(Kratos.ModelPartOperationUtilities.CheckValidityOfModelPartsForOperations(empty_model_part, self.model_parts_list, False))

        # check against same root model part
        self.assertTrue(Kratos.ModelPartOperationUtilities.CheckValidityOfModelPartsForOperations(self.model_part, self.model_parts_list, True))

        # check against different root model part
        temp_model_part = self.model.CreateModelPart("temp_test")
        Kratos.ConnectivityPreserveModeler().GenerateModelPart(self.model_part, temp_model_part, "Element2D4N", "LineCondition2D2N")
        self.assertTrue(Kratos.ModelPartOperationUtilities.CheckValidityOfModelPartsForOperations(temp_model_part, self.model_parts_list, True))

    @KratosUnittest.skipIf(Kratos.IsDistributedRun(), "only the test does not support MPI")
    def test_Order(self):
        merged_model_part = self.model_part.CreateSubModelPart("merge_order")
        Kratos.ModelPartOperationUtilities.Union(merged_model_part, self.model_part, [self.model_part], False)

        origin_node_ids_list = [entity.Id for entity in self.model_part.Nodes]
        origin_condition_ids_list = [entity.Id for entity in self.model_part.Conditions]
        origin_element_ids_list = [entity.Id for entity in self.model_part.Elements]

        merged_node_ids_list = [entity.Id for entity in merged_model_part.Nodes]
        merged_condition_ids_list = [entity.Id for entity in merged_model_part.Conditions]
        merged_element_ids_list = [entity.Id for entity in merged_model_part.Elements]

        self.assertEqual(origin_node_ids_list, merged_node_ids_list)
        self.assertEqual(origin_condition_ids_list, merged_condition_ids_list)
        self.assertEqual(origin_element_ids_list, merged_element_ids_list)

    @KratosUnittest.skipIf(Kratos.IsDistributedRun(), "only the test does not support MPI")
    def test_Union(self):
        merged_model_part = self.model_part.CreateSubModelPart("merge_1")
        Kratos.ModelPartOperationUtilities.Union(merged_model_part, self.model_part, self.model_parts_list, False)

        # check
        node_ids = []
        condition_ids = []
        element_ids = []
        for model_part in self.model_parts_list:
            node_ids.extend(TestModelPartOperationUtilities.__GenerateUniqueIdsList(model_part.Nodes))
            condition_ids.extend(TestModelPartOperationUtilities.__GenerateUniqueIdsList(model_part.Conditions))
            element_ids.extend(TestModelPartOperationUtilities.__GenerateUniqueIdsList(model_part.Elements))

        node_ids = list(set(node_ids))
        condition_ids = list(set(condition_ids))
        element_ids = list(set(element_ids))

        merged_node_ids = list(set(TestModelPartOperationUtilities.__GenerateUniqueIdsList(merged_model_part.Nodes)))
        merged_condition_ids = list(set(TestModelPartOperationUtilities.__GenerateUniqueIdsList(merged_model_part.Conditions)))
        merged_element_ids = list(set(TestModelPartOperationUtilities.__GenerateUniqueIdsList(merged_model_part.Elements)))

        self.assertEqual(sorted(node_ids), sorted(merged_node_ids))
        self.assertEqual(sorted(condition_ids), sorted(merged_condition_ids))
        self.assertEqual(sorted(element_ids), sorted(merged_element_ids))

        merged_with_neighbours_model_part = self.model_part.CreateSubModelPart("merge_2")
        Kratos.ModelPartOperationUtilities.Union(merged_with_neighbours_model_part, self.model_part, self.model_parts_list, True)
        self.__CheckNeighbours(node_ids, element_ids, merged_model_part, merged_with_neighbours_model_part)

    @KratosUnittest.skipIf(Kratos.IsDistributedRun(), "only the test does not support MPI")
    def test_Substract(self):
        substracted_model_part = self.model_part.CreateSubModelPart("substract_1")
        Kratos.ModelPartOperationUtilities.Substract(substracted_model_part, self.model_part, self.model_parts_list, False)

        # check
        node_ids = TestModelPartOperationUtilities.__GenerateUniqueIdsList(self.model_part.Nodes)
        condition_ids = TestModelPartOperationUtilities.__GenerateUniqueIdsList(self.model_part.Conditions)
        element_ids = TestModelPartOperationUtilities.__GenerateUniqueIdsList(self.model_part.Elements)

        nodes_to_be_removed = []
        conditions_to_be_removed = []
        elements_to_be_removed = []
        for model_part in self.model_parts_list:
            nodes_to_be_removed.extend(TestModelPartOperationUtilities.__GenerateUniqueIdsList(model_part.Nodes))
            conditions_to_be_removed.extend(TestModelPartOperationUtilities.__GenerateUniqueIdsList(model_part.Conditions))
            elements_to_be_removed.extend(TestModelPartOperationUtilities.__GenerateUniqueIdsList(model_part.Elements))

        node_ids = set(node_ids)
        condition_ids = set(condition_ids)
        element_ids = set(element_ids)

        nodes_to_be_removed = set(nodes_to_be_removed)
        conditions_to_be_removed = set(conditions_to_be_removed)
        elements_to_be_removed = set(elements_to_be_removed)

        condition_ids = list(condition_ids.difference(conditions_to_be_removed))
        element_ids = list(element_ids.difference(elements_to_be_removed))

        node_ids = node_ids.difference(nodes_to_be_removed)
        node_ids = list(node_ids)
        # now need to add back the boundary nodes
        for element_id in element_ids:
            node_ids.extend([node.Id for node in self.model_part.GetElement(element_id).GetGeometry()])
        node_ids = list(set(node_ids))

        merged_node_ids = list(set(TestModelPartOperationUtilities.__GenerateUniqueIdsList(substracted_model_part.Nodes)))
        merged_condition_ids = list(set(TestModelPartOperationUtilities.__GenerateUniqueIdsList(substracted_model_part.Conditions)))
        merged_element_ids = list(set(TestModelPartOperationUtilities.__GenerateUniqueIdsList(substracted_model_part.Elements)))

        self.assertEqual(sorted(node_ids), sorted(merged_node_ids))
        self.assertEqual(sorted(condition_ids), sorted(merged_condition_ids))
        self.assertEqual(sorted(element_ids), sorted(merged_element_ids))

        substract_with_neighbours_model_part = self.model_part.CreateSubModelPart("substract_2")
        Kratos.ModelPartOperationUtilities.Substract(substract_with_neighbours_model_part, self.model_part, self.model_parts_list, True)
        self.__CheckNeighbours(node_ids, element_ids, substracted_model_part, substract_with_neighbours_model_part)

    @KratosUnittest.skipIf(Kratos.IsDistributedRun(), "only the test does not support MPI")
    def test_Intersect(self):
        intersected_model_part = self.model_part.CreateSubModelPart("intersect_1")
        Kratos.ModelPartOperationUtilities.Intersect(intersected_model_part, self.model_part, self.model_parts_list, False)

        # check
        node_ids = TestModelPartOperationUtilities.__GenerateUniqueIdsList(self.model_parts_list[0].Nodes)
        condition_ids = TestModelPartOperationUtilities.__GenerateUniqueIdsList(self.model_parts_list[0].Conditions)
        element_ids = TestModelPartOperationUtilities.__GenerateUniqueIdsList(self.model_parts_list[0].Elements)
        for model_part in self.model_parts_list[1:]:
            node_ids = list(set(node_ids).intersection(TestModelPartOperationUtilities.__GenerateUniqueIdsList(model_part.Nodes)))
            condition_ids = list(set(condition_ids).intersection(TestModelPartOperationUtilities.__GenerateUniqueIdsList(model_part.Conditions)))
            element_ids = list(set(element_ids).intersection(TestModelPartOperationUtilities.__GenerateUniqueIdsList(model_part.Elements)))

        merged_node_ids = list(set(TestModelPartOperationUtilities.__GenerateUniqueIdsList(intersected_model_part.Nodes)))
        merged_condition_ids = list(set(TestModelPartOperationUtilities.__GenerateUniqueIdsList(intersected_model_part.Conditions)))
        merged_element_ids = list(set(TestModelPartOperationUtilities.__GenerateUniqueIdsList(intersected_model_part.Elements)))

        self.assertEqual(sorted(node_ids), sorted(merged_node_ids))
        self.assertEqual(sorted(condition_ids), sorted(merged_condition_ids))
        self.assertEqual(sorted(element_ids), sorted(merged_element_ids))

        intersect_with_neighbours_model_part = self.model_part.CreateSubModelPart("intersect_2")
        Kratos.ModelPartOperationUtilities.Intersect(intersect_with_neighbours_model_part, self.model_part, self.model_parts_list, True)
        self.__CheckNeighbours(node_ids, element_ids, intersected_model_part, intersect_with_neighbours_model_part)

    @KratosUnittest.skipIf(Kratos.IsDistributedRun(), "only the test does not support MPI")
    def test_HasIntersection(self):
        self.assertTrue(Kratos.ModelPartOperationUtilities.HasIntersection(self.model_parts_list))

    def test_Sum(self):
        # communicator check
        for node in self.model_part.Nodes:
            node.SetValue(Kratos.PRESSURE, node.Id)

        for condition in self.model_part.Conditions:
            condition.SetValue(Kratos.PRESSURE, condition.Id)

        for element in self.model_part.Elements:
            element.SetValue(Kratos.PRESSURE, element.Id)

        merged_model_part = self.model_part.GetSubModelPart("Parts_Solid_sub_domain_1").CreateSubModelPart("merge_3")
        Kratos.ModelPartOperationUtilities.Union(merged_model_part, self.model_part, self.model_parts_list, False)
        substract_model_part = self.model_part.GetSubModelPart("Parts_Solid_sub_domain_1").CreateSubModelPart("substract_3")
        Kratos.ModelPartOperationUtilities.Substract(substract_model_part, self.model_part, self.model_parts_list, False)
        intersect_model_part = self.model_part.GetSubModelPart("Parts_Solid_sub_domain_1").CreateSubModelPart("intersect_3")
        Kratos.ModelPartOperationUtilities.Intersect(intersect_model_part, self.model_part, self.model_parts_list, False)

        merged_sum = Kratos.VariableUtils().SumNonHistoricalNodeScalarVariable(Kratos.PRESSURE, merged_model_part)
        self.assertEqual(merged_sum, 13415.0)

        substracted_sum = Kratos.VariableUtils().SumNonHistoricalNodeScalarVariable(Kratos.PRESSURE, substract_model_part)
        self.assertEqual(substracted_sum, 34926.0)

        intersected_sum = Kratos.VariableUtils().SumNonHistoricalNodeScalarVariable(Kratos.PRESSURE, intersect_model_part)
        self.assertEqual(intersected_sum, 2894.0)

    @KratosUnittest.skipIf(Kratos.IsDistributedRun(), "only the test does not support MPI")
    def test_CheckParentModelPart(self):
        model = Kratos.Model()
        model_part = model.CreateModelPart("test")
        sub_model_part1 = model_part.CreateSubModelPart("sub_test")
        sub_sub_model_part1 = sub_model_part1.CreateSubModelPart("sub_sub_test1")
        sub_sub_model_part2 = sub_model_part1.CreateSubModelPart("sub_sub_test2")

        output_model_part = model_part.CreateSubModelPart("output")
        Kratos.ModelPartOperationUtilities.Intersect(output_model_part, model_part, [sub_sub_model_part1, sub_sub_model_part2], False)

        output_model_part_2 = model.CreateModelPart("test2")
        with self.assertRaises(RuntimeError):
            Kratos.ModelPartOperationUtilities.Intersect(output_model_part_2, model_part, [sub_sub_model_part1, sub_sub_model_part2], False)


if __name__ == '__main__':
    KratosUnittest.main()