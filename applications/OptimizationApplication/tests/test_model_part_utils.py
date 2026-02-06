import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as kratos_unittest
from KratosMultiphysics.testing.utilities import ReadModelPart

class TestOptAppModelPartUtils(kratos_unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")
        cls.model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 2
        ReadModelPart("model_part_utils_test/quads", cls.model_part)

        communicator: Kratos.Communicator = cls.model_part.GetCommunicator()

        cls.reference_local_node_ids = [node.Id for node in communicator.LocalMesh().Nodes]
        cls.reference_interface_node_ids = [node.Id for node in communicator.InterfaceMesh().Nodes]
        cls.reference_ghost_node_ids = [node.Id for node in communicator.GhostMesh().Nodes]

        cls.reference_local_condition_ids = [condition.Id for condition in communicator.LocalMesh().Conditions]
        cls.reference_interface_condition_ids = [condition.Id for condition in communicator.InterfaceMesh().Conditions]
        cls.reference_ghost_condition_ids = [condition.Id for condition in communicator.GhostMesh().Conditions]

        cls.reference_local_element_ids = [element.Id for element in communicator.LocalMesh().Elements]
        cls.reference_interface_element_ids = [element.Id for element in communicator.InterfaceMesh().Elements]
        cls.reference_ghost_element_ids = [element.Id for element in communicator.GhostMesh().Elements]

        cls.examined_model_parts = [
            cls.model["test.sensitivity_element_1"],
            cls.model["test.sensitivity_element_2"],
            cls.model["test.sensitivity_element_3"],
            cls.model["test.sensitivity_element_4"],
            cls.model["test.sensitivity_condition_1"]
        ]

        cls.reference_model_parts = [
            cls.model["test.evaluated_element_1"],
            cls.model["test.evaluated_element_2"],
            cls.model["test.evaluated_element_3"]
        ]

    @staticmethod
    def __GetAllIds(entities, communicator: Kratos.DataCommunicator):
        current_rank_node_ids = [entity.Id for entity in entities]
        max_nodes_size = communicator.MaxAll(len(current_rank_node_ids))
        for _ in range(max(max_nodes_size - len(current_rank_node_ids), 0)):
            current_rank_node_ids.append(0)
        result = list(set(communicator.AllGatherInts(current_rank_node_ids)))
        if len(result) > 0 and result[0] == 0:
            del result[0]
        return result

    def __CheckModelPart(self, model_part: Kratos.ModelPart):
        communicator: Kratos.Communicator = model_part.GetCommunicator()
        data_communicator: Kratos.DataCommunicator = communicator.GetDataCommunicator()

        self.assertEqual(data_communicator.IsDistributed(), self.model_part.GetCommunicator().GetDataCommunicator().IsDistributed())

        # check nodes
        for node in communicator.LocalMesh().Nodes:
            self.assertTrue(node.Id in self.reference_local_node_ids)

        for node in communicator.InterfaceMesh().Nodes:
            self.assertTrue(node.Id in self.reference_interface_node_ids)

        for node in communicator.GhostMesh().Nodes:
            self.assertTrue(node.Id in self.reference_ghost_node_ids)

        # check conditions
        for condition in communicator.LocalMesh().Conditions:
            self.assertTrue(condition.Id in self.reference_local_condition_ids)

        for condition in communicator.InterfaceMesh().Conditions:
            self.assertTrue(condition.Id in self.reference_interface_condition_ids)

        for condition in communicator.GhostMesh().Conditions:
            self.assertTrue(condition.Id in self.reference_ghost_condition_ids)

        # check elements
        for element in communicator.LocalMesh().Elements:
            self.assertTrue(element.Id in self.reference_local_element_ids)

        for element in communicator.InterfaceMesh().Elements:
            self.assertTrue(element.Id in self.reference_interface_element_ids)

        for element in communicator.GhostMesh().Elements:
            self.assertTrue(element.Id in self.reference_ghost_element_ids)

    def __CheckCommonEntities(self, model_parts: 'list[Kratos.ModelPart]', examined_model_parts: 'list[Kratos.ModelPart]', entity_container_getter):
        examined_entity_ids = []

        for model_part in examined_model_parts:
            examined_entity_ids.extend([entity.Id for entity in entity_container_getter(model_part)])

        examined_entity_ids = set(examined_entity_ids)

        for model_part in model_parts:
            current_entity_ids = sorted([entity.Id for entity in entity_container_getter(model_part)])

            reference_model_part = model_part.GetParentModelPart()
            reference_entity_ids = [entity.Id for entity in entity_container_getter(reference_model_part)]

            common_entity_ids = sorted(list(examined_entity_ids.intersection(reference_entity_ids)))

            for common_entity_id in common_entity_ids:
                self.assertTrue(common_entity_id in current_entity_ids)

    def __CheckParentEntities(self, model_parts: 'list[Kratos.ModelPart]', examined_model_parts: 'list[Kratos.ModelPart]', entity_container_getter):
        examined_node_ids = []
        for model_part in examined_model_parts:
            examined_node_ids.extend([node.Id for node in model_part.Nodes])

        examined_node_ids = set(examined_node_ids)

        for model_part in model_parts:
            parent_entities = []
            reference_model_part = model_part.GetParentModelPart()
            for entity in entity_container_getter(reference_model_part):
                for node in entity.GetGeometry():
                    if node.Id in examined_node_ids:
                        parent_entities.append(entity)
                        break

            parent_entity_nodes = []
            for parent_entity in parent_entities:
                for node in parent_entity.GetGeometry():
                    parent_entity_nodes.append(node.Id)


            current_node_ids = sorted([node.Id for node in model_part.Nodes])
            current_entity_ids = sorted([entity.Id for entity in entity_container_getter(model_part)])

            self.assertEqual(current_entity_ids, sorted([entity.Id for entity in parent_entities]))
            self.assertEqual(current_node_ids, sorted(list(set(parent_entity_nodes))))

    def __CheckModelParts(self, model_parts: 'list[Kratos.ModelPart]', ref_node_ids, ref_condition_ids, ref_element_ids):
        data_communicator: Kratos.DataCommunicator = self.model_part.GetCommunicator().GetDataCommunicator()

        all_node_ids = []
        all_condition_ids = []
        all_element_ids = []

        for model_part in model_parts:
            self.__CheckModelPart(model_part)

            all_node_ids.extend(TestOptAppModelPartUtils.__GetAllIds(model_part.Nodes, data_communicator))
            all_condition_ids.extend(TestOptAppModelPartUtils.__GetAllIds(model_part.Conditions, data_communicator))
            all_element_ids.extend(TestOptAppModelPartUtils.__GetAllIds(model_part.Elements, data_communicator))

        self.assertEqual(ref_node_ids, sorted(list(set(all_node_ids))))
        self.assertEqual(ref_condition_ids, sorted(list(set(all_condition_ids))))
        self.assertEqual(ref_element_ids, sorted(list(set(all_element_ids))))

    def test_GetSensitivityModelPartForDirectSensitivitiesError(self):
        with self.assertRaises(RuntimeError):
            Kratos.ModelPart = KratosOA.OptAppModelPartUtils.GetModelPartsWithCommonReferenceEntities(
            self.examined_model_parts, self.reference_model_parts, False, False, False, False)

    def test_GetSensitivityModelPartForDirectSensitivitiesNodes(self):
        model_parts: Kratos.ModelPart = KratosOA.OptAppModelPartUtils.GetModelPartsWithCommonReferenceEntities(
            self.examined_model_parts, self.reference_model_parts, True, False, False, False)
        self.__CheckModelParts(model_parts, [3, 4, 11, 12, 13, 16, 17, 18], [], [])
        self.__CheckCommonEntities(model_parts, self.examined_model_parts, lambda x: x.Nodes)

    def test_GetSensitivityModelPartForDirectSensitivitiesNodesConditions(self):
        model_parts: Kratos.ModelPart = KratosOA.OptAppModelPartUtils.GetModelPartsWithCommonReferenceEntities(
            self.examined_model_parts, self.reference_model_parts, True, True, False, False)
        self.__CheckModelParts(model_parts, [3, 4, 11, 12, 13, 16, 17, 18], [11], [])
        self.__CheckCommonEntities(model_parts, self.examined_model_parts, lambda x: x.Nodes)
        self.__CheckCommonEntities(model_parts, self.examined_model_parts, lambda x: x.Conditions)

    def test_GetSensitivityModelPartForDirectSensitivitiesNodesConditionsElements(self):
        model_parts: Kratos.ModelPart = KratosOA.OptAppModelPartUtils.GetModelPartsWithCommonReferenceEntities(
            self.examined_model_parts, self.reference_model_parts, True, True, True, False)
        self.__CheckModelParts(model_parts, [3, 4, 11, 12, 13, 16, 17, 18], [11], [5, 9, 10])
        self.__CheckCommonEntities(model_parts, self.examined_model_parts, lambda x: x.Nodes)
        self.__CheckCommonEntities(model_parts, self.examined_model_parts, lambda x: x.Conditions)
        self.__CheckCommonEntities(model_parts, self.examined_model_parts, lambda x: x.Elements)

    def test_GetSensitivityModelPartForDirectSensitivitiesNodesConditionsWithParents(self):
        model_parts: Kratos.ModelPart = KratosOA.OptAppModelPartUtils.GetModelPartsWithCommonReferenceEntities(
            self.examined_model_parts, [self.model["test"]], True, True, False, True)
        self.__CheckModelParts(model_parts, [1, 2, 3, 4, 5, 9, 10, 11, 12, 13, 15, 16, 17, 18, 20, 21, 22, 23, 24, 25], [1, 2, 3, 4, 6, 9, 10, 11, 12, 13, 14, 15, 16], [])
        self.__CheckCommonEntities(model_parts, self.examined_model_parts, lambda x: x.Nodes)
        self.__CheckCommonEntities(model_parts, self.examined_model_parts, lambda x: x.Conditions)

    def test_GetSensitivityModelPartForDirectSensitivitiesNodesConditionsElementsWithParents(self):
        communicator: Kratos.Communicator = self.model_part.GetCommunicator()
        if (communicator.GetDataCommunicator().IsDistributed() and communicator.TotalProcesses() == 3):
            self.skipTest("This is skipped until the bug in issue 10938 is fixed.")

        model_parts: Kratos.ModelPart = KratosOA.OptAppModelPartUtils.GetModelPartsWithCommonReferenceEntities(
            self.examined_model_parts, [self.model["test"]], True, True, True, True)

        self.__CheckModelParts(model_parts, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25], [1, 2, 3, 4, 6, 9, 10, 11, 12, 13, 14, 15, 16], [1, 2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16])
        self.__CheckCommonEntities(model_parts, self.examined_model_parts, lambda x: x.Nodes)
        self.__CheckCommonEntities(model_parts, self.examined_model_parts, lambda x: x.Conditions)
        self.__CheckCommonEntities(model_parts, self.examined_model_parts, lambda x: x.Elements)
        self.__CheckParentEntities(model_parts, self.examined_model_parts, lambda x: x.Elements)

    def test_GetSensitivityModelPartForDirectSensitivitiesNodesElementsWithParents(self):
        model_parts: Kratos.ModelPart = KratosOA.OptAppModelPartUtils.GetModelPartsWithCommonReferenceEntities(
            self.examined_model_parts, [self.model["test"]], False, False, True, True)
        self.__CheckModelParts(model_parts, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25], [], [1, 2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16])
        self.__CheckCommonEntities(model_parts, self.examined_model_parts, lambda x: x.Nodes)
        self.__CheckCommonEntities(model_parts, self.examined_model_parts, lambda x: x.Elements)
        self.__CheckParentEntities(model_parts, self.examined_model_parts, lambda x: x.Elements)

    def test_ClearEntitiesOfModelPartsWithCommonReferenceEntitiesBetweenReferenceListAndExaminedList(self):
        KratosOA.OptAppModelPartUtils.GetModelPartsWithCommonReferenceEntities(
            self.examined_model_parts, self.reference_model_parts, True, False, False, False)

        self.assertTrue(self.model.HasModelPart("test.evaluated_element_1.<OPTIMIZATION_APP_AUTO>_Nodes_NoConditions_NoElements_NoParents_ExaminedMPs_test>sensitivity_condition_1;test>sensitivity_element_1;test>sensitivity_element_2;test>sensitivity_element_3;test>sensitivity_element_4;"))
        self.assertTrue(self.model.HasModelPart("test.evaluated_element_2.<OPTIMIZATION_APP_AUTO>_Nodes_NoConditions_NoElements_NoParents_ExaminedMPs_test>sensitivity_condition_1;test>sensitivity_element_1;test>sensitivity_element_2;test>sensitivity_element_3;test>sensitivity_element_4;"))
        self.assertTrue(self.model.HasModelPart("test.evaluated_element_3.<OPTIMIZATION_APP_AUTO>_Nodes_NoConditions_NoElements_NoParents_ExaminedMPs_test>sensitivity_condition_1;test>sensitivity_element_1;test>sensitivity_element_2;test>sensitivity_element_3;test>sensitivity_element_4;"))

        KratosOA.OptAppModelPartUtils.RemoveModelPartsWithCommonReferenceEntitiesBetweenReferenceListAndExaminedList([self.model["test"]])

        self.assertFalse(self.model.HasModelPart("test.evaluated_element_1.<OPTIMIZATION_APP_AUTO>_Nodes_NoConditions_NoElements_NoParents_ExaminedMPs_test>sensitivity_condition_1;test>sensitivity_element_1;test>sensitivity_element_2;test>sensitivity_element_3;test>sensitivity_element_4;"))
        self.assertFalse(self.model.HasModelPart("test.evaluated_element_2.<OPTIMIZATION_APP_AUTO>_Nodes_NoConditions_NoElements_NoParents_ExaminedMPs_test>sensitivity_condition_1;test>sensitivity_element_1;test>sensitivity_element_2;test>sensitivity_element_3;test>sensitivity_element_4;"))
        self.assertFalse(self.model.HasModelPart("test.evaluated_element_3.<OPTIMIZATION_APP_AUTO>_Nodes_NoConditions_NoElements_NoParents_ExaminedMPs_test>sensitivity_condition_1;test>sensitivity_element_1;test>sensitivity_element_2;test>sensitivity_element_3;test>sensitivity_element_4;"))

        KratosOA.OptAppModelPartUtils.GetModelPartsWithCommonReferenceEntities(
            self.examined_model_parts, self.reference_model_parts, True, False, False, False)

        self.assertTrue(self.model.HasModelPart("test.evaluated_element_1.<OPTIMIZATION_APP_AUTO>_Nodes_NoConditions_NoElements_NoParents_ExaminedMPs_test>sensitivity_condition_1;test>sensitivity_element_1;test>sensitivity_element_2;test>sensitivity_element_3;test>sensitivity_element_4;"))
        self.assertTrue(self.model.HasModelPart("test.evaluated_element_2.<OPTIMIZATION_APP_AUTO>_Nodes_NoConditions_NoElements_NoParents_ExaminedMPs_test>sensitivity_condition_1;test>sensitivity_element_1;test>sensitivity_element_2;test>sensitivity_element_3;test>sensitivity_element_4;"))
        self.assertTrue(self.model.HasModelPart("test.evaluated_element_3.<OPTIMIZATION_APP_AUTO>_Nodes_NoConditions_NoElements_NoParents_ExaminedMPs_test>sensitivity_condition_1;test>sensitivity_element_1;test>sensitivity_element_2;test>sensitivity_element_3;test>sensitivity_element_4;"))

if __name__ == "__main__":
    kratos_unittest.main()