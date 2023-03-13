import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as kratos_unittest
from KratosMultiphysics.testing.utilities import ReadModelPart

class TestModelPartUtils(kratos_unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")
        cls.model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 2
        ReadModelPart("mdpas/quads", cls.model_part)

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

    def __CheckNodeIds(self, model_part:Kratos.ModelPart, ref_node_ids: 'list[int]'):
        mp_node_ids = []
        for node in model_part.Nodes:
            mp_node_ids.append(node.Id)
        mp_node_ids = sorted(mp_node_ids)
        self.assertEqual(ref_node_ids, mp_node_ids)

    def __CheckConditionIds(self, model_part:Kratos.ModelPart, ref_condition_ids: 'list[int]'):
        mp_condition_ids = []
        for condition in model_part.Conditions:
            mp_condition_ids.append(condition.Id)
        mp_condition_ids = sorted(mp_condition_ids)
        self.assertEqual(ref_condition_ids, mp_condition_ids)

    def __CheckElementIds(self, model_part:Kratos.ModelPart, ref_element_ids: 'list[int]'):
        mp_element_ids = []
        for element in model_part.Elements:
            mp_element_ids.append(element.Id)
        mp_element_ids = sorted(mp_element_ids)
        self.assertEqual(ref_element_ids, mp_element_ids)

    @staticmethod
    def __GetValidEntityIdsList(entities, ref_ids):
        return sorted([entity.Id for entity in entities if entity.Id in ref_ids])

    def __CheckModelParts(self, model_parts: 'list[Kratos.ModelPart]', ref_node_ids, ref_condition_ids, ref_element_ids):
        # create local refs in mpi tests
        local_ref_node_ids = TestModelPartUtils.__GetValidEntityIdsList(self.model_part.Nodes, ref_node_ids)
        local_ref_condition_ids = TestModelPartUtils.__GetValidEntityIdsList(self.model_part.Conditions, ref_condition_ids)
        local_ref_element_ids = TestModelPartUtils.__GetValidEntityIdsList(self.model_part.Elements, ref_element_ids)

        for model_part in model_parts:
            self.__CheckModelPart(model_part)
            self.__CheckNodeIds(model_part, TestModelPartUtils.__GetValidEntityIdsList(model_part.GetParentModelPart().Nodes, local_ref_node_ids))
            self.__CheckConditionIds(model_part, TestModelPartUtils.__GetValidEntityIdsList(model_part.GetParentModelPart().Conditions, local_ref_condition_ids))
            self.__CheckElementIds(model_part, TestModelPartUtils.__GetValidEntityIdsList(model_part.GetParentModelPart().Elements, local_ref_element_ids))

    def test_GetSensitivityModelPartForDirectSensitivitiesError(self):
        with self.assertRaises(RuntimeError):
            Kratos.ModelPart = KratosOA.ModelPartUtils.GetModelPartsWithCommonReferenceEntitiesBetweenReferenceListAndExaminedList(
            self.examined_model_parts, self.reference_model_parts, False, False, False, False)

    def test_GetSensitivityModelPartForDirectSensitivitiesNodes(self):
        model_parts: Kratos.ModelPart = KratosOA.ModelPartUtils.GetModelPartsWithCommonReferenceEntitiesBetweenReferenceListAndExaminedList(
            self.examined_model_parts, self.reference_model_parts, True, False, False, False)
        self.__CheckModelParts(model_parts, [3, 4, 11, 12, 13, 16, 17, 18], [], [])

    def test_GetSensitivityModelPartForDirectSensitivitiesNodesConditions(self):
        model_parts: Kratos.ModelPart = KratosOA.ModelPartUtils.GetModelPartsWithCommonReferenceEntitiesBetweenReferenceListAndExaminedList(
            self.examined_model_parts, self.reference_model_parts, True, True, False, False)
        self.__CheckModelParts(model_parts, [3, 4, 11, 12, 13, 16, 17, 18], [11], [])

    def test_GetSensitivityModelPartForDirectSensitivitiesNodesConditionsElements(self):
        model_parts: Kratos.ModelPart = KratosOA.ModelPartUtils.GetModelPartsWithCommonReferenceEntitiesBetweenReferenceListAndExaminedList(
            self.examined_model_parts, self.reference_model_parts, True, True, True, False)
        self.__CheckModelParts(model_parts, [3, 4, 11, 12, 13, 16, 17, 18], [11], [5, 9, 10])

    def test_GetSensitivityModelPartForDirectSensitivitiesNodesConditionsWithParents(self):
        model_parts: Kratos.ModelPart = KratosOA.ModelPartUtils.GetModelPartsWithCommonReferenceEntitiesBetweenReferenceListAndExaminedList(
            self.examined_model_parts, [self.model["test"]], True, True, False, True)
        self.__CheckModelParts(model_parts, [1, 2, 3, 4, 5, 9, 10, 11, 12, 13, 15, 16, 17, 18, 20, 21, 22, 23, 24, 25], [1, 2, 3, 4, 6, 9, 10, 11, 12, 13, 14, 15, 16], [])

    def test_GetSensitivityModelPartForDirectSensitivitiesNodesConditionsElementsWithParents(self):
        model_parts: Kratos.ModelPart = KratosOA.ModelPartUtils.GetModelPartsWithCommonReferenceEntitiesBetweenReferenceListAndExaminedList(
            self.examined_model_parts, [self.model["test"]], True, True, True, True)
        self.__CheckModelParts(model_parts, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25], [1, 2, 3, 4, 6, 9, 10, 11, 12, 13, 14, 15, 16], [1, 2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16])

    def test_GetSensitivityModelPartForDirectSensitivitiesNodesElementsWithParents(self):
        model_parts: Kratos.ModelPart = KratosOA.ModelPartUtils.GetModelPartsWithCommonReferenceEntitiesBetweenReferenceListAndExaminedList(
            self.examined_model_parts, [self.model["test"]], False, False, True, True)
        self.__CheckModelParts(model_parts, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25], [], [1, 2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16])

    def test_ClearEntitiesOfModelPartsWithCommonReferenceEntitiesBetweenReferenceListAndExaminedList(self):
        KratosOA.ModelPartUtils.GetModelPartsWithCommonReferenceEntitiesBetweenReferenceListAndExaminedList(
            self.examined_model_parts, self.reference_model_parts, True, False, False, False)

        self.assertTrue(self.model.HasModelPart("test.evaluated_element_1.<OPTIMIZATION_APP_AUTO>_Nodes_NoConditions_NoElements_NoParents_ExaminedMPs_test>sensitivity_condition_1;test>sensitivity_element_1;test>sensitivity_element_2;test>sensitivity_element_3;test>sensitivity_element_4;"))
        self.assertTrue(self.model.HasModelPart("test.evaluated_element_2.<OPTIMIZATION_APP_AUTO>_Nodes_NoConditions_NoElements_NoParents_ExaminedMPs_test>sensitivity_condition_1;test>sensitivity_element_1;test>sensitivity_element_2;test>sensitivity_element_3;test>sensitivity_element_4;"))
        self.assertTrue(self.model.HasModelPart("test.evaluated_element_3.<OPTIMIZATION_APP_AUTO>_Nodes_NoConditions_NoElements_NoParents_ExaminedMPs_test>sensitivity_condition_1;test>sensitivity_element_1;test>sensitivity_element_2;test>sensitivity_element_3;test>sensitivity_element_4;"))

        KratosOA.ModelPartUtils.RemoveModelPartsWithCommonReferenceEntitiesBetweenReferenceListAndExaminedList([self.model["test"]])

        self.assertFalse(self.model.HasModelPart("test.evaluated_element_1.<OPTIMIZATION_APP_AUTO>_Nodes_NoConditions_NoElements_NoParents_ExaminedMPs_test>sensitivity_condition_1;test>sensitivity_element_1;test>sensitivity_element_2;test>sensitivity_element_3;test>sensitivity_element_4;"))
        self.assertFalse(self.model.HasModelPart("test.evaluated_element_2.<OPTIMIZATION_APP_AUTO>_Nodes_NoConditions_NoElements_NoParents_ExaminedMPs_test>sensitivity_condition_1;test>sensitivity_element_1;test>sensitivity_element_2;test>sensitivity_element_3;test>sensitivity_element_4;"))
        self.assertFalse(self.model.HasModelPart("test.evaluated_element_3.<OPTIMIZATION_APP_AUTO>_Nodes_NoConditions_NoElements_NoParents_ExaminedMPs_test>sensitivity_condition_1;test>sensitivity_element_1;test>sensitivity_element_2;test>sensitivity_element_3;test>sensitivity_element_4;"))

        KratosOA.ModelPartUtils.GetModelPartsWithCommonReferenceEntitiesBetweenReferenceListAndExaminedList(
            self.examined_model_parts, self.reference_model_parts, True, False, False, False)

        self.assertTrue(self.model.HasModelPart("test.evaluated_element_1.<OPTIMIZATION_APP_AUTO>_Nodes_NoConditions_NoElements_NoParents_ExaminedMPs_test>sensitivity_condition_1;test>sensitivity_element_1;test>sensitivity_element_2;test>sensitivity_element_3;test>sensitivity_element_4;"))
        self.assertTrue(self.model.HasModelPart("test.evaluated_element_2.<OPTIMIZATION_APP_AUTO>_Nodes_NoConditions_NoElements_NoParents_ExaminedMPs_test>sensitivity_condition_1;test>sensitivity_element_1;test>sensitivity_element_2;test>sensitivity_element_3;test>sensitivity_element_4;"))
        self.assertTrue(self.model.HasModelPart("test.evaluated_element_3.<OPTIMIZATION_APP_AUTO>_Nodes_NoConditions_NoElements_NoParents_ExaminedMPs_test>sensitivity_condition_1;test>sensitivity_element_1;test>sensitivity_element_2;test>sensitivity_element_3;test>sensitivity_element_4;"))

if __name__ == "__main__":
    Kratos.Tester.SetVerbosity(Kratos.Tester.Verbosity.PROGRESS)  # TESTS_OUTPUTS
    kratos_unittest.main()