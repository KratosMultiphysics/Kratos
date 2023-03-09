import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as kratos_unittest

class TestResponseUtils(kratos_unittest.TestCase):
    @staticmethod
    def __get_l2_norm_square(value: Kratos.Array3) -> float:
        return value[0]**2 + value[1]**2 + value[2]**2

    @staticmethod
    def __get_node_id(model_part: Kratos.ModelPart, node_pos: Kratos.Array3, id_offset = 0) -> int:
        node: Kratos.Node
        for node in model_part.Nodes:
            current_node_pos = Kratos.Array3([node.X, node.Y, node.Z])
            if TestResponseUtils.__get_l2_norm_square(current_node_pos - node_pos) < 1e-12:
                return node.Id

        node = model_part.CreateNewNode(model_part.NumberOfNodes() + 1 + id_offset, node_pos[0], node_pos[1], node_pos[2])
        return node.Id

    @staticmethod
    def __generate_quad_elements(model_part: Kratos.ModelPart, rows: int, cols: int, id_offset = 0):
        properties = model_part.CreateNewProperties(1)
        for i in range(rows):
            for j in range(cols):
                center = Kratos.Array3([i, j, 0.0])
                node_id_1 = TestResponseUtils.__get_node_id(model_part, center + Kratos.Array3([-0.5, -0.5, 0.0]), id_offset)
                node_id_2 = TestResponseUtils.__get_node_id(model_part, center + Kratos.Array3([-0.5, +0.5, 0.0]), id_offset)
                node_id_3 = TestResponseUtils.__get_node_id(model_part, center + Kratos.Array3([+0.5, +0.5, 0.0]), id_offset)
                node_id_4 = TestResponseUtils.__get_node_id(model_part, center + Kratos.Array3([+0.5, -0.5, 0.0]), id_offset)

                model_part.CreateNewElement("Element2D4N", model_part.NumberOfElements() + 1 + id_offset, [node_id_1, node_id_2, node_id_3, node_id_4], properties)

    @staticmethod
    def __generate_line_conditions(model_part: Kratos.ModelPart, rows: int, cols: int, id_offset = 0):
        properties = model_part.CreateNewProperties(2)
        for i in range(rows):
            # left conditions
            center = Kratos.Array3([-0.5, i, 0.0])
            node_id_1 = TestResponseUtils.__get_node_id(model_part, center + Kratos.Array3([0.0, -0.5, 0.0]), id_offset)
            node_id_2 = TestResponseUtils.__get_node_id(model_part, center + Kratos.Array3([0.0, +0.5, 0.0]), id_offset)
            model_part.CreateNewCondition("LineCondition2D2N", model_part.NumberOfConditions() + 1 + id_offset, [node_id_1, node_id_2], properties)

            # right conditions
            center = Kratos.Array3([cols - 0.5, i, 0.0])
            node_id_1 = TestResponseUtils.__get_node_id(model_part, center + Kratos.Array3([0.0, -0.5, 0.0]), id_offset)
            node_id_2 = TestResponseUtils.__get_node_id(model_part, center + Kratos.Array3([0.0, +0.5, 0.0]), id_offset)
            model_part.CreateNewCondition("LineCondition2D2N", model_part.NumberOfConditions() + 1 + id_offset, [node_id_1, node_id_2], properties)

        for i in range(cols):
            # bottom conditions
            center = Kratos.Array3([i, -0.5, 0.0])
            node_id_1 = TestResponseUtils.__get_node_id(model_part, center + Kratos.Array3([-0.5, 0.0, 0.0]), id_offset)
            node_id_2 = TestResponseUtils.__get_node_id(model_part, center + Kratos.Array3([+0.5, 0.0, 0.0]), id_offset)
            model_part.CreateNewCondition("LineCondition2D2N", model_part.NumberOfConditions() + 1 + id_offset, [node_id_1, node_id_2], properties)

            # top conditions
            center = Kratos.Array3([i, rows - 0.5, 0.0])
            node_id_1 = TestResponseUtils.__get_node_id(model_part, center + Kratos.Array3([-0.5, 0.0, 0.0]), id_offset)
            node_id_2 = TestResponseUtils.__get_node_id(model_part, center + Kratos.Array3([+0.5, 0.0, 0.0]), id_offset)
            model_part.CreateNewCondition("LineCondition2D2N", model_part.NumberOfConditions() + 1 + id_offset, [node_id_1, node_id_2], properties)

    @staticmethod
    def __add_elements(model_part: Kratos.ModelPart, element_ids: 'list[int]'):
        model_part.AddElements(element_ids)
        node_ids = []
        for element_id in element_ids:
            for node in model_part.GetElement(element_id).GetGeometry():
                node_ids.append(node.Id)

        node_ids = list(set(node_ids))
        model_part.AddNodes(node_ids)

    @staticmethod
    def __add_conditions(model_part: Kratos.ModelPart, condition_ids: 'list[int]'):
        model_part.AddConditions(condition_ids)
        node_ids = []
        for condition_id in condition_ids:
            for node in model_part.GetCondition(condition_id).GetGeometry():
                node_ids.append(node.Id)

        node_ids = list(set(node_ids))
        model_part.AddNodes(node_ids)

    @classmethod
    def setUpClass(cls) -> None:
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")

        TestResponseUtils.__generate_quad_elements(cls.model_part, 4, 4)
        TestResponseUtils.__generate_line_conditions(cls.model_part, 4, 4)

        # now create the sub model parts
        model_part = cls.model_part.CreateSubModelPart("evaluated_element_1")
        TestResponseUtils.__add_elements(model_part, [5, 9, 7])
        TestResponseUtils.__add_conditions(model_part, [11, 13])
        model_part = cls.model_part.CreateSubModelPart("evaluated_element_2")
        TestResponseUtils.__add_elements(model_part, [6, 10, 11])
        model_part = cls.model_part.CreateSubModelPart("evaluated_element_3")
        TestResponseUtils.__add_elements(model_part, [6, 10])

        model_part = cls.model_part.CreateSubModelPart("sensitivity_element_1")
        TestResponseUtils.__add_elements(model_part, [1, 5])
        TestResponseUtils.__add_conditions(model_part, [11])
        model_part = cls.model_part.CreateSubModelPart("sensitivity_element_2")
        TestResponseUtils.__add_elements(model_part, [1, 5, 9])
        model_part = cls.model_part.CreateSubModelPart("sensitivity_element_3")
        TestResponseUtils.__add_elements(model_part, [9, 13])
        model_part = cls.model_part.CreateSubModelPart("sensitivity_element_4")
        TestResponseUtils.__add_elements(model_part, [10, 14])

        model_part = cls.model_part.CreateSubModelPart("sensitivity_condition_1")
        TestResponseUtils.__add_conditions(model_part, [10, 12])
        model_part = cls.model_part.CreateSubModelPart("sensitivity_condition_2")
        TestResponseUtils.__add_conditions(model_part, [10, 12, 14])
        model_part = cls.model_part.CreateSubModelPart("sensitivity_condition_3")
        TestResponseUtils.__add_conditions(model_part, [14, 16])
        model_part = cls.model_part.CreateSubModelPart("sensitivity_condition_4")
        TestResponseUtils.__add_conditions(model_part, [8, 6])

        cls.sensitivity_model_parts_list = [
            cls.model["test.sensitivity_element_1"],
            cls.model["test.sensitivity_element_2"],
            cls.model["test.sensitivity_element_3"],
            cls.model["test.sensitivity_element_4"],
            cls.model["test.sensitivity_condition_1"],
            cls.model["test.sensitivity_condition_2"],
            cls.model["test.sensitivity_condition_3"],
            cls.model["test.sensitivity_condition_4"]
        ]
        cls.evaluated_model_parts_list = [
            cls.model["test.evaluated_element_1"],
            cls.model["test.evaluated_element_2"],
            cls.model["test.evaluated_element_3"]
        ]

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

    def test_GetSensitivityModelPartForDirectSensitivitiesNone(self):
        with self.assertRaises(RuntimeError):
            KratosOA.ResponseUtils.GetSensitivityModelPartForDirectSensitivities(self.sensitivity_model_parts_list, self.evaluated_model_parts_list, False, False, False)

    def test_GetSensitivityModelPartForDirectSensitivitiesNodes(self):
        model_part: Kratos.ModelPart = KratosOA.ResponseUtils.GetSensitivityModelPartForDirectSensitivities(self.sensitivity_model_parts_list, self.evaluated_model_parts_list, True, False, False)
        self.__CheckNodeIds(model_part, [3, 4, 11, 12, 13, 16, 17, 18])
        self.__CheckConditionIds(model_part, [])
        self.__CheckElementIds(model_part, [])

    def test_GetSensitivityModelPartForDirectSensitivitiesConditions(self):
        model_part: Kratos.ModelPart = KratosOA.ResponseUtils.GetSensitivityModelPartForDirectSensitivities(self.sensitivity_model_parts_list, self.evaluated_model_parts_list, False, True, False)
        self.assertEqual(model_part.NumberOfElements(), 0)

        self.__CheckNodeIds(model_part, [4, 12])
        self.__CheckConditionIds(model_part, [11])
        self.__CheckElementIds(model_part, [])

        model_part: Kratos.ModelPart = KratosOA.ResponseUtils.GetSensitivityModelPartForDirectSensitivities(self.sensitivity_model_parts_list, self.evaluated_model_parts_list, True, True, False)
        self.assertEqual(model_part.NumberOfElements(), 0)

        self.__CheckNodeIds(model_part, [3, 4, 11, 12, 13, 16, 17, 18])
        self.__CheckConditionIds(model_part, [11])
        self.__CheckElementIds(model_part, [])

    def test_GetSensitivityModelPartForDirectSensitivitiesElements(self):
        model_part_1: Kratos.ModelPart = KratosOA.ResponseUtils.GetSensitivityModelPartForDirectSensitivities(self.sensitivity_model_parts_list, self.evaluated_model_parts_list, False, False, True)
        self.__CheckNodeIds(model_part_1, [3, 4, 11, 12, 13, 16, 17, 18])
        self.__CheckConditionIds(model_part_1, [])
        self.__CheckElementIds(model_part_1, [5, 9, 10])

        model_part_2: Kratos.ModelPart = KratosOA.ResponseUtils.GetSensitivityModelPartForDirectSensitivities(self.sensitivity_model_parts_list, self.evaluated_model_parts_list, False, True, True)
        self.__CheckNodeIds(model_part_2, [3, 4, 11, 12, 13, 16, 17, 18])
        self.__CheckConditionIds(model_part_2, [11])
        self.__CheckElementIds(model_part_2, [5, 9, 10])

        model_part_3: Kratos.ModelPart = KratosOA.ResponseUtils.GetSensitivityModelPartForDirectSensitivities(self.sensitivity_model_parts_list, self.evaluated_model_parts_list, True, True, True)
        self.__CheckNodeIds(model_part_3, [3, 4, 11, 12, 13, 16, 17, 18])
        self.__CheckConditionIds(model_part_3, [11])
        self.__CheckElementIds(model_part_3, [5, 9, 10])

        self.assertNotEqual(model_part_1.FullName(), model_part_2.FullName())
        self.assertNotEqual(model_part_1.FullName(), model_part_3.FullName())
        self.assertNotEqual(model_part_2.FullName(), model_part_3.FullName())

    def test_GetSensitivityModelPartForAdjointSensitivitiesNone(self):
        with self.assertRaises(RuntimeError):
            KratosOA.ResponseUtils.GetSensitivityModelPartForAdjointSensitivities(self.sensitivity_model_parts_list, self.model["test"], False, False, False)

    def test_GetSensitivityModelPartForAdjointSensitivitiesSensitivityEntities(self):
        model_part = KratosOA.ResponseUtils.GetSensitivityModelPartForAdjointSensitivities(self.sensitivity_model_parts_list, self.model["test"], False, True, False)
        self.__CheckNodeIds(model_part, [1, 2, 3, 4, 9, 10, 11, 12, 13, 15, 16, 17, 18, 20, 21, 22, 23, 24, 25])
        self.__CheckConditionIds(model_part, [6, 8, 10, 11, 12, 14, 16])
        self.__CheckElementIds(model_part, [1, 5, 9, 10, 13, 14])

    def test_GetSensitivityModelPartForAdjointSensitivitiesParentEntities(self):
        model_part = KratosOA.ResponseUtils.GetSensitivityModelPartForAdjointSensitivities(self.sensitivity_model_parts_list, self.model["test"], True, False, False)
        self.__CheckElementIds(model_part, [1, 2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16])
        self.__CheckConditionIds(model_part, [1, 2, 3, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16])
        self.__CheckNodeIds(model_part, list(range(1, 26, 1)))

    def test_GetSensitivityModelPartForAdjointSensitivitiesForcedEntities1(self):
        duplicate_model_part = self.model.CreateModelPart("duplicated_model_part1")
        for node in self.model["test"].Nodes:
            duplicate_model_part.AddNode(node)
        TestResponseUtils.__generate_quad_elements(duplicate_model_part, 4, 4, 100)
        TestResponseUtils.__generate_line_conditions(duplicate_model_part, 4, 4, 100)

        model_part = KratosOA.ResponseUtils.GetSensitivityModelPartForAdjointSensitivities(self.sensitivity_model_parts_list, duplicate_model_part, False, True, False)
        self.__CheckElementIds(model_part, [101, 105, 109, 110, 113, 114])
        self.__CheckConditionIds(model_part, [106, 108, 110, 111, 112, 114, 116])
        self.__CheckNodeIds(model_part, [1, 2, 3, 4, 9, 10, 11, 12, 13, 15, 16, 17, 18, 20, 21, 22, 23, 24, 25])

    def test_GetSensitivityModelPartForAdjointSensitivitiesForcedEntities2(self):
        duplicate_model_part = self.model.CreateModelPart("duplicated_model_part2")
        for node in self.model["test"].Nodes:
            duplicate_model_part.AddNode(node)
        TestResponseUtils.__generate_quad_elements(duplicate_model_part, 4, 4, 100)
        TestResponseUtils.__generate_line_conditions(duplicate_model_part, 4, 4, 100)

        model_part = KratosOA.ResponseUtils.GetSensitivityModelPartForAdjointSensitivities(self.sensitivity_model_parts_list, duplicate_model_part, True, False, False)
        self.__CheckElementIds(model_part, [101, 102, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116])
        self.__CheckConditionIds(model_part, [101, 102, 103, 104, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116])
        self.__CheckNodeIds(model_part, list(range(1, 26, 1)))

if __name__ == "__main__":
    Kratos.Tester.SetVerbosity(Kratos.Tester.Verbosity.PROGRESS)  # TESTS_OUTPUTS
    kratos_unittest.main()