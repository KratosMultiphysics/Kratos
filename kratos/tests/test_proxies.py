# --- Core Imports ---
import KratosMultiphysics
from KratosMultiphysics.container_proxy import ContainerProxy
import KratosMultiphysics.KratosUnittest


class TestEntityProxy(KratosMultiphysics.KratosUnittest.TestCase):

    def setUp(self) -> None:
        self.model = KratosMultiphysics.Model()
        self.model_part = self.model.CreateModelPart("root")
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)
        properties = KratosMultiphysics.Properties(1)

        self.model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        self.model_part.CreateNewNode(2, 0.0, 0.0, 0.0)
        self.model_part.CreateNewNode(3, 0.0, 0.0, 0.0)
        self.model_part.CreateNewNode(4, 0.0, 0.0, 0.0)
        for node in self.model_part.Nodes:
            node[KratosMultiphysics.PRESSURE] = node.Id
            node.SetSolutionStepValue(KratosMultiphysics.PRESSURE, node.Id)

        self.model_part.CreateNewElement("Element2D3N", 1, [1, 2, 3], properties)
        self.model_part.CreateNewElement("Element2D3N", 2, [3, 4, 1], properties)
        for element in self.model_part.Elements:
            element[KratosMultiphysics.PRESSURE] = element.Id

        self.model_part.CreateNewCondition("LineCondition2D2N", 1, [1, 2], properties)
        self.model_part.CreateNewCondition("LineCondition2D2N", 2, [3, 4], properties)
        for condition in self.model_part.Conditions:
            condition[KratosMultiphysics.PRESSURE] = condition.Id

        self.model_part.ProcessInfo[KratosMultiphysics.PRESSURE] = 2.0
        self.model_part[KratosMultiphysics.PRESSURE] = 2.0

    def test_HasValue(self) -> None:
        proxy = KratosMultiphysics.EntityProxy(KratosMultiphysics.Globals.DataLocation.ProcessInfo, self.model_part)
        self.assertTrue(proxy.HasValue(KratosMultiphysics.PRESSURE))
        self.assertFalse(proxy.HasValue(KratosMultiphysics.VELOCITY))

        proxy = KratosMultiphysics.EntityProxy(KratosMultiphysics.Globals.DataLocation.ProcessInfo, self.model_part.ProcessInfo)
        self.assertTrue(proxy.HasValue(KratosMultiphysics.PRESSURE))
        self.assertFalse(proxy.HasValue(KratosMultiphysics.VELOCITY))

        for entity_type in (KratosMultiphysics.Globals.DataLocation.NodeHistorical,
                            KratosMultiphysics.Globals.DataLocation.NodeNonHistorical,
                            KratosMultiphysics.Globals.DataLocation.Element,
                            KratosMultiphysics.Globals.DataLocation.Condition):
            for proxy in ContainerProxy(entity_type, self.model_part):
                self.assertTrue(proxy.HasValue(KratosMultiphysics.PRESSURE))
                self.assertFalse(proxy.HasValue(KratosMultiphysics.VELOCITY))

    def test_GetValue(self) -> None:
        proxy = KratosMultiphysics.EntityProxy(KratosMultiphysics.Globals.DataLocation.ProcessInfo, self.model_part)
        self.assertEqual(proxy.GetValue(KratosMultiphysics.PRESSURE), 2.0)

        proxy = KratosMultiphysics.EntityProxy(KratosMultiphysics.Globals.DataLocation.ProcessInfo, self.model_part.ProcessInfo)
        self.assertEqual(proxy.GetValue(KratosMultiphysics.PRESSURE), 2.0)

        for entity_type in (KratosMultiphysics.Globals.DataLocation.NodeHistorical,
                            KratosMultiphysics.Globals.DataLocation.NodeNonHistorical,
                            KratosMultiphysics.Globals.DataLocation.Element,
                            KratosMultiphysics.Globals.DataLocation.Condition):
            for i_entity, proxy in enumerate(ContainerProxy(entity_type, self.model_part)):
                id_entity = i_entity + 1
                self.assertEqual(proxy.GetValue(KratosMultiphysics.PRESSURE), id_entity)

    def test_SetValue(self) -> None:
        process_info_proxy = KratosMultiphysics.EntityProxy(KratosMultiphysics.Globals.DataLocation.ProcessInfo, self.model_part)
        process_info_proxy.SetValue(KratosMultiphysics.PRESSURE, process_info_proxy.GetValue(KratosMultiphysics.PRESSURE) * 2)

        model_part_proxy = KratosMultiphysics.EntityProxy(KratosMultiphysics.Globals.DataLocation.ProcessInfo, self.model_part)
        model_part_proxy.SetValue(KratosMultiphysics.PRESSURE, model_part_proxy.GetValue(KratosMultiphysics.PRESSURE) * 2)

        self.assertAlmostEqual(process_info_proxy.GetValue(KratosMultiphysics.PRESSURE), 4.0)
        self.assertAlmostEqual(model_part_proxy.GetValue(KratosMultiphysics.PRESSURE), 4.0)

        process_info_proxy.SetValue(KratosMultiphysics.PRESSURE, 2.0)
        model_part_proxy.SetValue(KratosMultiphysics.PRESSURE, 2.0)

        self.assertAlmostEqual(process_info_proxy.GetValue(KratosMultiphysics.PRESSURE), 2.0)
        self.assertAlmostEqual(model_part_proxy.GetValue(KratosMultiphysics.PRESSURE), 2.0)

        for entity_type in (KratosMultiphysics.Globals.DataLocation.NodeHistorical,
                            KratosMultiphysics.Globals.DataLocation.NodeNonHistorical,
                            KratosMultiphysics.Globals.DataLocation.Element,
                            KratosMultiphysics.Globals.DataLocation.Condition):
            for i_entity, proxy in enumerate(ContainerProxy(entity_type, self.model_part)):
                id_entity = i_entity + 1
                self.assertEqual(proxy.GetValue(KratosMultiphysics.PRESSURE), id_entity)
                proxy.SetValue(KratosMultiphysics.PRESSURE, 2 * id_entity)
                self.assertEqual(proxy.GetValue(KratosMultiphysics.PRESSURE), 2 * id_entity)
                proxy.SetValue(KratosMultiphysics.PRESSURE, id_entity)
                self.assertEqual(proxy.GetValue(KratosMultiphysics.PRESSURE), id_entity)


if __name__ == "__main__":
    KratosMultiphysics.KratosUnittest.main()
