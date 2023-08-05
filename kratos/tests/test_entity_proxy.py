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
            id = node.Id
            node[KratosMultiphysics.PRESSURE] = id
            node.SetSolutionStepValue(KratosMultiphysics.PRESSURE, id)

        self.model_part.CreateNewElement("Element2D3N", 1, [1, 2, 3], properties)
        self.model_part.CreateNewElement("Element2D3N", 2, [3, 4, 1], properties)
        for element in self.model_part.Elements:
            id = element.Id
            element[KratosMultiphysics.PRESSURE] = id

        self.model_part.CreateNewCondition("LineCondition2D2N", 1, [1, 2], properties)
        self.model_part.CreateNewCondition("LineCondition2D2N", 2, [3, 4], properties)
        for condition in self.model_part.Conditions:
            id = condition.Id
            condition[KratosMultiphysics.PRESSURE] = id

    def test_HasValue(self) -> None:
        for entity_type in (KratosMultiphysics.Globals.DataLocation.NodeHistorical,
                            KratosMultiphysics.Globals.DataLocation.NodeNonHistorical,
                            KratosMultiphysics.Globals.DataLocation.Element,
                            KratosMultiphysics.Globals.DataLocation.Condition):
            for proxy in ContainerProxy(self.model_part, entity_type):
                self.assertTrue(proxy.HasValue(KratosMultiphysics.PRESSURE))
                self.assertFalse(proxy.HasValue(KratosMultiphysics.VELOCITY))

    def test_GetValue(self) -> None:
        for entity_type in (KratosMultiphysics.Globals.DataLocation.NodeHistorical,
                            KratosMultiphysics.Globals.DataLocation.NodeNonHistorical,
                            KratosMultiphysics.Globals.DataLocation.Element,
                            KratosMultiphysics.Globals.DataLocation.Condition):
            for i_entity, proxy in enumerate(ContainerProxy(self.model_part, entity_type)):
                id = i_entity + 1
                self.assertEqual(proxy.GetValue(KratosMultiphysics.PRESSURE), id)

    def test_SetValue(self) -> None:
        for entity_type in (KratosMultiphysics.Globals.DataLocation.NodeHistorical,
                            KratosMultiphysics.Globals.DataLocation.NodeNonHistorical,
                            KratosMultiphysics.Globals.DataLocation.Element,
                            KratosMultiphysics.Globals.DataLocation.Condition):
            for i_entity, proxy in enumerate(ContainerProxy(self.model_part, entity_type)):
                id = i_entity + 1
                self.assertEqual(proxy.GetValue(KratosMultiphysics.PRESSURE), id)
                proxy.SetValue(KratosMultiphysics.PRESSURE, 2 * id)
                self.assertEqual(proxy.GetValue(KratosMultiphysics.PRESSURE), 2 * id)
                proxy.SetValue(KratosMultiphysics.PRESSURE, id)
                self.assertEqual(proxy.GetValue(KratosMultiphysics.PRESSURE), id)


if __name__ == "__main__":
    KratosMultiphysics.KratosUnittest.main()
