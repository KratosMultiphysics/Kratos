import json
import math
from collections import Counter

import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.multipatch_utilities import get_gauss_points
from KratosMultiphysics.model_parameters_factory import (
    KratosModelParametersFactory,
)
from KratosMultiphysics.modelers.multipatch_modeler import (
    FemMultipatchModeler,
)


class TestFemMultipatchModeler(KratosUnittest.TestCase):

    @staticmethod
    def _settings(patches, background_type="triangle"):
        return KM.Parameters(json.dumps({
            "background": {"origin": [-1.0, -0.5],
                           "dimensions": [6.0, 4.0], "divisions": [24, 16],
                           "element_type": background_type},
            "internal_patches": patches,
        }))

    @staticmethod
    def _patches():
        return [
            {"origin": [0.0, 0.5], "dimensions": [1.5, 2.0],
             "divisions": [7, 9]},
            {"origin": [2.5, 0.25], "dimensions": [1.0, 1.5],
             "divisions": [5, 8]},
            {"origin": [3.75, 2.0], "dimensions": [0.5, 0.75],
             "divisions": [3, 4]},
        ]

    def _build(self, patches, background_type="triangle"):
        model = KM.Model()
        modeler = FemMultipatchModeler(
            model, self._settings(patches, background_type))
        modeler.SetupGeometryModel()
        modeler.PrepareGeometryModel()
        modeler.SetupModelPart()
        return model, modeler

    def _check_partition(self, patches, background_type="triangle"):
        model, modeler = self._build(patches, background_type)
        background = modeler.patches[0]
        background_factor = 2 if background_type == "triangle" else 1
        self.assertEqual(background.NumberOfElements(),
                         384 * background_factor)
        self.assertIsInstance(modeler.background_utility,
                              KM.FemMultipatchUtilities)
        self.assertEqual(background.NumberOfNodes(), 425)
        self.assertEqual(len(modeler.patches), len(patches) + 1)
        node_ids = set()
        element_ids = set()
        active_area = 0.0
        for index, patch in enumerate(modeler.patches):
            element_type = modeler.patch_settings[index][
                "element_type"].GetString()
            number_of_nodes = 3 if element_type == "triangle" else 4
            ids = {node.Id for node in patch.Nodes}
            self.assertTrue(node_ids.isdisjoint(ids))
            node_ids.update(ids)
            ids = {element.Id for element in patch.Elements}
            self.assertTrue(element_ids.isdisjoint(ids))
            element_ids.update(ids)
            self.assertEqual(patch.NumberOfConditions(), 0)
            for element in patch.Elements:
                self.assertEqual(element.GetGeometry().PointsNumber(),
                                 number_of_nodes)
                self.assertGreater(element.GetGeometry().Area(), 0.0)
                if element.Is(KM.ACTIVE):
                    active_area += element.GetGeometry().Area()
        self.assertAlmostEqual(active_area, 24.0)
        self.assertEqual(modeler.model_part.NumberOfConditions(), 0)

        for index, patch_settings in enumerate(patches, 1):
            x, y = patch_settings["origin"]
            width, height = patch_settings["dimensions"]
            nx, ny = patch_settings["divisions"]
            patch = modeler.patches[index]
            factor = (2 if patch_settings.get("element_type", "triangle") ==
                      "triangle" else 1)
            self.assertEqual(patch.NumberOfElements(), factor * nx * ny)
            self.assertEqual(patch.NumberOfNodes(), (nx + 1) * (ny + 1))
            patch_area = sum(e.GetGeometry().Area() for e in patch.Elements)
            self.assertAlmostEqual(patch_area, width * height)
            inactive = background.GetSubModelPart(f"Inactive_Patch_{index}")
            self.assertEqual(inactive.NumberOfElements(),
                             round(background_factor * width * height /
                                   0.25 ** 2))
            inactive_area = sum(e.GetGeometry().Area()
                                for e in inactive.Elements)
            self.assertAlmostEqual(inactive_area, width * height)
            self.assertTrue(all(not e.Is(KM.ACTIVE)
                                for e in inactive.Elements))

            # On this face-aligned mesh an active element's bounding box can
            # only touch a hole, never intersect its interior with positive
            # area.
            for element in background.Elements:
                if not element.Is(KM.ACTIVE):
                    continue
                coordinates = [(n.X, n.Y) for n in element.GetGeometry()]
                overlap_x = (
                    min(max(p[0] for p in coordinates), x + width) -
                    max(min(p[0] for p in coordinates), x))
                overlap_y = (
                    min(max(p[1] for p in coordinates), y + height) -
                    max(min(p[1] for p in coordinates), y))
                self.assertTrue(overlap_x < 1e-12 or overlap_y < 1e-12)

            interface = model[f"FemMultipatch.Interfaces.Patch_{index}"]
            true = interface.GetSubModelPart("TrueBoundary")
            surrogate = interface.GetSubModelPart("SurrogateBoundary")
            self.assertEqual(true.NumberOfGeometries(), 4)
            self.assertEqual(surrogate.NumberOfGeometries(),
                             round(2 * (width + height) / 0.25))
            self.assertTrue({n.Id for n in true.Nodes}.isdisjoint(node_ids))
            self.assertTrue({n.Id for n in surrogate.Nodes}.issubset(
                {n.Id for n in background.Nodes}))
            for boundary in (true, surrogate):
                perimeter = sum(g.Length() for g in boundary.Geometries)
                self.assertAlmostEqual(perimeter, 2 * (width + height))
                degree = Counter(n.Id for g in boundary.Geometries for n in g)
                self.assertTrue(all(value == 2 for value in degree.values()))
                for geometry in boundary.Geometries:
                    sides = [(0, x), (0, x + width), (1, y), (1, y + height)]
                    self.assertTrue(any(
                        all(abs((n.X, n.Y)[axis] - value) < 1e-12
                            for n in geometry)
                        for axis, value in sides))

            # Verify every stored surrogate edge has exactly two neighbours:
            # one active triangle and one inactive triangle of this hole.
            for geometry in surrogate.Geometries:
                edge_ids = {node.Id for node in geometry}
                neighbours = [
                    e for e in background.Elements
                    if edge_ids.issubset({n.Id for n in e.GetGeometry()})]
                self.assertEqual(len(neighbours), 2)
                self.assertEqual(sum(e.Is(KM.ACTIVE) for e in neighbours), 1)
                removed = next(e for e in neighbours if not e.Is(KM.ACTIVE))
                self.assertTrue(removed.Id in inactive.Elements)

    def test_background_only(self):
        for element_type in ("triangle", "quadrilateral"):
            with self.subTest(element_type=element_type):
                self._check_partition([], element_type)

    def test_one_patch(self):
        for element_type in ("triangle", "quadrilateral"):
            with self.subTest(element_type=element_type):
                patches = self._patches()[:1]
                patches[0]["element_type"] = element_type
                self._check_partition(patches, element_type)

    def test_two_patches(self):
        self._check_partition(self._patches()[:2])

    def test_mixed_element_types(self):
        for background_type in ("triangle", "quadrilateral"):
            with self.subTest(background_type=background_type):
                patches = self._patches()[:2]
                patches[0]["element_type"] = "quadrilateral"
                patches[1]["element_type"] = "triangle"
                self._check_partition(patches, background_type)

    def test_arbitrary_patch_count(self):
        self._check_partition(self._patches())

    def test_actual_gauss_points(self):
        patches = self._patches()[:2]
        patches[0]["element_type"] = "quadrilateral"
        _, modeler = self._build(patches)
        for patch in modeler.patches:
            points = get_gauss_points(patch)
            grouped_points = {}
            for element_id, coordinates in points:
                grouped_points.setdefault(element_id, []).append(coordinates)
            for element in patch.Elements:
                if not element.Is(KM.ACTIVE):
                    self.assertNotIn(element.Id, grouped_points)
                    continue
                geometry = element.GetGeometry()
                triangle = geometry.PointsNumber() == 3
                method = (KM.GeometryData.IntegrationMethod.GI_GAUSS_1
                          if triangle else
                          KM.GeometryData.IntegrationMethod.GI_GAUSS_2)
                self.assertEqual(geometry.GetDefaultIntegrationMethod(),
                                 method)
                actual = grouped_points[element.Id]
                self.assertEqual(len(actual), 1 if triangle else 4)
                for point, reference in zip(
                        actual, element.GetIntegrationPoints()):
                    self.assertVectorAlmostEqual(point, reference)
                if triangle:
                    expected = [sum((n.X, n.Y, n.Z)[axis]
                                    for n in geometry) / 3.0
                                for axis in range(3)]
                    self.assertVectorAlmostEqual(actual[0], expected)
                else:
                    # Tensor product GI_GAUSS_2: two points per direction.
                    xmin = min(n.X for n in geometry)
                    xmax = max(n.X for n in geometry)
                    ymin = min(n.Y for n in geometry)
                    ymax = max(n.Y for n in geometry)
                    for point in actual:
                        self.assertAlmostEqual(
                            abs(point[0] - (xmin + xmax) / 2),
                            (xmax - xmin) / (2 * math.sqrt(3)))
                        self.assertAlmostEqual(
                            abs(point[1] - (ymin + ymax) / 2),
                            (ymax - ymin) / (2 * math.sqrt(3)))
        all_points = get_gauss_points(modeler.patches[0], active_only=False)
        self.assertEqual(len(all_points), 768)

    def test_reactivation_preserves_elements(self):
        _, modeler = self._build(self._patches()[:2])
        background = modeler.patches[0]
        utility = modeler.background_utility
        ids_before = [e.Id for e in background.Elements]
        retained = next(e for e in background.Elements if not e.Is(KM.ACTIVE))
        # Activation is independently reusable; no meshes or interfaces move.
        owners = utility.UpdateBackgroundActivation([], modeler.tolerance)
        self.assertEqual(owners, {})
        self.assertTrue(all(e.Is(KM.ACTIVE) for e in background.Elements))
        self.assertTrue(retained.Is(KM.ACTIVE))
        self.assertEqual(ids_before, [e.Id for e in background.Elements])
        self.assertEqual(utility.ExtractSurrogateEdges(owners), {})
        owners = utility.UpdateBackgroundActivation(
            modeler.boxes[1:], modeler.tolerance)
        self.assertFalse(retained.Is(KM.ACTIVE))
        self.assertEqual(set(utility.ExtractSurrogateEdges(owners)),
                         {1, 2})

    def test_nodes_are_independent(self):
        _, modeler = self._build(self._patches()[:1])
        background, internal = modeler.patches
        node = internal.GetNode(min(n.Id for n in internal.Nodes))
        coincident = next(
            n for n in background.Nodes
            if abs(n.X - node.X) < 1e-12 and abs(n.Y - node.Y) < 1e-12)
        original_x = coincident.X
        node.X += 0.01
        self.assertEqual(coincident.X, original_x)

    def test_invalid_boxes(self):
        invalid_patches = [
            ({"origin": [0.1, 0.5]}, "align"),
            ({"origin": [-1.0, 0.5]}, "strictly inside"),
            ({"origin": [-1.25, 0.5]}, "strictly inside"),
            ({"dimensions": [0.0, 1.0]}, "dimensions"),
            ({"divisions": [0, 4]}, "positive integers"),
            ({"origin": [0.0]}, "two components"),
            ({"element_type": "hexahedron"}, "element_type"),
        ]
        for change, message in invalid_patches:
            with self.subTest(change=change):
                patch = self._patches()[0].copy()
                patch.update(change)
                model = KM.Model()
                with self.assertRaisesRegex(ValueError, message):
                    FemMultipatchModeler(model, self._settings([patch]))
                self.assertFalse(model.HasModelPart("FemMultipatch"))

    def test_overlapping_or_touching_boxes(self):
        for origin in ([0.5, 0.5], [1.5, 0.5], [1.5, 2.5]):
            with self.subTest(origin=origin):
                patches = self._patches()[:2]
                patches[1]["origin"] = origin
                patches[1]["dimensions"] = [0.5, 0.5]
                with self.assertRaisesRegex(ValueError, "overlap or touch"):
                    FemMultipatchModeler(KM.Model(), self._settings(patches))

    def test_factory_and_unchanged_settings(self):
        model = KM.Model()
        settings = self._settings(self._patches()[:1])
        original = settings.WriteJsonString()
        item = KM.Parameters('''{
            "name": "",
            "parameters": {}
        }''')
        item["name"].SetString(
            "KratosMultiphysics.modelers.multipatch_modeler."
            "FemMultipatchModeler")
        item["parameters"] = settings
        modeler = KratosModelParametersFactory(model).ConstructItem(item)
        modeler.SetupModelPart()
        self.assertEqual(settings.WriteJsonString(), original)
        with self.assertRaisesRegex(RuntimeError, "only run once"):
            modeler.SetupModelPart()


if __name__ == "__main__":
    KratosUnittest.main()
