#     |  /           |
#     ' /   __| _` | __|  _ \   __|
#     . \  |   (   | |   (   |\__ `
#    _|\_\_|  \__,_|\__|\___/ ____/
#                    Multi-Physics
#
#   License:         BSD License
#                    Kratos default license: kratos/license.txt
#
#   Main author:     Nicolò Antonelli
#

"""Geometry and polynomial reproduction for conservative Q1 multipatch gaps."""
import json
import math
import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.modelers.multipatch_modeler import FemMultipatchModeler


def gap_settings(angle=23.0, divisions=12, origin=(0.33, 0.37), size=0.32):
    return KM.Parameters(json.dumps({
        "model_part_name": "Fluid", "interface_mode": "gap_sbm",
        "background": {"origin": [0.0, 0.0], "dimensions": [1.0, 1.0],
                       "divisions": [divisions, divisions], "element_type": "quadrilateral"},
        "internal_patches": [{"origin": list(origin), "dimensions": [size, size],
                              "divisions": [5, 7], "element_type": "quadrilateral",
                              "rotation_angle": angle}]}))


class TestFemMultipatchGapModeler(KratosUnittest.TestCase):
    def _build(self, settings):
        model = KM.Model()
        modeler = FemMultipatchModeler(model, settings)
        modeler.SetupModelPart()
        return model, modeler

    def test_partition_rotations_and_fitted_limit(self):
        cases = [(0.0, (0.33, 0.37), 0.32), (23.0, (0.33, 0.37), 0.32),
                 (45.0, (0.33, 0.37), 0.32), (137.0, (0.33, 0.37), 0.32), (0.0, (0.25, 0.25), 0.5)]
        for angle, origin, size in cases:
            with self.subTest(angle=angle, origin=origin):
                _, modeler = self._build(gap_settings(angle, origin=origin, size=size))
                bg, inner = modeler.patches
                root = modeler.model_part
                interface = root.GetSubModelPart("Interfaces.Patch_1")
                gap_area = sum(g.Area() for g in interface.GetSubModelPart("GapTriangles").Geometries)
                active_area = sum(e.GetGeometry().Area() for e in bg.Elements if e.Is(KM.ACTIVE))
                inner_area = sum(e.GetGeometry().Area() for e in inner.Elements)
                self.assertAlmostEqual(active_area + inner_area + gap_area, 1.0, delta=1e-10)
                length = sum(g.Length() for g in interface.GetSubModelPart("PhysicalInterfaces").Geometries)
                self.assertAlmostEqual(length, 4 * size, delta=1e-10)
                self.assertEqual(root.NumberOfNodes(), bg.NumberOfNodes() + inner.NumberOfNodes())
                self.assertEqual(root.NumberOfConditions(), 0)
                self.assertEqual(bg.NumberOfElements(), 144)
                inner_ids = {n.Id for n in inner.Nodes}
                self.assertTrue(inner_ids.isdisjoint(n.Id for n in bg.Nodes))
                # QP geometries carry exactly the original parent's nodes.
                for g in interface.GetSubModelPart("GapVolumes").Geometries:
                    self.assertEqual(len(g), 4)
                    self.assertTrue({n.Id for n in g}.issubset({n.Id for n in bg.Nodes}))
                if origin == (0.25, 0.25):
                    self.assertLess(gap_area, 1e-12)

    def test_coincident_edge_projection_roundoff(self):
        # At n=40 the right square edge x=0.65 is grid-aligned. Its
        # projected zero-area gap previously acquired signed area 1.4e-18.
        settings = gap_settings(angle=0.0, divisions=40)
        settings["internal_patches"][0]["divisions"].SetVector(KM.Vector([14, 14]))
        _, modeler = self._build(settings)
        bg, inner = modeler.patches
        interface = modeler.model_part.GetSubModelPart("Interfaces.Patch_1")
        area = sum(e.GetGeometry().Area() for e in bg.Elements if e.Is(KM.ACTIVE))
        area += sum(e.GetGeometry().Area() for e in inner.Elements)
        area += sum(g.Area() for g in interface.GetSubModelPart("GapTriangles").Geometries)
        self.assertAlmostEqual(area, 1.0, delta=1e-10)
        length = sum(g.Length() for g in interface.GetSubModelPart("PhysicalInterfaces").Geometries)
        self.assertAlmostEqual(length, 4 * 0.32, delta=1e-10)

    def test_small_corner_cut_near_grid_vertex(self):
        # At 47 degrees, cell 104 intersects the square in only 3.8e-12
        # area, but the penetration is about 2e-6. It must be inactive:
        # skipping it makes the adjacent projected gap polygons invalid.
        settings = gap_settings(angle=47.0, divisions=19)
        settings["internal_patches"][0]["divisions"].SetVector(KM.Vector([7, 7]))
        _, modeler = self._build(settings)
        bg, inner = modeler.patches
        self.assertFalse(bg.GetElement(104).Is(KM.ACTIVE))
        interface = modeler.model_part.GetSubModelPart("Interfaces.Patch_1")
        area = sum(e.GetGeometry().Area() for e in bg.Elements if e.Is(KM.ACTIVE))
        area += sum(e.GetGeometry().Area() for e in inner.Elements)
        area += sum(g.Area() for g in interface.GetSubModelPart("GapTriangles").Geometries)
        self.assertAlmostEqual(area, 1.0, delta=1e-10)
        length = sum(g.Length() for g in interface.GetSubModelPart("PhysicalInterfaces").Geometries)
        self.assertAlmostEqual(length, 4 * 0.32, delta=1e-10)

    def test_separated_patches(self):
        settings = gap_settings(size=0.17, origin=(0.19, 0.2))
        settings["internal_patches"].Append(KM.Parameters(json.dumps({
            "origin": [0.63, 0.61], "dimensions": [0.17, 0.17],
            "divisions": [3, 4], "element_type": "quadrilateral", "rotation_angle": -17.0})))
        _, modeler = self._build(settings)
        area = sum(e.GetGeometry().Area() for p in modeler.patches for e in p.Elements if e.Is(KM.ACTIVE))
        for i in (1, 2):
            area += sum(g.Area() for g in modeler.model_part.GetSubModelPart(
                f"Interfaces.Patch_{i}.GapTriangles").Geometries)
        self.assertAlmostEqual(area, 1.0, delta=1e-10)

    def test_reject_touching_surrogate_holes(self):
        settings = gap_settings(size=0.2, origin=(0.3, 0.3))
        settings["internal_patches"].Append(KM.Parameters(json.dumps({
            "origin": [0.52, 0.3], "dimensions": [0.2, 0.2],
            "divisions": [3, 4], "element_type": "quadrilateral"})))
        with self.assertRaisesRegex(RuntimeError, "overlap|touch"):
            self._build(settings)

    def test_reject_outer_boundary_and_unsupported_meshes(self):
        settings = gap_settings(origin=(0.0, 0.3))
        with self.assertRaisesRegex(ValueError, "strictly inside"):
            self._build(settings)
        settings = gap_settings()
        settings["background"]["element_type"].SetString("triangle")
        with self.assertRaisesRegex(ValueError, "Q1"):
            self._build(settings)


if __name__ == "__main__":
    KratosUnittest.main()
