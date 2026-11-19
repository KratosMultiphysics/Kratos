import json
import math
import os
import tempfile

import numpy as np

import KratosMultiphysics as KM
import KratosMultiphysics.GeoMechanicsApplication  # noqa: F401 (registers the Geo variables)
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.GeoMechanicsApplication import bishop_slope_stability as bishop
from KratosMultiphysics.GeoMechanicsApplication import bishop_slope_stability_process

# Homogeneous slope: toe at (4, 4), crest at (10, 10), base at y = 0, sides at x = 0 and x = 20
SURFACE = np.array([[0.0, 4.0], [4.0, 4.0], [10.0, 10.0], [20.0, 10.0]])
UNIT_WEIGHT = 20.0e3
CIRCLE = ((5.0, 13.0), 9.5)


def surface_y(x):
    return np.interp(x, SURFACE[:, 0], SURFACE[:, 1])


def structured_slope_mesh(columns_per_meter, rows):
    """Triangles of a structured mesh whose top boundary is exactly SURFACE."""
    xs = np.unique(np.concatenate([np.linspace(0.0, 20.0, int(20 * columns_per_meter) + 1), SURFACE[:, 0]]))
    coordinates = np.array([[x, y] for x in xs for y in np.linspace(0.0, surface_y(x), rows + 1)])
    triangles = []
    for i in range(len(xs) - 1):
        for j in range(rows):
            a, b = i * (rows + 1) + j, (i + 1) * (rows + 1) + j
            triangles += [[a, b, b + 1], [a, b + 1, a + 1]]
    forbidden = (coordinates[:, 1] < 1.0e-12) | (coordinates[:, 0] < 1.0e-12) | (coordinates[:, 0] > 20.0 - 1.0e-12)
    return coordinates, np.array(triangles), forbidden


def sliding_mass_polygon(center, radius, entry, exit_, points_on_arc=20000):
    """Exact geometry of the sliding mass: surface between entry and exit closed by the arc."""
    surface = [entry] + [p for p in SURFACE if entry[0] < p[0] < exit_[0]] + [exit_]
    # Lower branch from the exit back to the entry (both angles in (-pi, 0))
    theta = np.linspace(math.atan2(exit_[1] - center[1], exit_[0] - center[0]),
                        math.atan2(entry[1] - center[1], entry[0] - center[0]), points_on_arc)
    arc = np.column_stack([center[0] + radius * np.cos(theta), center[1] + radius * np.sin(theta)])
    return np.vstack([np.array(surface), arc[1:-1]])


def polygon_area_and_centroid_x(polygon):
    x, y = polygon[:, 0], polygon[:, 1]
    x1, y1 = np.roll(x, -1), np.roll(y, -1)
    cross = x * y1 - x1 * y
    area = 0.5 * np.sum(cross)
    return abs(area), np.sum((x + x1) * cross) / (6.0 * area)


def constant_base(cohesion, friction_angle_deg):
    def base_properties(x_mid, y_mid):
        return np.full(len(x_mid), cohesion), np.full(len(x_mid), math.radians(friction_angle_deg))
    return base_properties


def meshless_column_weight(xs, y_base):
    return UNIT_WEIGHT * np.clip(surface_y(xs) - y_base, 0.0, None)


class KratosGeoMechanicsBishopSlopeStabilityTests(KratosUnittest.TestCase):
    def evaluate_on_mesh(self, mesh, forbidden, base_properties, number_of_slices=30, circle=CIRCLE):
        gamma = np.full(len(mesh.triangles), UNIT_WEIGHT)
        material = np.zeros(len(mesh.triangles), dtype=int)
        surface_breaks = bishop.surface_breakpoints(mesh, forbidden)

        def extra_breakpoints(arc):
            return np.concatenate([surface_breaks, bishop.base_breakpoints(mesh, arc, material)])

        return bishop.evaluate_circle(mesh, circle[0], circle[1], number_of_slices,
                                      lambda xs, yb: bishop.column_weights(mesh, gamma, xs, yb), base_properties,
                                      forbidden_nodes=forbidden, extra_breakpoints=extra_breakpoints)

    def test_cohesive_soil_matches_moment_equilibrium(self):
        # For phi = 0 Bishop reduces to FoS = c R^2 dtheta / M, with M the moment of the weight about the center
        coordinates, triangles, forbidden = structured_slope_mesh(2, 8)
        mesh = bishop.build_triangle_mesh(coordinates, triangles)
        result, slices, arc = self.evaluate_on_mesh(mesh, forbidden, constant_base(20.0e3, 0.0), number_of_slices=200)

        area, centroid_x = polygon_area_and_centroid_x(sliding_mass_polygon(arc.center, arc.radius, arc.entry, arc.exit))
        moment = UNIT_WEIGHT * area * abs(centroid_x - arc.center[0])
        theta = [math.atan2(p[1] - arc.center[1], p[0] - arc.center[0]) for p in (arc.entry, arc.exit)]
        expected = 20.0e3 * arc.radius ** 2 * abs(theta[1] - theta[0]) / moment
        self.assertAlmostEqual(result.factor_of_safety, expected, delta=1.0e-4 * expected)
        self.assertAlmostEqual(result.factor_of_safety, result.fellenius_factor_of_safety, delta=1.0e-3 * expected)

    def test_slice_weights_on_mesh_equal_meshless_weights_and_total_weight(self):
        coordinates, triangles, forbidden = structured_slope_mesh(2, 8)
        mesh = bishop.build_triangle_mesh(coordinates, triangles)
        base = constant_base(20.0e3, 30.0)
        result, slices, arc = self.evaluate_on_mesh(mesh, forbidden, base)

        meshless = bishop.build_slices(arc.center, arc.radius, np.append(slices.x_left, slices.x_right[-1]),
                                       meshless_column_weight, base)
        np.testing.assert_allclose(slices.weight, meshless.weight, rtol=1.0e-10)
        self.assertAlmostEqual(result.factor_of_safety, bishop.bishop_factor_of_safety(meshless).factor_of_safety,
                               delta=1.0e-10)

        area, _ = polygon_area_and_centroid_x(sliding_mass_polygon(arc.center, arc.radius, arc.entry, arc.exit))
        self.assertAlmostEqual(np.sum(slices.weight), UNIT_WEIGHT * area, delta=1.0e-6 * UNIT_WEIGHT * area)

    def test_factor_of_safety_does_not_depend_on_the_mesh(self):
        factors = []
        for columns_per_meter, rows in [(1, 4), (4, 16)]:
            coordinates, triangles, forbidden = structured_slope_mesh(columns_per_meter, rows)
            mesh = bishop.build_triangle_mesh(coordinates, triangles)
            gamma = np.full(len(mesh.triangles), UNIT_WEIGHT)
            surface_breaks = bishop.surface_breakpoints(mesh, forbidden)

            def evaluate(xc, yc, radius):
                evaluation = bishop.evaluate_circle(
                    mesh, (xc, yc), radius, 30, lambda xs, yb: bishop.column_weights(mesh, gamma, xs, yb),
                    constant_base(20.0e3, 30.0), forbidden_nodes=forbidden, extra_breakpoints=lambda arc: surface_breaks)
                return None if evaluation is None else evaluation[0].factor_of_safety

            factors.append(bishop.grid_search(evaluate, np.linspace(2.0, 8.0, 7), np.linspace(10.0, 16.0, 7),
                                              np.linspace(5.0, 12.0, 8))[0])
        self.assertAlmostEqual(factors[0], factors[1], delta=0.005 * factors[1])



    def test_inadmissible_circles_are_discarded(self):
        coordinates, triangles, forbidden = structured_slope_mesh(2, 8)
        mesh = bishop.build_triangle_mesh(coordinates, triangles)
        self.assertIsNone(bishop.find_slip_arc(mesh, (5.0, 30.0), 2.0, forbidden))    # does not cut the mesh
        self.assertIsNone(bishop.find_slip_arc(mesh, (5.0, 13.0), 12.9, forbidden))   # cuts the base
        self.assertIsNone(bishop.find_slip_arc(mesh, (7.0, 4.0), 3.0, forbidden))     # crossings above the center
        self.assertIsNotNone(bishop.find_slip_arc(mesh, CIRCLE[0], CIRCLE[1], forbidden))

    def test_process_on_a_model_part(self):
        coordinates, triangles, forbidden = structured_slope_mesh(2, 8)
        model = KM.Model()
        model_part = model.CreateModelPart("Slope")
        model_part.AddNodalSolutionStepVariable(KM.DISTANCE)
        for i, (x, y) in enumerate(coordinates, 1):
            model_part.CreateNewNode(i, x, y, 0.0)
        properties = model_part.GetProperties()[1]
        properties.SetValue(KM.KratosGlobals.GetVariable("GEO_COHESION"), 20.0e3)
        properties.SetValue(KM.KratosGlobals.GetVariable("GEO_FRICTION_ANGLE"), 30.0)
        properties.SetValue(KM.KratosGlobals.GetVariable("DENSITY_SOLID"), UNIT_WEIGHT / 9.81)
        for i, triangle in enumerate(triangles, 1):
            model_part.CreateNewElement("Element2D3N", i, [int(n) + 1 for n in triangle], properties)
        boundary = model_part.CreateSubModelPart("Base_And_Sides")
        boundary.AddNodes([i + 1 for i in np.flatnonzero(forbidden)])

        with tempfile.TemporaryDirectory() as directory:
            output_file = os.path.join(directory, "bishop.json")
            settings = KM.Parameters(json.dumps({"Parameters": {
                "model_part_name": "Slope",
                "number_of_slices": 30,
                "search": {"center_x_range": [CIRCLE[0][0], CIRCLE[0][0]], "n_centers_x": 1,
                           "center_y_range": [CIRCLE[0][1], CIRCLE[0][1]], "n_centers_y": 1,
                           "radius_range": [CIRCLE[1], CIRCLE[1]], "n_radii": 1},
                "forbidden_boundary_model_parts": ["Slope.Base_And_Sides"],
                "output_file": output_file}}))
            process = bishop_slope_stability_process.Factory(settings, model)
            process.ExecuteFinalizeSolutionStep()
            with open(output_file) as f:
                output = json.load(f)

        mesh = bishop.build_triangle_mesh(coordinates, triangles)
        expected, slices, arc = self.evaluate_on_mesh(mesh, forbidden, constant_base(20.0e3, 30.0))
        self.assertAlmostEqual(output["factor_of_safety"], expected.factor_of_safety, delta=1.0e-6)
        np.testing.assert_allclose([s["W"] for s in output["slices"]], slices.weight, rtol=1.0e-9)
        self.assertEqual(output["sliding_direction"], "negative x")

        for node in model_part.Nodes:
            expected_distance = math.hypot(node.X - arc.center[0], node.Y - arc.center[1]) - arc.radius
            self.assertAlmostEqual(node.GetSolutionStepValue(KM.DISTANCE), expected_distance)


if __name__ == "__main__":
    KratosUnittest.main()
