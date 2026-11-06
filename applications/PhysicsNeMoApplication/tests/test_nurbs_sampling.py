import math

import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils

from KratosMultiphysics.PhysicsNeMoApplication.bridges.mesh_bridge import nurbs_sampling

try:
    import torch  # noqa: F401
    import physicsnemo.mesh  # noqa: F401
    have_physicsnemo_mesh = True
except ImportError:
    have_physicsnemo_mesh = False


def _MakeNodesVector(model_part, coordinates, first_id=1):
    nodes = Kratos.NodesVector()
    for offset, (x, y, z) in enumerate(coordinates):
        nodes.append(model_part.CreateNewNode(first_id + offset, float(x), float(y), float(z)))
    return nodes


def CreateFlatSurface(model_part):
    """Degree-1 unit square in the z = 0 plane (2x2 control points)."""
    nodes = _MakeNodesVector(
        model_part, [(0, 0, 0), (1, 0, 0), (0, 1, 0), (1, 1, 0)])
    knots = Kratos.Vector(2)
    knots[0], knots[1] = 0.0, 1.0
    return Kratos.NurbsSurfaceGeometry3D(nodes, 1, 1, knots, knots)


def CreateQuarterCylinderSurface(model_part, first_id=1):
    """Rational quadratic quarter circle (radius 1) extruded to z = 1.

    The classical exact-arc NURBS: weights (1, sqrt(2)/2, 1) put every
    evaluated point exactly on the circle - which is the property the
    sampling tests pin.
    """
    arc = [(1, 0), (1, 1), (0, 1)]
    nodes = _MakeNodesVector(
        model_part,
        [(x, y, z) for z in (0.0, 1.0) for x, y in arc],
        first_id=first_id)
    knots_u = Kratos.Vector(4)
    knots_u[0], knots_u[1], knots_u[2], knots_u[3] = 0.0, 0.0, 1.0, 1.0
    knots_v = Kratos.Vector(2)
    knots_v[0], knots_v[1] = 0.0, 1.0
    weights = Kratos.Vector(6)
    for i, w in enumerate([1.0, math.sqrt(2.0) / 2.0, 1.0] * 2):
        weights[i] = w
    return Kratos.NurbsSurfaceGeometry3D(nodes, 2, 1, knots_u, knots_v, weights)


def CreateUnitCubeVolume(model_part):
    """Degree-1 unit cube (2x2x2 control points)."""
    corners = [(x, y, z) for z in (0.0, 1.0) for y in (0.0, 1.0) for x in (0.0, 1.0)]
    nodes = _MakeNodesVector(model_part, corners)
    knots = Kratos.Vector(2)
    knots[0], knots[1] = 0.0, 1.0
    return Kratos.NurbsVolumeGeometry(nodes, 1, 1, 1, knots, knots, knots)


def _TetVolumes(points, cells):
    a, b, c, d = (points[cells[:, i]] for i in range(4))
    return numpy.einsum("ij,ij->i", b - a, numpy.cross(c - a, d - a)) / 6.0


def _BoundaryFaceCounts(cells):
    counts = {}
    for tet in cells:
        for corners in ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)):
            key = tuple(sorted(int(tet[i]) for i in corners))
            counts[key] = counts.get(key, 0) + 1
    return counts


class TestNurbsSampling(KratosUnittest.TestCase):
    """Exact NURBS sampling against core geometries - no IgaApplication, no
    torch: the geometry evaluation and the isogeometric gather are pure
    Kratos + numpy."""

    def setUp(self):
        self.model = Kratos.Model()
        self.model_part = self.model.CreateModelPart("Main")
        self.model_part.AddNodalSolutionStepVariable(Kratos.TEMPERATURE)
        self.model_part.AddNodalSolutionStepVariable(Kratos.VELOCITY)

    def test_FlatSurfaceIsSampledExactly(self):
        surface = CreateFlatSurface(self.model_part)
        points, cells, local_coordinates = nurbs_sampling.SampleNurbsGeometry(surface, 4)
        self.assertEqual(len(points), 25)
        self.assertEqual(len(cells), 32)
        # the degree-1 patch IS the unit square: samples land on the lattice
        for point, local in zip(points, local_coordinates):
            self.assertVectorAlmostEqual(point, [local[0], local[1], 0.0], places=12)
        # the two triangles per cell tile the square exactly
        a, b, c = (points[cells[:, i]] for i in range(3))
        area = float(numpy.linalg.norm(numpy.cross(b - a, c - a), axis=1).sum() / 2.0)
        self.assertAlmostEqual(area, 1.0, places=12)

    def test_RationalArcSamplesLieExactlyOnTheCircle(self):
        # the property a chordal/linear approximation cannot have: EVERY
        # sampled point of the rational patch is on the radius-1 cylinder
        surface = CreateQuarterCylinderSurface(self.model_part)
        self.assertEqual(surface.NumberOfControlPointsU(), 3)
        points, _, _ = nurbs_sampling.SampleNurbsGeometry(surface, [16, 2])
        radii = numpy.hypot(points[:, 0], points[:, 1])
        self.assertAlmostEqual(float(numpy.abs(radii - 1.0).max()), 0.0, places=12)
        self.assertAlmostEqual(float(points[:, 2].min()), 0.0, places=12)
        self.assertAlmostEqual(float(points[:, 2].max()), 1.0, places=12)

    def test_VolumeIsWatertightAndPositivelyOriented(self):
        volume = CreateUnitCubeVolume(self.model_part)
        points, cells, _ = nurbs_sampling.SampleNurbsGeometry(volume, 3)
        self.assertEqual(len(points), 4 ** 3)
        self.assertEqual(len(cells), 6 * 3 ** 3)
        volumes = _TetVolumes(points, cells)
        self.assertGreater(float(volumes.min()), 0.0)
        self.assertAlmostEqual(float(volumes.sum()), 1.0, places=12)
        # watertight: every interior face shared by exactly two tetrahedra,
        # boundary faces by one, and the boundary area is the cube's
        counts = _BoundaryFaceCounts(cells)
        self.assertLessEqual(set(counts.values()), {1, 2})
        boundary = numpy.array([face for face, n in counts.items() if n == 1])
        a, b, c = (points[boundary[:, i]] for i in range(3))
        area = float(numpy.linalg.norm(numpy.cross(b - a, c - a), axis=1).sum() / 2.0)
        self.assertAlmostEqual(area, 6.0, places=12)

    def test_IsogeometricGatherReproducesAffineFields(self):
        surface = CreateFlatSurface(self.model_part)
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(Kratos.TEMPERATURE, 1.0 + 2.0 * node.X + 3.0 * node.Y)
            node.SetSolutionStepValue(Kratos.VELOCITY, [node.X, -node.Y, 0.5])
        points, _, local_coordinates = nurbs_sampling.SampleNurbsGeometry(surface, 3)

        temperature = nurbs_sampling.EvaluateNodalFieldOnLattice(
            surface, local_coordinates, Kratos.TEMPERATURE)
        self.assertEqual(temperature.shape, (len(points),))
        self.assertVectorAlmostEqual(
            temperature, 1.0 + 2.0 * points[:, 0] + 3.0 * points[:, 1], places=12)

        velocity = nurbs_sampling.EvaluateNodalFieldOnLattice(
            surface, local_coordinates, Kratos.VELOCITY)
        self.assertEqual(velocity.shape, (len(points), 3))
        self.assertVectorAlmostEqual(velocity[:, 0], points[:, 0], places=12)
        self.assertVectorAlmostEqual(velocity[:, 1], -points[:, 1], places=12)

    def test_GatherOnTheRationalPatchIsTheGeometrys(self):
        # constants must come back exactly - the rational basis is a
        # partition of unity, and the basis used is the C++ geometry's own
        surface = CreateQuarterCylinderSurface(self.model_part)
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(Kratos.TEMPERATURE, 4.25)
        _, _, local_coordinates = nurbs_sampling.SampleNurbsGeometry(surface, [8, 2])
        temperature = nurbs_sampling.EvaluateNodalFieldOnLattice(
            surface, local_coordinates, Kratos.TEMPERATURE)
        self.assertAlmostEqual(float(numpy.abs(temperature - 4.25).max()), 0.0, places=12)

    def test_ModelPartGeometriesAreFound(self):
        if not have_physicsnemo_mesh:
            self.skipTest("Missing required python modules: torch, physicsnemo.")
        surface = CreateQuarterCylinderSurface(self.model_part)
        self.model_part.AddGeometry(surface)
        meshes = nurbs_sampling.SampleModelPartNurbsGeometries(self.model_part, 4)
        self.assertEqual(len(meshes), 1)
        mesh = next(iter(meshes.values()))
        self.assertEqual(int(mesh.cells.shape[1]), 3)

        empty = self.model.CreateModelPart("Empty")
        with self.assertRaisesRegex(ValueError, "carries no"):
            nurbs_sampling.SampleModelPartNurbsGeometries(empty, 4)

    def test_Validation(self):
        surface = CreateFlatSurface(self.model_part)
        with self.assertRaisesRegex(ValueError, "divisions"):
            nurbs_sampling.SampleNurbsGeometry(surface, [4])
        with self.assertRaisesRegex(ValueError, "divisions"):
            nurbs_sampling.SampleNurbsGeometry(surface, 0)
        with self.assertRaisesRegex(TypeError, "Nurbs"):
            nurbs_sampling.NurbsParametricLattice("not a geometry", 4)


@KratosUnittest.skipUnless(kratos_utils.CheckIfApplicationsAvailable("IgaApplication"),
                           "IgaApplication is not available.")
class TestNurbsSamplingWithIgaModeler(KratosUnittest.TestCase):
    """The modeler workflow: NurbsGeometryModeler builds the geometry into
    ModelPart.Geometries, and the sampler finds it there - the actual
    IgaApplication route, exercised when the application is compiled."""

    def test_ModelerBuiltSurfaceIsSampled(self):
        import KratosMultiphysics.IgaApplication  # noqa: F401 - registers the modeler
        from KratosMultiphysics.modeler_factory import KratosModelerFactory

        model = Kratos.Model()
        modelers = Kratos.Parameters("""[{
            "modeler_name": "NurbsGeometryModeler",
            "Parameters": {
                "model_part_name"      : "Mesh",
                "lower_point_xyz"      : [0.0, 0.0, 0.0],
                "upper_point_xyz"      : [2.0, 1.0, 0.0],
                "lower_point_uvw"      : [0.0, 0.0, 0.0],
                "upper_point_uvw"      : [1.0, 1.0, 0.0],
                "polynomial_order"     : [3, 2],
                "number_of_knot_spans" : [3, 2]
            }
        }]""")
        for modeler in KratosModelerFactory().ConstructListOfModelers(model, modelers):
            modeler.SetupGeometryModel()
            modeler.PrepareGeometryModel()
            modeler.SetupModelPart()

        model_part = model.GetModelPart("Mesh")
        geometries = [geometry for geometry in model_part.Geometries
                      if nurbs_sampling._IsNurbsSurface(geometry)]
        self.assertEqual(len(geometries), 1)
        points, cells, _ = nurbs_sampling.SampleNurbsGeometry(geometries[0], [6, 3])
        # the modeler's patch IS the 2 x 1 rectangle: exact corners, flat plane
        self.assertAlmostEqual(float(numpy.abs(points[:, 2]).max()), 0.0, places=12)
        self.assertAlmostEqual(float(points[:, 0].min()), 0.0, places=12)
        self.assertAlmostEqual(float(points[:, 0].max()), 2.0, places=12)
        self.assertAlmostEqual(float(points[:, 1].max()), 1.0, places=12)
        a, b, c = (points[cells[:, i]] for i in range(3))
        area = float(numpy.linalg.norm(numpy.cross(b - a, c - a), axis=1).sum() / 2.0)
        self.assertAlmostEqual(area, 2.0, places=12)


@KratosUnittest.skipUnless(have_physicsnemo_mesh,
                           "Missing required python modules: torch, physicsnemo.")
class TestNurbsMeshBuild(KratosUnittest.TestCase):
    """BuildNurbsMesh: the sampled geometry as a physicsnemo Mesh with the
    isogeometric gather attached as point_data."""

    def test_MeshCarriesGatheredPointData(self):
        model = Kratos.Model()
        model_part = model.CreateModelPart("Main")
        model_part.AddNodalSolutionStepVariable(Kratos.TEMPERATURE)
        surface = CreateFlatSurface(model_part)
        for node in model_part.Nodes:
            node.SetSolutionStepValue(Kratos.TEMPERATURE, 2.0 * node.X + node.Y)

        mesh = nurbs_sampling.BuildNurbsMesh(
            surface, 3, field_specs=[(Kratos.TEMPERATURE, "node_historical")])
        points = mesh.points.detach().cpu().numpy()
        values = mesh.point_data["TEMPERATURE"].detach().cpu().numpy()
        self.assertVectorAlmostEqual(values, 2.0 * points[:, 0] + points[:, 1], places=12)

        with self.assertRaisesRegex(ValueError, "nodal"):
            nurbs_sampling.BuildNurbsMesh(
                surface, 2, field_specs=[(Kratos.TEMPERATURE, "element")])


if __name__ == '__main__':
    KratosUnittest.main()
