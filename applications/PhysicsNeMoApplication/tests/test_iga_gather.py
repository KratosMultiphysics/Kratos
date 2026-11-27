"""A real isogeometric solve feeding the isogeometric gather.

`nurbs_sampling` could already sample an exact NURBS geometry and gather
control-point fields through the geometry's own basis, but it had only ever
been fed geometries built by hand. An IGA ANALYSIS produces Brep wrappers
instead, which the sampler used to reject on type - so the one thing the
bridge was written for, real isogeometric results, was the one thing it
could not consume.
"""

import sys
from pathlib import Path

import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils

_CASES_DIR = Path(__file__).parent / "kratos_solver_cases"
sys.path.insert(0, str(_CASES_DIR))

have_iga = kratos_utils.CheckIfApplicationsAvailable(
    "IgaApplication", "StructuralMechanicsApplication", "LinearSolversApplication")


@KratosUnittest.skipUnless(have_iga, "Missing Iga/StructuralMechanics/LinearSolvers.")
class TestIgaAnalysisGather(KratosUnittest.TestCase):
    def setUp(self):
        import iga_case

        if not iga_case.IsAvailable():
            self.skipTest("IgaApplication test data is not present.")
        self.model = Kratos.Model()
        self.control_displacements, self.model_part = iga_case.SolveIga(self.model)
        self.surface = [geometry for geometry in self.model_part.Geometries][0]

    def test_TheAnalysisSolves(self):
        """Not a given: the application's own single_patch case is
        structurally singular and returns NaN through its configured
        iterative solver without complaining, which is why this case is the
        Scordelis roof instead."""
        self.assertEqual(self.model_part.NumberOfNodes(), 36)   # control points
        self.assertTrue(numpy.isfinite(self.control_displacements).all())
        self.assertGreater(numpy.abs(self.control_displacements).max(), 0.0)

    def test_TheSolveesGeometryIsSampledExactly(self):
        from KratosMultiphysics.PhysicsNeMoApplication.bridges.mesh_bridge import (
            nurbs_sampling)

        points, cells, local_coordinates = nurbs_sampling.SampleNurbsGeometry(
            self.surface, 3)
        points = numpy.asarray(points)
        self.assertEqual(points.shape[1], 3)
        self.assertEqual(numpy.asarray(cells).shape[1], 3)      # triangles
        self.assertTrue(numpy.isfinite(points).all())
        # every sampled point is the geometry's own evaluation, so it lies
        # on the surface rather than on a facet through its control points
        for row, local in enumerate(local_coordinates[:5]):
            numpy.testing.assert_allclose(
                points[row],
                numpy.array(self.surface.GlobalCoordinates(Kratos.Vector(list(local)))),
                atol=1e-10)

    def test_TheGatherRunsThroughTheNurbsBasis(self):
        """The control points are not on the surface, so the field at a
        surface point is a basis-weighted combination of control values -
        which is why this is a gather and not a lookup."""
        from KratosMultiphysics.PhysicsNeMoApplication.bridges.mesh_bridge import (
            nurbs_sampling)

        local_coordinates, _ = nurbs_sampling.NurbsParametricLattice(self.surface, 3)
        gathered = numpy.asarray(nurbs_sampling.EvaluateNodalFieldOnLattice(
            self.surface, local_coordinates, Kratos.DISPLACEMENT))

        self.assertEqual(gathered.shape, (len(local_coordinates), 3))
        self.assertTrue(numpy.isfinite(gathered).all())
        self.assertGreater(numpy.abs(gathered).max(), 0.0)
        # a B-spline evaluation lies in the convex hull of its control
        # values, so the surface field cannot exceed the control field
        self.assertLessEqual(numpy.abs(gathered).max(),
                             numpy.abs(self.control_displacements).max() + 1e-9)

    def test_TheSamplerAcceptsWhatAnAnalysisProduces(self):
        """The type check used to demand a concrete NurbsSurfaceGeometry3D,
        so geometries built by hand were sampleable and the ones a solve
        produces were not."""
        from KratosMultiphysics.PhysicsNeMoApplication.bridges.mesh_bridge import (
            nurbs_sampling)

        self.assertNotIsInstance(self.surface, Kratos.NurbsSurfaceGeometry3D)
        self.assertTrue(nurbs_sampling._IsNurbsSurface(self.surface))
        self.assertFalse(nurbs_sampling._IsNurbsVolume(self.surface))

    def test_AMeshIsBuiltFromTheSolvedGeometry(self):
        from KratosMultiphysics.PhysicsNeMoApplication.bridges.mesh_bridge import (
            nurbs_sampling)

        try:
            mesh = nurbs_sampling.BuildNurbsMesh(self.surface, 3)
        except ImportError:
            self.skipTest("BuildNurbsMesh needs torch and physicsnemo.")
        self.assertEqual(int(mesh.points.shape[1]), 3)
        self.assertGreater(int(mesh.cells.shape[0]), 0)


if __name__ == '__main__':
    KratosUnittest.main()
