"""Tests for the geometry guardrail: out-of-distribution detection on the
SHAPE rather than on the field values."""

from pathlib import Path

import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils

from KratosMultiphysics.PhysicsNeMoApplication.deployment import geometry_guard_utils

from test_grid_bridge import CreateStructuredTetModelPart

try:
    import torch
    have_torch = True
except ImportError:
    have_torch = False

try:
    import warnings as _warnings
    with _warnings.catch_warnings():
        _warnings.simplefilter("ignore")
        from physicsnemo.experimental.guardrails.geometry import (  # noqa: F401
            GeometryGuardrail)
    have_guardrail = True
except ImportError:
    have_guardrail = False


def _ScaledBox(model, name, scale, divisions=3):
    """A unit tet cube scaled per axis - a one-parameter geometry family."""
    model_part = CreateStructuredTetModelPart(
        model, name, divisions=divisions, historical_variables=(Kratos.PRESSURE,))
    for node in model_part.Nodes:
        node.X0 *= scale[0]
        node.Y0 *= scale[1]
        node.Z0 *= scale[2]
        node.X, node.Y, node.Z = node.X0, node.Y0, node.Z0
    return model_part


def _Family(model, count=24, spread=0.15, seed=0, prefix="family"):
    rng = numpy.random.default_rng(seed)
    return [_ScaledBox(model, f"{prefix}_{seed}_{i}", 1.0 + spread * rng.standard_normal(3))
            for i in range(count)]


@KratosUnittest.skipUnless(have_torch and have_guardrail,
                           "Missing required python modules: torch, physicsnemo (>= 2.2).")
class TestGeometryGuardUtils(KratosUnittest.TestCase):
    def setUp(self):
        self.model = Kratos.Model()
        self.sidecar = Path("test_geometry_guard.npz")

    def tearDown(self):
        kratos_utils.DeleteFileIfExisting(str(self.sidecar))

    def test_SurfaceOfModelPartIsAClosedTriangleSurface(self):
        surface = geometry_guard_utils.SurfaceOfModelPart(
            _ScaledBox(self.model, "one", (1.0, 1.0, 1.0)))
        self.assertEqual(surface.cells.shape[1], 3)   # triangles
        self.assertGreater(surface.points.shape[0], 0)
        # every boundary facet of a cube, and no interior ones
        self.assertEqual(surface.cells.shape[0], 6 * 2 * 3 * 3)

    def test_TheFamilyIsAcceptedAndAnAlienShapeIsRejected(self):
        family = _Family(self.model)
        guard = geometry_guard_utils.CreateGeometryGuard(Kratos.Parameters("{}"))
        geometry_guard_utils.FitGeometryGuard(guard, family)

        near = geometry_guard_utils.QueryGeometry(guard, family[0])
        self.assertEqual(near["status"], "OK")

        slab = geometry_guard_utils.QueryGeometry(
            guard, _ScaledBox(self.model, "slab", (5.0, 1.0, 0.2)))
        self.assertEqual(slab["status"], "REJECT")
        self.assertGreater(slab["percentile"], near["percentile"])

    def test_AShapeMoreRegularThanTheFamilyIsAlsoOutOfDistribution(self):
        """Worth knowing before trusting a verdict: the guardrail measures
        distance in DESCRIPTOR space, not "is this a sensible shape". A
        family of randomly perturbed boxes rejects the perfect cube, whose
        exactly equal side lengths no training member has."""
        family = _Family(self.model)
        guard = geometry_guard_utils.CreateGeometryGuard(Kratos.Parameters("{}"))
        geometry_guard_utils.FitGeometryGuard(guard, family)
        perfect = geometry_guard_utils.QueryGeometry(
            guard, _ScaledBox(self.model, "perfect", (1.0, 1.0, 1.0)))
        self.assertEqual(perfect["status"], "REJECT")

    def test_ASidecarMustBeNamedNpz(self):
        family = _Family(self.model, count=4)
        guard = geometry_guard_utils.CreateGeometryGuard(Kratos.Parameters("{}"))
        geometry_guard_utils.FitGeometryGuard(guard, family)
        with self.assertRaisesRegex(ValueError, "npz"):
            geometry_guard_utils.SaveGeometryGuard(guard, "guard.pt")

    def test_TheFieldGuardCannotSeeThis(self):
        """The point of the shape guardrail: a geometry far from the
        training family whose nodal VALUES are perfectly ordinary."""
        family = _Family(self.model)
        guard = geometry_guard_utils.CreateGeometryGuard(Kratos.Parameters("{}"))
        geometry_guard_utils.FitGeometryGuard(guard, family)

        slab = _ScaledBox(self.model, "ordinary_values", (5.0, 1.0, 0.2))
        for node in slab.Nodes:
            node.SetSolutionStepValue(Kratos.PRESSURE, 0.5)  # in range for any field guard
        self.assertEqual(
            geometry_guard_utils.QueryGeometry(guard, slab)["status"], "REJECT")

    def test_SaveLoadRoundTripKeepsTheVerdict(self):
        family = _Family(self.model)
        guard = geometry_guard_utils.CreateGeometryGuard(Kratos.Parameters("{}"))
        geometry_guard_utils.FitGeometryGuard(guard, family)
        slab = _ScaledBox(self.model, "slab", (5.0, 1.0, 0.2))
        before = geometry_guard_utils.QueryGeometry(guard, slab)

        geometry_guard_utils.SaveGeometryGuard(guard, self.sidecar)
        self.assertTrue(self.sidecar.exists())
        reloaded = geometry_guard_utils.LoadGeometryGuard(self.sidecar)
        after = geometry_guard_utils.QueryGeometry(reloaded, slab)
        self.assertEqual(after["status"], before["status"])
        self.assertAlmostEqual(after["percentile"], before["percentile"], places=6)

    def test_ASingleGeometryIsNotAFamily(self):
        with self.assertRaisesRegex(ValueError, "FAMILY"):
            geometry_guard_utils.FitGeometryGuard(
                geometry_guard_utils.CreateGeometryGuard(Kratos.Parameters("{}")),
                [_ScaledBox(self.model, "lonely", (1.0, 1.0, 1.0))])

    def test_AnUnderDeterminedFitIsReported(self):
        """Upstream's descriptor vector is 22-wide, and a GMM fitted on
        fewer geometries is under-determined - its verdicts stop meaning
        anything, with no complaint from upstream. Only the size check can
        warn about it; the resulting verdicts are too unstable to assert."""
        family = _Family(self.model, count=6, spread=0.05, seed=1)
        meshes = [geometry_guard_utils.SurfaceOfModelPart(part) for part in family]
        self.assertEqual(geometry_guard_utils.FeatureWidth(meshes[0]), 22)
        complaint = geometry_guard_utils.CheckFamilySize(meshes)
        self.assertIsNotNone(complaint)
        self.assertIn("under-determined", complaint)

        # a family large enough draws no complaint
        self.assertIsNone(geometry_guard_utils.CheckFamilySize(
            [geometry_guard_utils.SurfaceOfModelPart(part)
             for part in _Family(self.model, count=24, seed=2)]))

        # fitting anyway still works mechanically, which is the problem
        guard = geometry_guard_utils.CreateGeometryGuard(Kratos.Parameters("{}"))
        geometry_guard_utils.FitGeometryGuard(guard, family)
        self.assertIn(
            geometry_guard_utils.QueryGeometry(guard, family[0])["status"],
            ("OK", "WARN", "REJECT"))


@KratosUnittest.skipUnless(have_torch and have_guardrail,
                           "Missing required python modules: torch, physicsnemo (>= 2.2).")
class TestGeometryGuardCheck(KratosUnittest.TestCase):
    def setUp(self):
        self.model = Kratos.Model()
        self.sidecar = Path("test_geometry_guard_check.npz")
        self.family = _Family(self.model)
        guard = geometry_guard_utils.CreateGeometryGuard(Kratos.Parameters("{}"))
        geometry_guard_utils.FitGeometryGuard(guard, self.family)
        geometry_guard_utils.SaveGeometryGuard(guard, self.sidecar)

    def tearDown(self):
        kratos_utils.DeleteFileIfExisting(str(self.sidecar))

    def _Check(self, policy, model_part):
        check = geometry_guard_utils.GeometryGuardCheck(Kratos.Parameters("""{
            "guard_file" : "%s",
            "policy"     : "%s"
        }""" % (self.sidecar, policy)))
        return check, check.Check(model_part, "Test")

    def test_AdvisoryFlagsAndContinues(self):
        check, flagged = self._Check("advisory", _ScaledBox(self.model, "slab", (5.0, 1.0, 0.2)))
        self.assertTrue(flagged)
        self.assertEqual(check.last_status, "REJECT")

    def test_StrictRaises(self):
        with self.assertRaisesRegex(RuntimeError, "strict"):
            self._Check("strict", _ScaledBox(self.model, "slab2", (5.0, 1.0, 0.2)))

    def test_IgnoreAndEmptyFileAreDisabled(self):
        alien = _ScaledBox(self.model, "disabled_probe", (5.0, 1.0, 0.2))
        for settings in ('{ "guard_file" : "%s", "policy" : "ignore" }' % self.sidecar,
                         '{ "policy" : "advisory" }'):
            check = geometry_guard_utils.GeometryGuardCheck(Kratos.Parameters(settings))
            self.assertFalse(check.enabled)
            self.assertFalse(check.Check(alien, "Test"))

    def test_AnInDistributionShapeIsNotFlagged(self):
        check, flagged = self._Check("advisory", self.family[0])
        self.assertFalse(flagged)
        self.assertEqual(check.last_status, "OK")

    def test_TheVerdictIsCachedUntilTheMeshChanges(self):
        part = self.family[1]
        check = geometry_guard_utils.GeometryGuardCheck(Kratos.Parameters("""{
            "guard_file" : "%s"
        }""" % self.sidecar))
        self.assertFalse(check.Check(part, "Test"))
        # a second call on the same mesh must not rebuild the surface
        check._guard = None  # would raise on a reload attempt with no file set
        check.guard_file = "does_not_exist.npz"
        self.assertFalse(check.Check(part, "Test"))

    def test_UnknownPolicyRaises(self):
        with self.assertRaisesRegex(ValueError, "policy"):
            geometry_guard_utils.GeometryGuardCheck(
                Kratos.Parameters('{ "guard_file" : "x.npz", "policy" : "abort" }'))


if __name__ == '__main__':
    KratosUnittest.main()
