"""The shipped deformation layer against ShapeOptimizationApplication's
vertex morphing.

`ShapeOptimizationApplication` is not compiled in the reference environment,
and compiling it costs a full libKratosCore.so relink. Instead, the reference
fields in `shape_optimization_cases/vertex_morphing_reference.npz` were
produced once in a WHEEL-ONLY environment (KratosMultiphysics +
KratosShapeOptimizationApplication from PyPI, both GCC-built and therefore
internally consistent) by `shape_optimization_cases/generate_reference.py`.
Mixing a PyPI wheel with the locally compiled core would not work at all:
the wheels are GCC-built while this core is Clang-built, and pybind11 keys
its type registry on compiler identity.

**They are close relatives, not the same operator.** Both are compactly
supported, linear, kernel-weighted scattered-data maps, but they normalize
differently:

  Kratos MapperVertexMorphing   u = sum(w_j d_j) / sum(w_j)
  DeformPoints(..., "morph")    u = sum(a_j d_j) / (1 + sum(a_j))

Kratos is a true partition of unity, so a uniform control field maps to an
exact translation. physicsnemo's is a *regularized* (compact Shepard) field:
the `1 +` damps it, so it approaches a translation only as the weight sum
grows. With a control point at every node the two agree to floating-point
precision; with sparse control points the damping is plainly visible. The
kernels differ too - Kratos offers linear/gaussian/constant/cosine/quartic/
green, physicsnemo implements wendland_c2 only.

This file asserts what genuinely corresponds - support, linearity and the
dense-control limit - and pins what does not, so a change to either side
shows up rather than passing silently.
"""

from pathlib import Path

import numpy

import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.PhysicsNeMoApplication.bridges.mesh_bridge import deformation

try:
    import torch
    from physicsnemo.mesh import deformation as _pn_deformation  # noqa: F401
    have_deformation = True
except ImportError:
    have_deformation = False

_REFERENCE = Path(__file__).parent / "shape_optimization_cases" / "vertex_morphing_reference.npz"


@KratosUnittest.skipUnless(have_deformation,
                           "Missing required python modules: torch, physicsnemo >= 2.2.")
class TestVertexMorphingComparison(KratosUnittest.TestCase):

    def setUp(self):
        data = numpy.load(_REFERENCE)
        self.coordinates = data["coordinates"]
        self.node_ids = data["node_ids"]
        self.radius = float(data["radius"][0])
        self.centre_id = int(data["centre_node_id"][0])
        self.reference = {name: data[name] for name in data.files
                          if name.endswith(("_impulse", "_uniform"))}
        self.centre_row = int(numpy.where(self.node_ids == self.centre_id)[0][0])
        self.points = torch.as_tensor(self.coordinates, dtype=torch.float64)

    def _Morph(self, control_rows, control_points=None):
        control_points = self.coordinates if control_points is None else control_points
        deformed = deformation.DeformPoints(
            self.points,
            torch.as_tensor(control_rows, dtype=torch.float64),
            "morph",
            control_points=torch.as_tensor(control_points, dtype=torch.float64),
            radius=self.radius)
        return (deformed - self.points).numpy()

    # ---- what genuinely corresponds -------------------------------------

    def test_ReferenceItselfIsLocalAndTheShippedMorphMatchesThatSupport(self):
        # An impulse at one node moves only nodes within the filter radius.
        # Both operators must agree on WHICH nodes those are, even though
        # they disagree on the values.
        distance = numpy.linalg.norm(
            self.coordinates - self.coordinates[self.centre_row], axis=1)
        inside = distance < self.radius - 1e-12
        outside = distance > self.radius + 1e-12

        control = numpy.zeros_like(self.coordinates)
        control[self.centre_row, 2] = 1.0
        ours = self._Morph(control)

        for name, field in self.reference.items():
            if not name.endswith("_impulse"):
                continue
            with self.subTest(reference=name):
                self.assertEqual(numpy.abs(field[outside]).max(), 0.0)
                self.assertGreater(numpy.abs(field[inside]).max(), 0.0)
        self.assertEqual(numpy.abs(ours[outside]).max(), 0.0)
        self.assertGreater(numpy.abs(ours[inside]).max(), 0.0)

    def test_BothOperatorsAreLinear(self):
        first = numpy.zeros_like(self.coordinates)
        first[self.centre_row, 2] = 1.0
        second = numpy.zeros_like(self.coordinates)
        second[0, 0] = 0.5

        a = self._Morph(first)
        b = self._Morph(second)
        both = self._Morph(first + second)
        self.assertLess(numpy.abs(both - a - b).max(), 1e-12)

        # the reference's linearity is visible in the fixture: the uniform
        # case is the sum of impulses at every node, which for a normalized
        # operator equals the uniform translation
        for filter_name in ("linear", "gaussian"):
            with self.subTest(filter=filter_name):
                uniform = self.reference[f"{filter_name}_uniform"]
                self.assertLess(numpy.abs(uniform[:, 2] - 1.0).max(), 1e-9)

    # ---- what does NOT correspond, pinned ------------------------------

    def test_NormalizationDiffersButAgreesInTheDenseControlLimit(self):
        # Kratos divides by sum(w) - an exact partition of unity, so a
        # uniform control field is an exact translation.
        for filter_name in ("linear", "gaussian"):
            uniform = self.reference[f"{filter_name}_uniform"]
            self.assertLess(numpy.abs(uniform[:, 2] - 1.0).max(), 1e-9)

        # physicsnemo divides by (1 + sum(a)). With a control point at every
        # node the weight sum is large and the damping vanishes: the two
        # normalizations agree to floating-point precision.
        dense = numpy.zeros_like(self.coordinates)
        dense[:, 2] = 1.0
        self.assertLess(numpy.abs(self._Morph(dense)[:, 2] - 1.0).max(), 1e-9)

        # With SPARSE controls the damping is visible - the same uniform
        # field no longer maps to a translation. This is the case for which a
        # user must not assume vertex-morphing semantics.
        sparse = self._Morph(numpy.array([[0.0, 0.0, 1.0]]),
                             control_points=self.coordinates[[self.centre_row]])
        self.assertGreater(numpy.abs(sparse[:, 2] - 1.0).max(), 0.1)
        # ...but AT the control point it is exact: a coincident control has
        # unbounded weight, so a/(1+a) -> 1. Vertex morphing damps even there
        # (its impulse peak below is well under 1), which is the sharpest
        # statement of the difference: morph interpolates, morphing filters.
        self.assertAlmostEqual(sparse[self.centre_row, 2], 1.0, places=12)
        self.assertLess(self.reference["linear_impulse"][self.centre_row, 2], 0.5)

    def test_KernelsDifferSoTheImpulseResponsesDiffer(self):
        # Kratos offers linear/gaussian/constant/cosine/quartic/green;
        # physicsnemo's morph implements only wendland_c2. The two impulse
        # responses therefore have the same support but different shapes.
        control = numpy.zeros_like(self.coordinates)
        control[self.centre_row, 2] = 1.0
        ours = self._Morph(control)

        for filter_name in ("linear", "gaussian"):
            with self.subTest(filter=filter_name):
                field = self.reference[f"{filter_name}_impulse"]
                peak_ratio = ours[self.centre_row, 2] / field[self.centre_row, 2]
                self.assertGreater(abs(peak_ratio - 1.0), 1e-3)

    def test_RbfIsTheInterpolatingAnalogue(self):
        # "rbf" reproduces its control displacements AT the control points,
        # which vertex morphing does not; reporting which of the two shipped
        # methods is nearer keeps the guidance in the docs honest.
        # rbf carries a polynomial term and so needs at least D + 1 controls
        rows = [self.centre_row, 0, 1, 2]
        control_points = self.coordinates[rows]
        control = numpy.zeros((len(rows), 3))
        control[0, 2] = 1.0
        # polynomial=False is required here: the reference surface is planar,
        # so the polynomial term (1, x, y, z) is rank-deficient for any subset
        # of its nodes and the RBF solve is singular.
        deformed = deformation.DeformPoints(
            self.points, torch.as_tensor(control, dtype=torch.float64), "rbf",
            control_points=torch.as_tensor(control_points, dtype=torch.float64),
            polynomial=False)
        rbf = (deformed - self.points).numpy()
        # exact AT the control point - interpolation, not filtering
        self.assertAlmostEqual(rbf[self.centre_row, 2], 1.0, places=6)

        # vertex morphing is a filter: the peak is damped by the kernel sum
        field = self.reference["linear_impulse"]
        self.assertLess(field[self.centre_row, 2], 1.0)


if __name__ == '__main__':
    KratosUnittest.main()
