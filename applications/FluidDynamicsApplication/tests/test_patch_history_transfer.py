#     |  /           |
#     ' /   __| _` | __|  _ \   __|
#     . \  |   (   | |   (   |\__ `
#    _|\_\_|  \__,_|\__|\___/ ____/
#                    Multi-Physics
#
#   Fluid Dynamics Application
#
#   License:         BSD License
#                    Kratos default license: kratos/license.txt
#
#   Main author:     Nicolò Antonelli
#

"""History validity, inner interpolation and uncovered points for BDF startup."""
import math
import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.FluidDynamicsApplication.patch_history_transfer import (
    GapHistory, OuterHistoryValidity, transfer_background_histories)


class BoundedHistoryTransferTest(KratosUnittest.TestCase):
    def setUp(self):
        self.model = KM.Model()
        self.root = self.model.CreateModelPart('Fluid')
        self.root.AddNodalSolutionStepVariable(KM.VELOCITY)
        self.root.SetBufferSize(3)
        self.inner = self.root.CreateSubModelPart('Inner')
        self.speed, self.dt = (0.01, -0.02), 0.01
        self.c, self.s = math.cos(0.4), math.sin(0.4)
        for i, (x, y) in enumerate(((0, 0), (1, 0), (1, 1), (0, 1)), 1):
            node = self.inner.CreateNewNode(i, self.c*x-self.s*y, self.s*x+self.c*y, 0.0)
            for k in (1, 2):
                hx, hy = node.X-self.speed[0]*k*self.dt, node.Y-self.speed[1]*k*self.dt
                node.SetSolutionStepValue(KM.VELOCITY, k, [2*hx+3*hy+k, -hx+4*hy-k, 0.0])
        self.inner.CreateNewElement('Element2D4N', 1, [1, 2, 3, 4], self.root.CreateNewProperties(0))

    def target(self, node_id, x, y):
        return self.root.CreateNewNode(node_id, self.c*x-self.s*y, self.s*x+self.c*y, 0.0)

    def test_mixed_outer_and_inner_history(self):
        target = self.target(10, 0.5, 0.5)
        target.SetSolutionStepValue(KM.VELOCITY, 2, [9.0, 8.0, 0.0])
        validity = OuterHistoryValidity((target.Id,), (target.Id,))
        validity.advance()
        validity.slots[1].clear()
        counts = transfer_background_histories(self.inner, [target], 0.2, self.dt, self.speed, validity)
        self.assertEqual((counts['outer'], counts['inner']), (1, 1))
        self.assertEqual(target.GetSolutionStepValue(KM.VELOCITY, 2)[0], 9.0)
        self.assertAlmostEqual(target.GetSolutionStepValue(KM.VELOCITY, 1)[0], 2*target.X+3*target.Y+1, places=10)

    def test_uncovered_point_fails_without_partial_writes(self):
        inside, outside = self.target(10, 0.5, 0.5), self.target(11, 1.01, 0.5)
        inside.SetSolutionStepValue(KM.VELOCITY, 1, [9.0, 8.0, 0.0])
        with self.assertRaisesRegex(RuntimeError, 'extrapolation is disabled'):
            transfer_background_histories(self.inner, [inside, outside], 0.2, self.dt, self.speed)
        self.assertEqual(inside.GetSolutionStepValue(KM.VELOCITY, 1)[0], 9.0)

    def test_bdf1_startup_requires_only_first_slot(self):
        target = self.target(10, 0.5, 0.5)
        counts = transfer_background_histories(self.inner, [target], 0.01, self.dt, self.speed,
                                               required_steps=(1,))
        self.assertEqual(counts['inner'], 1)

    def test_gap_parent_q1_reconstruction_and_time(self):
        # The target lies outside the rotated inner mesh. The gap triangle
        # carries a Cartesian parent's full Q1 polynomial, including xy.
        target = self.root.CreateNewNode(10, 2.25, 0.25, 0.0)
        parent = tuple((i, x, y, (x*y + 2*x, y, 0.0)) for i, (x, y) in
                       enumerate(((2., 0.), (3., 0.), (3., 1.), (2., 1.)), 20))
        cells = [(((2., 0.), (3., 0.), (2., 1.)), parent)]
        history = GapHistory()
        history.slots[1] = (0.19, cells)
        counts = transfer_background_histories(self.inner, [target], 0.2, self.dt, self.speed,
                                                required_steps=(1,), gap_history=history)
        self.assertEqual(counts['gap'], 1)
        self.assertAlmostEqual(target.GetSolutionStepValue(KM.VELOCITY, 1)[0],
                               target.X*target.Y+2*target.X, places=12)
        with self.assertRaisesRegex(RuntimeError, 'time mismatch'):
            history.interpolate(1, 0.18, target.X, target.Y)

    def test_gap_interface_with_distinct_parent_extensions_is_rejected(self):
        parent_a = tuple((i, x, y, (0., 0., 0.)) for i, (x, y) in
                         enumerate(((0., 0.), (1., 0.), (1., 1.), (0., 1.)), 20))
        parent_b = tuple((i, x, y, (1., 0., 0.)) for i, (x, y) in
                         enumerate(((0., 0.), (1., 0.), (1., 1.), (0., 1.)), 30))
        history = GapHistory()
        history.slots[1] = (0.0, [(((0., 0.), (1., 0.), (0., 1.)), parent_a),
                                   (((1., 1.), (0., 1.), (1., 0.)), parent_b)])
        with self.assertRaisesRegex(RuntimeError, 'Ambiguous historical gap interface'):
            history.interpolate(1, 0.0, 0.5, 0.5)

    def test_snapshot_uses_gap_geometry_and_freezes_parent_values(self):
        from test_fluid_patch_gap_coupling import create_gap_problem
        _, modeler, _, _ = create_gap_problem()
        root = modeler.model_part
        for node in modeler.patches[0].Nodes:
            node.SetSolutionStepValue(KM.VELOCITY, [node.X*node.Y, node.X+node.Y, 0.0])
        history = GapHistory()
        history.capture(root.GetSubModelPart('Interfaces.Patch_1'), 0.0)
        triangle, _ = history.slots[0][1][0]
        x = sum(p[0] for p in triangle)/3
        y = sum(p[1] for p in triangle)/3
        value = history.interpolate(0, 0.0, x, y)
        self.assertAlmostEqual(value[0], x*y, places=10)
        for node in modeler.patches[0].Nodes:
            node.SetSolutionStepValue(KM.VELOCITY, [0.0, 0.0, 0.0])
        self.assertAlmostEqual(history.interpolate(0, 0.0, x, y)[0], x*y, places=10)


if __name__ == '__main__':
    KratosUnittest.main()
