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

"""Focused tests of history classification and diagnostic output."""
import csv
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
from unittest.mock import patch

import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.FluidDynamicsApplication.history_diagnostics import (
    HistoryDiagnostics, classify_geometry)
from KratosMultiphysics.FluidDynamicsApplication import patch_history_transfer


class TestHistoryDiagnostics(KratosUnittest.TestCase):
    def test_multiple_patch_history_uses_each_translation(self):
        class Node:
            Id, X, Y = 7, 0., 0.
            def __init__(self):
                self.values = {}
            def SetSolutionStepValue(self, variable, *args):
                self.values[args[0] if len(args) == 2 else 0] = args[-1]
            def GetSolutionStepValue(self, variable, slot):
                return self.values[slot]
        class Source:
            def GetSolutionStepValue(self, variable, slot):
                return [2., 3., 0.]
        class Element:
            Id = 42
            def GetGeometry(self):
                return [Source() for _ in range(4)]
        class Locator:
            def __init__(self, patch):
                self.patch = patch
            def UpdateSearchDatabase(self):
                pass
            def FindPointOnMesh(self, point):
                calls.append(self.patch.name)
                if self.patch.name == 2 and abs(point[1] - .2) < 1e-12:
                    return True, [.25]*4, Element()
                return False, None, None
        node = Node()
        calls = []
        first = SimpleNamespace(name=1, Elements=[], Nodes=[SimpleNamespace(X=10., Y=10.)])
        # A rotated square's bounding box includes its extrema, including the query.
        second = SimpleNamespace(name=2, Elements=[], Nodes=[SimpleNamespace(X=x, Y=y) for x, y in
                                               ((0., .2), (.1, .1), (0., 0.), (-.1, .1))])
        with patch.object(patch_history_transfer.KM, 'BinBasedFastPointLocator2D', Locator):
            counts = patch_history_transfer.transfer_background_histories(
                first, [node], 1., .1, (1., 0.), required_steps=(1,),
                inner_patches=[first, second], velocities=[(1., 0.), (0., 2.)],
                gap_histories=[None, None])
        self.assertEqual(counts['inner'], 1)
        self.assertEqual(calls, [2])
        self.assertEqual(node.values[1], [2., 3., 0.])

    def test_historical_displacement_overrides_constant_velocity(self):
        class Node:
            Id, X, Y = 7, 0., 0.
            def __init__(self):
                self.values = {}
            def SetSolutionStepValue(self, variable, *args):
                self.values[args[0] if len(args) == 2 else 0] = args[-1]
            def GetSolutionStepValue(self, variable, slot):
                return self.values[slot]
        class Source:
            def GetSolutionStepValue(self, variable, slot):
                return [1., 2., 0.]
        class Element:
            Id = 5
            def GetGeometry(self):
                return [Source()] * 4
        class Locator:
            def __init__(self, part): pass
            def UpdateSearchDatabase(self): pass
            def FindPointOnMesh(self, point):
                self_point.append(point)
                return True, [.25]*4, Element()
        self_point = []
        part = SimpleNamespace(Elements=[], Nodes=[SimpleNamespace(X=0., Y=0.),
                                                   SimpleNamespace(X=1., Y=1.)])
        with patch.object(patch_history_transfer.KM, 'BinBasedFastPointLocator2D', Locator):
            patch_history_transfer.transfer_background_histories(
                part, [Node()], 1., .1, (10., 0.), required_steps=(1,),
                displacement_offsets=[{1: (.3, .4)}])
        self.assertEqual(self_point, [[.3, .4, 0.]])

    def test_step_observer_records_historical_location_and_validity(self):
        class Node:
            Id, X, Y = 7, 0.5, 0.5
        class Validity:
            slots = [set(), {7}, set()]
        square = ((0., 0.), (1., 0.), (1., 1.), (0., 1.))
        geometry = {'inner': [(3, square)], 'active_outer': [], 'gap': []}
        with TemporaryDirectory() as directory:
            diagnostics = HistoryDiagnostics(directory)
            samples = []
            observe = diagnostics.make_observer(
                2, 1.2, .1, {7}, Validity(), [None, geometry, None],
                [set(), {7}, None], [set(), {7}, None], samples)
            observe(Node(), 1, 'outer', None, None)
            self.assertEqual(samples[0]['location'], 'inner')
            self.assertAlmostEqual(samples[0]['historical_time'], 1.1)
            self.assertTrue(samples[0]['valid_previous'])
            self.assertTrue(samples[0]['computed_previous'])
            self.assertTrue(samples[0]['suspicious'])

    def test_transfer_observer_reports_outer_and_gap_without_changing_values(self):
        class Node:
            Id, X, Y = 7, .25, .25
            def __init__(self):
                self.values = {1: [9., 8., 0.], 2: [0., 0., 0.]}
            def SetSolutionStepValue(self, variable, *args):
                if len(args) == 2:
                    self.values[args[0]] = args[1]
                else:
                    self.values[0] = args[0]
            def GetSolutionStepValue(self, variable, slot):
                return self.values[slot]
        class Locator:
            def __init__(self, inner): pass
            def UpdateSearchDatabase(self): pass
            def FindPointOnMesh(self, point): return False, None, None
        class Gap:
            def interpolate(self, slot, time, x, y): return [1., 2., 0.]
        class Validity:
            slots = [set(), {7}, set()]
        node = Node()
        observations = []
        with patch.object(patch_history_transfer.KM, 'BinBasedFastPointLocator2D', Locator):
            counts = patch_history_transfer.transfer_background_histories(
                SimpleNamespace(Elements=[], Nodes=[SimpleNamespace(X=.25, Y=.25)]),
                [node], 1., .1, (0., 0.), Validity(), (1, 2), Gap(),
                observer=lambda node, slot, source, parent, shape: observations.append((slot, source)))
        self.assertEqual(observations, [(1, 'outer'), (2, 'gap')])
        self.assertEqual(counts['outer'], 1)
        self.assertEqual(counts['gap'], 1)
        self.assertEqual(node.values[1], [9., 8., 0.])
        self.assertEqual(node.values[2], [1., 2., 0.])

    def test_location_is_independent_and_unknown_is_explicit(self):
        square = ((0., 0.), (1., 0.), (1., 1.), (0., 1.))
        snapshot = {'inner': [(7, square)], 'active_outer': [], 'gap': []}
        self.assertEqual(classify_geometry(None, .5, .5), ('unknown', ''))
        self.assertEqual(classify_geometry(snapshot, .5, .5), ('inner', '7'))
        self.assertEqual(classify_geometry(snapshot, 2., 2.), ('uncovered', ''))

    def test_shared_boundary_is_not_arbitrarily_classified(self):
        left = ((0., 0.), (1., 0.), (1., 1.), (0., 1.))
        right = ((1., 0.), (2., 0.), (2., 1.), (1., 1.))
        snapshot = {'inner': [(1, left)], 'active_outer': [(2, right)], 'gap': []}
        self.assertEqual(classify_geometry(snapshot, 1., .5), ('boundary', ''))

    def test_summary_distinguishes_node_events_and_slot_sources(self):
        with TemporaryDirectory() as directory:
            output = HistoryDiagnostics(directory)
            sample = dict.fromkeys(output.columns, '')
            sample.update(step=1, node_id=5, history_slot=1, source='outer',
                          location='inner', newly_active=True, suspicious=True)
            output.write_step(1, .1, [sample], {'activated_nodes': 1, 'deactivated_nodes': 0,
                              'activated_elements': 3, 'deactivated_elements': 2})
            summary = json.loads((Path(directory)/'history_summary.json').read_text())
            self.assertEqual(summary['cumulative']['activated_nodes'], 1)
            self.assertEqual(summary['cumulative']['activated_elements'], 3)
            self.assertEqual(summary['cumulative']['sources']['outer'], 1)
            self.assertEqual(summary['cumulative']['locations']['inner'], 1)
            self.assertEqual(summary['cumulative']['suspicious'], 1)
            with (Path(directory)/'history_samples.csv').open(newline='') as stream:
                self.assertEqual(len(list(csv.DictReader(stream))), 1)


if __name__ == '__main__':
    KratosUnittest.main()
