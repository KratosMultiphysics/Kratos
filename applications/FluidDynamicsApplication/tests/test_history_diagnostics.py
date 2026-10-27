"""Focused tests of history classification and diagnostic output."""
import csv
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.FluidDynamicsApplication.history_diagnostics import (
    HistoryDiagnostics, classify_geometry)
from KratosMultiphysics.FluidDynamicsApplication import patch_history_transfer


class TestHistoryDiagnostics(KratosUnittest.TestCase):
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
                object(), [node], 1., .1, (0., 0.), Validity(), (1, 2), Gap(),
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
