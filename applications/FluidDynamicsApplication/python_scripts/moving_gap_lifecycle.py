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

"""Reusable lifecycle for rigidly translating two-dimensional fluid patches.

The caller owns initial fluid values, boundary conditions and the linear
strategy. Clear any strategy before PrepareSolutionStep: rebuilding removes
its computing entities. Call FinalizeSolutionStep only after a successful solve.
"""
import math

import KratosMultiphysics as KM
from KratosMultiphysics.FluidDynamicsApplication.rebuild_fluid_patch_gap_process import RebuildFluidPatchGapProcess
from KratosMultiphysics.FluidDynamicsApplication.patch_history_transfer import (
    transfer_background_histories, OuterHistoryValidity, GapHistory)
from KratosMultiphysics.FluidDynamicsApplication.history_diagnostics import HistoryDiagnostics


class MovingGapLifecycle:
    """Coordinate BDF mesh motion, gap rebuilding and numerical history transfer.

    The model must already contain the initial computing part and coupling.
    Pass one settings/velocity pair or lists of pairs. This lifecycle supports
    rigidly translating inner patches, a fixed
    background, fixed time step and a three-step buffer. The caller must set
    initial nodal fluid values and the DELTA_TIME before creating the object.
    """

    def __init__(self, model, coupling_settings, velocity, dt, *,
                 tolerance=1e-10, bdf1_startup=True, history_diagnostics=None,
                 geometry_name='Interfaces', computing_name='Computing', motions=None, cylinder=None):
        settings = coupling_settings if isinstance(coupling_settings, (list, tuple)) else [coupling_settings]
        velocities = list(velocity) if velocity and isinstance(velocity[0], (list, tuple)) else [velocity]
        if len(settings) != len(velocities) or not all(
            len(v) == 2 and all(math.isfinite(x) for x in v) for v in velocities
        ):
            raise ValueError('Each inner patch requires two finite translation velocity components.')
        if not math.isfinite(dt) or dt <= 0:
            raise ValueError('time step must be finite and positive')
        self.velocities = [tuple(v) for v in velocities]
        self.motions = list(motions) if motions is not None else [None] * len(settings)
        if len(self.motions) != len(settings):
            raise ValueError('Each patch requires a motion specification.')
        for motion in self.motions:
            if motion is not None:
                if motion['type'] != 'sinusoidal' or not all(
                        math.isfinite(motion.get(key, 0.0)) for key in ('amplitude', 'frequency', 'phase')) or motion['frequency'] < 0.0:
                    raise ValueError('Sinusoidal motion requires finite amplitude and non-negative frequency.')
        self.dt = dt
        self.bdf1_startup = bdf1_startup
        self.rebuild = RebuildFluidPatchGapProcess(
            model, settings, geometry_name=geometry_name,
            computing_name=computing_name, tolerance=tolerance)
        self.root = self.rebuild.root
        self.background = self.rebuild.background
        self.inner_patches = self.rebuild.inner
        self.cylinder = cylinder
        self.computing = self.root.GetSubModelPart(computing_name)
        self._initial_time = self.root.ProcessInfo[KM.TIME]
        self._last_time = self._initial_time
        # Preserve the initial position and displacement, including restart offsets.
        self._initial_motion = {
            node.Id: (node.X, node.Y, tuple(node.GetSolutionStepValue(KM.DISPLACEMENT)))
            for patch in self.inner_patches for node in patch.Nodes
        }
        self.inner_ids = {node.Id for patch in self.inner_patches for node in patch.Nodes}
        previous = {node.Id for node in self.computing.Nodes} - self.inner_ids
        self.validity = OuterHistoryValidity(previous, () if bdf1_startup else previous)
        self.gap_histories = [GapHistory() for patch in self.inner_patches]
        for i, history in enumerate(self.gap_histories):
            history.capture(self._get_gap_geometry(i), self._initial_time)
        self.bdf = KM.TimeDiscretization.BDF(2)
        self.previous_nodes = {node.Id for node in self.computing.Nodes}
        self.previous_elements = {e.Id for e in self.background.Elements if e.Is(KM.ACTIVE)}
        self.previous_triangles = sum(self._get_gap_geometry(i).GetSubModelPart('GapTriangles').NumberOfGeometries()
                                      for i in range(len(self.inner_patches)))
        self.step = 0
        self._pending = False
        self._failed = False
        self.diagnostics = HistoryDiagnostics(history_diagnostics) if history_diagnostics is not None else None
        if self.diagnostics:
            self._geometry_slots = [self._snapshot_geometry(), None, None]
            self._activity_slots = [set(previous), None, None]
            # Initial values may be valid supplied data; they were not solved
            # by this lifecycle. Older computation provenance is unavailable.
            self._computed_slots = [set(), None, None]

    def _get_gap_geometry(self, index):
        return self.rebuild.GetPatchGeometry(index)

    def _translation(self, index, time):
        motion = self.motions[index]
        if motion is None:
            elapsed = time - self._initial_time
            return tuple(v * elapsed for v in self.velocities[index])
        omega = 2.0 * math.pi * motion['frequency']
        phase = motion.get('phase', 0.0)
        return (motion['amplitude'] * (math.sin(omega * time + phase) - math.sin(omega * self._initial_time + phase)), 0.0)

    def _snapshot_geometry(self):
        """Copy geometry before the next rebuild deletes its gap entities."""
        def vertices(geometry):
            return tuple((float(node.X), float(node.Y)) for node in geometry)
        interfaces = [self._get_gap_geometry(i) for i in range(len(self.inner_patches))]
        return {
            'inner': [(e.Id, vertices(e.GetGeometry())) for patch in self.inner_patches for e in patch.Elements],
            'active_outer': [(e.Id, vertices(e.GetGeometry())) for e in self.background.Elements if e.Is(KM.ACTIVE)],
            'gap': [(g.Id, vertices(g)) for interface in interfaces for g in interface.GetSubModelPart('GapTriangles').Geometries]
        }

    def prepare_solution_step(self, time, *, exact_history=None):
        """Move and rebuild; optionally supply exact_history(node, historical_time).

        If exact_history is None, numerical transfer uses Kratos history sources.
        Returns computing part and motion/transfer diagnostics. Callers must set
        forcing and boundary data, solve, then finalize the step.
        """
        if self._pending:
            raise RuntimeError('Finalize the previous successful step first.')
        if self._failed:
            raise RuntimeError('The previous preparation failed after changing the model; restore a checkpoint before retrying.')
        if not math.isfinite(time) or not math.isclose(time-self._last_time, self.dt, rel_tol=1e-10, abs_tol=1e-12):
            raise ValueError('The moving gap lifecycle requires successive fixed time steps.')
        # Rebuilding mutates Kratos containers and cannot be rolled back here.
        # Keep the lifecycle unusable after an exception until the model is restored.
        self._failed = True
        self.step += 1
        self.root.CloneTimeStep(time)
        self.validity.advance()
        for history in self.gap_histories:
            history.advance()
        if self.diagnostics:
            self._geometry_slots = [None, self._geometry_slots[0], self._geometry_slots[1]]
            self._activity_slots = [None, self._activity_slots[0], self._activity_slots[1]]
            self._computed_slots = [set(), self._computed_slots[0], self._computed_slots[1]]
        self.bdf.ComputeAndSaveBDFCoefficients(self.root.ProcessInfo)
        if self.bdf1_startup and self.step == 1:
            self.root.ProcessInfo[KM.BDF_COEFFICIENTS] = KM.Vector([1/self.dt, -1/self.dt, 0.])
        coefficients = self.root.ProcessInfo[KM.BDF_COEFFICIENTS]
        for index, patch in enumerate(self.inner_patches):
            translation = self._translation(index, time)
            for node in patch.Nodes:
                x, y, displacement = self._initial_motion[node.Id]
                node.X = x + translation[0]
                node.Y = y + translation[1]
                node.SetSolutionStepValue(KM.DISPLACEMENT, [displacement[0]+translation[0],
                    displacement[1]+translation[1], displacement[2]])
                mesh_velocity = [sum(coefficients[k]*node.GetSolutionStepValue(KM.DISPLACEMENT, k)[d]
                                     for k in range(3)) for d in range(2)]
                node.SetSolutionStepValue(KM.MESH_VELOCITY, [*mesh_velocity, 0.])
            if self.cylinder is not None and self.cylinder.patch.FullName() == patch.FullName():
                motion = self.motions[index]
                if motion is None:
                    body_velocity = self.velocities[index]
                else:
                    omega = 2.0 * math.pi * motion['frequency']
                    body_velocity = (motion['amplitude'] * omega * math.cos(omega*time + motion.get('phase', 0.0)), 0.0)
                self.cylinder.update(translation, body_velocity)
        self.rebuild.Execute()
        self.computing = self.rebuild.computing
        current_nodes = {node.Id for node in self.computing.Nodes}
        current_elements = {e.Id for e in self.background.Elements if e.Is(KM.ACTIVE)}
        required = (1,) if coefficients[2] == 0.0 else (1, 2)
        targets = [self.root.GetNode(i) for i in sorted(current_nodes-self.inner_ids)
                   if any(i not in self.validity.slots[k] for k in required)]
        counts = {'outer': 0, 'inner': 0, 'gap': 0, 'extrapolation': 0, 'missing': 0}
        samples = [] if self.diagnostics else None
        if self.diagnostics:
            self._activity_slots[0] = set(current_nodes-self.inner_ids)
            self._geometry_slots[0] = self._snapshot_geometry()
            new_outer = current_nodes-self.previous_nodes-self.inner_ids
            observe = self.diagnostics.make_observer(
                self.step, time, self.dt, new_outer, self.validity,
                self._geometry_slots, self._activity_slots, self._computed_slots, samples)
            events = {'activated_nodes': len(new_outer),
                      'deactivated_nodes': len((self.previous_nodes-self.inner_ids)-current_nodes),
                      'activated_elements': len(current_elements-self.previous_elements),
                      'deactivated_elements': len(self.previous_elements-current_elements)}
        if exact_history is not None:
            for node in targets:
                for k in required:
                    if node.Id not in self.validity.slots[k]:
                        node.SetSolutionStepValue(KM.VELOCITY, k, exact_history(node, time-k*self.dt))
                        if self.diagnostics:
                            observe(node, k, 'exact', None, None)
                    elif self.diagnostics:
                        observe(node, k, 'outer', None, None)
        elif targets:
            try:
                offsets = []
                for i in range(len(self.inner_patches)):
                    current = self._translation(i, time)
                    offsets.append({k: tuple(current[d] - self._translation(i, time-k*self.dt)[d]
                                             for d in range(2)) for k in required})
                counts = transfer_background_histories(
                    self.inner_patches[0], targets, time, self.dt, self.velocities[0], self.validity,
                    required, observer=observe if self.diagnostics else None,
                    inner_patches=self.inner_patches, velocities=self.velocities, gap_histories=self.gap_histories,
                    displacement_offsets=offsets)
            except Exception:
                if self.diagnostics:
                    self.diagnostics.write_step(self.step, time, samples, events, committed=False)
                raise
        if self.diagnostics:
            # Nodes already holding every required valid slot bypass transfer;
            # record their reuse without modifying the numerical target set.
            target_ids = {node.Id for node in targets}
            for node_id in sorted(new_outer-target_ids):
                node = self.root.GetNode(node_id)
                for slot in required:
                    observe(node, slot, 'outer', None, None)
            self.diagnostics.write_step(self.step, time, samples, events)
        triangles = [g for i in range(len(self.inner_patches))
                     for g in self._get_gap_geometry(i).GetSubModelPart('GapTriangles').Geometries]
        result = {'newly_active_nodes': len(current_nodes-self.previous_nodes-self.inner_ids),
                  'newly_active_elements': len(current_elements-self.previous_elements),
                  'deactivated_elements': len(self.previous_elements-current_elements),
                  'gap_triangles': len(triangles),
                  'gap_triangle_change': len(triangles)-self.previous_triangles,
                  'gap_area': sum(g.Area() for g in triangles),
                  'missing_history_samples': counts['missing'],
                  'transferred_history_samples': counts['inner'],
                  'interpolated_history_samples': counts['inner'],
                  'outer_history_samples': counts['outer'], 'gap_history_samples': counts['gap'],
                  'extrapolated_history_samples': counts['extrapolation']}
        self._current_nodes = current_nodes
        self._current_elements = current_elements
        self._current_triangles = len(triangles)
        self._pending = True
        self._failed = False
        return self.computing, result

    def finalize_solution_step(self, time):
        """Mark solved background values and capture gap support for future BDF steps."""
        if not self._pending:
            raise RuntimeError('Prepare and solve a step before finalizing it.')
        if not math.isclose(time, self.root.ProcessInfo[KM.TIME], rel_tol=1e-12, abs_tol=1e-12):
            raise ValueError('Finalization time must match the prepared time step.')
        for i, history in enumerate(self.gap_histories):
            history.capture(self._get_gap_geometry(i), time)
        self.validity.mark_solved(self._current_nodes-self.inner_ids)
        if self.diagnostics:
            self._computed_slots[0] = set(self._current_nodes-self.inner_ids)
        self.previous_nodes = self._current_nodes
        self.previous_elements = self._current_elements
        self.previous_triangles = self._current_triangles
        self._last_time = time
        self._pending = False
