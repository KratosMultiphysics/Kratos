#     |  /           |
#     ' /   __| _` | __|  _ \   __|
#     . \  |   (   | |   (   |\__ `
#    _|\_\_|  \__,_|\__|\___/ ____/
#                    Multi-Physics
#
#   License: BSD License (kratos/license.txt)
#
"""Same-phase Q1 history snapshots for fixed Cartesian two-fluid patches."""
import math

import KratosMultiphysics as KM


def _sample_parent(parent, x, y):
    a, b, _, d = parent
    ex, ey = b[1] - a[1], b[2] - a[2]
    fx, fy = d[1] - a[1], d[2] - a[2]
    determinant = ex * fy - ey * fx
    if abs(determinant) <= 1e-30:
        raise RuntimeError("Degenerate phase-history parent")
    s = ((x - a[1]) * fy - (y - a[2]) * fx) / determinant
    t = (ex * (y - a[2]) - ey * (x - a[1])) / determinant
    shape = ((1 - s) * (1 - t), s * (1 - t), s * t, (1 - s) * t)
    return [sum(shape[i] * parent[i][3][d] for i in range(4)) for d in range(3)]


def _distance_to_polygon(polygon, x, y, tolerance):
    signs, distance = [], math.inf
    for a, b in zip(polygon, polygon[1:] + polygon[:1]):
        dx, dy = b[0] - a[0], b[1] - a[1]
        length2 = dx * dx + dy * dy
        if length2 <= 0.0:
            raise RuntimeError("Degenerate historical phase polygon")
        signs.append(dx * (y - a[1]) - dy * (x - a[0]))
        fraction = max(0.0, min(1.0, ((x - a[0]) * dx + (y - a[1]) * dy) / length2))
        distance = min(distance, math.hypot(x - a[0] - fraction * dx, y - a[1] - fraction * dy))
    if all(sign >= -tolerance for sign in signs) or all(sign <= tolerance for sign in signs):
        return 0.0
    return distance


class PhaseHistory:
    """Store physical bulk/gap support and Q1 velocities at prior times."""

    def __init__(self, phase_name):
        self.phase_name = phase_name
        self.slots = [None, None, None]

    def advance(self):
        self.slots = [None, self.slots[0], self.slots[1]]

    def capture(self, phase, gap_geometry, time):
        if not math.isfinite(time):
            raise ValueError("History time must be finite")

        def parent_snapshot(nodes):
            if len(nodes) != 4:
                raise ValueError("Phase history requires Q1 parent elements")
            return tuple((node.Id, node.X, node.Y, tuple(node.GetSolutionStepValue(KM.VELOCITY))) for node in nodes)

        cells = []
        for element in phase.Elements:
            if element.IsDefined(KM.ACTIVE) and element.IsNot(KM.ACTIVE):
                continue
            nodes = list(element.GetGeometry())
            cells.append((tuple((node.X, node.Y) for node in nodes), parent_snapshot(nodes), "bulk"))
        volumes = {geometry.Id: geometry for geometry in gap_geometry.GetSubModelPart("GapVolumes").Geometries}
        triangles = list(gap_geometry.GetSubModelPart("GapTriangles").Geometries)
        if len(volumes) != len(triangles):
            raise RuntimeError("Gap triangle and quadrature counts differ")
        for triangle in triangles:
            if triangle.Id - 1 not in volumes:
                raise RuntimeError("Gap triangle has no parent quadrature geometry")
            parent = parent_snapshot(list(volumes[triangle.Id - 1]))
            cells.append((tuple((node.X, node.Y) for node in triangle), parent, "gap"))
        self.slots[0] = (float(time), cells)

    def sample(self, slot, time, x, y, tolerance=1e-10):
        if slot not in (1, 2) or self.slots[slot] is None:
            raise RuntimeError(f"Missing {self.phase_name} history in slot {slot}")
        recorded_time, cells = self.slots[slot]
        if not math.isclose(recorded_time, time, rel_tol=1e-12, abs_tol=1e-12):
            raise RuntimeError(f"{self.phase_name} history time mismatch: {time} versus {recorded_time}")
        candidates = []
        for polygon, parent, source in cells:
            distance = _distance_to_polygon(polygon, x, y, tolerance * tolerance)
            h = min(math.hypot(parent[(i + 1) % 4][1] - parent[i][1],
                               parent[(i + 1) % 4][2] - parent[i][2]) for i in range(4))
            if distance <= h + tolerance:
                candidates.append((distance, tuple(node[0] for node in parent), source,
                                   _sample_parent(parent, x, y)))
        if not candidates:
            raise RuntimeError(f"No same-phase history for {self.phase_name}, slot {slot}, "
                               f"time {time}, point ({x}, {y}) within one cell width")
        candidates.sort(key=lambda candidate: (candidate[0], candidate[1]))
        nearest = candidates[0]
        for candidate in candidates[1:]:
            if candidate[0] > nearest[0] + tolerance:
                break
            if max(abs(a - b) for a, b in zip(candidate[3], nearest[3])) > tolerance:
                raise RuntimeError(f"Ambiguous {self.phase_name} history at ({x}, {y}), time {time}")
        return nearest[3], nearest[2] if nearest[0] <= tolerance else "extrapolation"


def make_current_phase_velocity_sampler(water_history, air_history, tolerance=1e-10):
    """Sample captured physical phase support for level-set convection.

    At an internal seam the extension has two traces; use their arithmetic
    mean, counting each parent only once. The transport process separately
    averages the water and air values on the physical interface.
    """
    if not math.isfinite(tolerance) or tolerance <= 0.0:
        raise ValueError("tolerance must be finite and positive")
    histories = (water_history, air_history)
    if any(history.slots[0] is None for history in histories):
        raise ValueError("Capture both phase fields before constructing a velocity sampler")
    if not math.isclose(water_history.slots[0][0], air_history.slots[0][0], rel_tol=1e-12, abs_tol=1e-12):
        raise ValueError("Water and air velocity snapshots must have the same time")
    snapshots = tuple(history.slots[0][1] for history in histories)

    def sample(x, y, phase):
        if phase not in (0, 1) or not math.isfinite(x) or not math.isfinite(y):
            raise ValueError("Expected a finite point and phase 0 (water) or 1 (air)")
        parents = {}
        for polygon, parent, _ in snapshots[phase]:
            if _distance_to_polygon(polygon, x, y, tolerance * tolerance) <= tolerance:
                key = tuple(node[0] for node in parent)
                parents[key] = parent
        if not parents:
            raise RuntimeError(f"No physical {histories[phase].phase_name} velocity support at ({x}, {y})")
        values = [_sample_parent(parent, x, y) for parent in parents.values()]
        return [sum(value[d] for value in values) / len(values) for d in range(3)]

    return sample


def transfer_physical_velocity_histories(histories, phase_nodes, time, dt, required_slots, tolerance=1e-10):
    """Transfer the previous physical velocity at each current phase node.

    At a material interface both physical traces are averaged. No phase
    extrapolation is used: the source must be in the previous physical tank.
    All samples are prepared before modifying any nodal history.
    """
    if (len(histories) != 2 or len(phase_nodes) != 2 or not math.isfinite(time)
            or not math.isfinite(dt) or dt <= 0.0 or not math.isfinite(tolerance) or tolerance <= 0.0):
        raise ValueError("Expected two phase histories, finite time and positive dt and tolerance")
    for slot in required_slots:
        if slot not in (1, 2):
            raise ValueError("History slot must be 1 or 2")
        for history in histories:
            snapshot = history.slots[slot]
            if snapshot is None or not math.isclose(snapshot[0], time - slot * dt, rel_tol=1e-12, abs_tol=1e-12):
                raise RuntimeError(f"Missing {history.phase_name} physical velocity at time {time - slot * dt}")

    values = []
    counts = {"same_phase": 0, "other_phase": 0, "interface": 0}
    for target_phase, nodes in enumerate(phase_nodes):
        for node in nodes:
            for slot in required_slots:
                phase_values = []
                for history in histories:
                    parents = {}
                    for polygon, parent, _ in history.slots[slot][1]:
                        if _distance_to_polygon(polygon, node.X, node.Y, tolerance * tolerance) <= tolerance:
                            parents[tuple(item[0] for item in parent)] = parent
                    if parents:
                        samples = [_sample_parent(parent, node.X, node.Y) for parent in parents.values()]
                        phase_values.append([sum(sample[d] for sample in samples) / len(samples)
                                             for d in range(3)])
                    else:
                        phase_values.append(None)
                available = [value for value in phase_values if value is not None]
                if not available:
                    raise RuntimeError(f"No previous physical velocity for node {node.Id}, phase "
                                       f"{histories[target_phase].phase_name}, time {time - slot * dt}, "
                                       f"point ({node.X}, {node.Y})")
                source = "interface" if len(available) == 2 else (
                    "same_phase" if phase_values[target_phase] is not None else "other_phase")
                counts[source] += 1
                values.append((node, slot, [sum(value[d] for value in available) / len(available)
                                             for d in range(3)]))
    for node, slot, velocity in values:
        node.SetSolutionStepValue(KM.VELOCITY, slot, velocity)
    return counts


def transfer_phase_histories(history, nodes, time, dt, required_slots, valid_slots):
    """Fill missing same-phase velocity slots atomically; return source counts."""
    if not math.isfinite(time) or not math.isfinite(dt) or dt <= 0.0:
        raise ValueError("History transfer requires finite time and positive dt")
    values, counts = {}, {"existing": 0, "bulk": 0, "gap": 0, "extrapolation": 0}
    for node in nodes:
        for slot in required_slots:
            if slot not in (1, 2):
                raise ValueError("History slot must be 1 or 2")
            if node.Id in valid_slots[slot]:
                counts["existing"] += 1
                continue
            try:
                velocity, source = history.sample(slot, time - slot * dt, node.X, node.Y)
            except RuntimeError as exc:
                raise RuntimeError(f"Cannot transfer node {node.Id}, phase {history.phase_name}, slot {slot}: {exc}") from exc
            values[node.Id, slot] = velocity
            counts[source] += 1
    for node in nodes:
        for slot in required_slots:
            if (node.Id, slot) in values:
                node.SetSolutionStepValue(KM.VELOCITY, slot, values[node.Id, slot])
    return counts
