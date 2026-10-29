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

"""BDF2 velocity history transfer for rigidly translating Q1 patches."""
import math
import KratosMultiphysics as KM


# Disabled Q1 extrapolation retained for future work.
# def _disabled_q1_extrapolation_reference(inner, nodes, time, dt, speed):
#     """Interpolate or bound Q1 extrapolation at the two historical positions.
#
#     Rigid translation maps each historical query to the current mesh. Nearest
#     boundary edges determine distance and the Q1 parent for extrapolation.
#     The permitted distance is 10% of that parent's shortest edge length.
#     """
#     if isinstance(speed, (int, float)):
#         speed = (speed, 0.0)
#     locator = KM.BinBasedFastPointLocator2D(inner)
#     locator.UpdateSearchDatabase()
#     edges = {}
#     for element in inner.Elements:
#         geometry = list(element.GetGeometry())
#         for i in range(4):
#             a, b = geometry[i], geometry[(i+1) % 4]
#             key = tuple(sorted((a.Id, b.Id)))
#             edges.setdefault(key, []).append((element, a, b))
#     boundary = [entries[0] for entries in edges.values() if len(entries) == 1]
#     if not boundary or any(len(entries) > 2 for entries in edges.values()):
#         raise RuntimeError("Inner Q1 mesh has empty or non-manifold boundary.")
#     values, interpolated, extrapolated, max_distance = {}, 0, 0, 0.0
#     warnings = []
#     for node in nodes:
#         for k in (1, 2):
#             x, y = node.X + speed[0]*k*dt, node.Y + speed[1]*k*dt
#             found, shape, element = locator.FindPointOnMesh([x, y, 0.0])
#             if found:
#                 interpolated += 1
#             else:
#                 nearest = None
#                 for parent, a, b in boundary:
#                     dx, dy = b.X-a.X, b.Y-a.Y
#                     length2 = dx*dx+dy*dy
#                     if length2 == 0.0:
#                         raise RuntimeError("Zero-length inner boundary edge.")
#                     projection = max(0.0, min(1.0, ((x-a.X)*dx+(y-a.Y)*dy)/length2))
#                     distance = math.hypot(x-a.X-projection*dx, y-a.Y-projection*dy)
#                     if nearest is None or (distance, parent.Id) < (nearest[0], nearest[1].Id):
#                         nearest = (distance, parent)
#                 distance, element = nearest
#                 parent_nodes = list(element.GetGeometry())
#                 edge_lengths = [math.hypot(parent_nodes[(i+1)%4].X-parent_nodes[i].X,
#                                             parent_nodes[(i+1)%4].Y-parent_nodes[i].Y) for i in range(4)]
#                 h = min(edge_lengths)
#                 if h <= 0.0 or distance > 0.10*h:
#                     raise RuntimeError(f"History extrapolation exceeds 0.10*h: node={node.Id}, "
#                         f"step={k}, time={time-k*dt}, point=({node.X}, {node.Y}), "
#                         f"distance={distance}, h={h}, parent={element.Id}.")
#                 # Affine Q1 coordinates of the unprojected point. Never clamp
#                 # s or t: values outside [0, 1] are the intended extrapolation.
#                 a, b, d = parent_nodes[0], parent_nodes[1], parent_nodes[3]
#                 ex, ey, fx, fy = b.X-a.X, b.Y-a.Y, d.X-a.X, d.Y-a.Y
#                 determinant = ex*fy-ey*fx
#                 if abs(determinant) <= 1e-30:
#                     raise RuntimeError("Degenerate inner Q1 parent.")
#                 s = ((x-a.X)*fy-(y-a.Y)*fx)/determinant
#                 t = (ex*(y-a.Y)-ey*(x-a.X))/determinant
#                 shape = [(1-s)*(1-t), s*(1-t), s*t, (1-s)*t]
#                 extrapolated += 1
#                 max_distance = max(max_distance, distance)
#                 warnings.append((node.Id, k, time-k*dt, distance, h, element.Id))
#             source = list(element.GetGeometry())
#             values[node.Id, k] = [sum(shape[i]*source[i].GetSolutionStepValue(KM.VELOCITY, k)[d]
#                                        for i in range(4)) for d in range(3)]
#     # Report accepted extrapolations only after every sample has passed.
#     for node_id, step, historical_time, distance, h, parent_id in warnings:
#         KM.Logger.PrintWarning("PatchHistoryTransfer", f"Extrapolated node {node_id}, step {step}, "
#             f"time {historical_time}, distance {distance}, h {h}, parent {parent_id}.")
#     # Commit only after every sample has passed the distance check.
#     for node in nodes:
#         for k in (1, 2):
#             node.SetSolutionStepValue(KM.VELOCITY, k, values[node.Id, k])
#         node.SetSolutionStepValue(KM.VELOCITY, values[node.Id, 1])
#     return interpolated, extrapolated, max_distance
#


class OuterHistoryValidity:
    """Track solved outer-node support alongside a three-slot Kratos buffer.

    Call advance immediately after CloneTimeStep and mark_solved only after a
    successful solve. Initial slots must only include supplied valid data.
    """
    def __init__(self, current_nodes=(), previous_nodes=()):
        self.slots = [set(current_nodes), set(previous_nodes), set()]

    def advance(self):
        self.slots = [set(), self.slots[0], self.slots[1]]

    def mark_solved(self, node_ids):
        self.slots[0] = set(node_ids)


class GapHistory:
    """Two historical snapshots of gap triangles and their Q1 parent fields."""
    def __init__(self):
        self.slots = [None, None, None]

    def advance(self):
        self.slots = [None, self.slots[0], self.slots[1]]

    def capture(self, interface, time):
        volumes = {g.Id: g for g in interface.GetSubModelPart("GapVolumes").Geometries}
        triangles = list(interface.GetSubModelPart("GapTriangles").Geometries)
        if len(volumes) != len(triangles):
            raise RuntimeError("Gap triangle and quadrature geometry counts differ.")
        cells = []
        for triangle in triangles:
            # CreateGeometry assigns each quadrature geometry immediately
            # before its matching visualization triangle.
            if triangle.Id-1 not in volumes:
                raise RuntimeError("Gap triangle has no matching parent quadrature geometry.")
            parent = list(volumes[triangle.Id-1])
            if len(parent) != 4:
                raise RuntimeError("Gap history requires Cartesian Q1 parent nodes.")
            cells.append((tuple((n.X, n.Y) for n in triangle),
                          tuple((n.Id, n.X, n.Y, tuple(n.GetSolutionStepValue(KM.VELOCITY))) for n in parent)))
        self.slots[0] = (time, cells)

    def interpolate(self, slot, time, x, y):
        snapshot = self.slots[slot]
        if snapshot is None:
            return None
        recorded_time, cells = snapshot
        if not math.isclose(recorded_time, time, rel_tol=1e-12, abs_tol=1e-12):
            raise RuntimeError(f"Gap history time mismatch: requested {time}, stored {recorded_time}.")
        candidates = []
        for triangle, parent in cells:
            def cross(a, b, p):
                return (b[0]-a[0])*(p[1]-a[1])-(b[1]-a[1])*(p[0]-a[0])
            area = cross(triangle[0], triangle[1], triangle[2])
            if area <= 0.0:
                raise RuntimeError("Invalid historical gap triangle orientation.")
            bary = [cross(triangle[(i+1)%3], triangle[(i+2)%3], (x, y))/area for i in range(3)]
            if min(bary) < -1e-12:
                continue
            a, b, d = parent[0], parent[1], parent[3]
            ex, ey, fx, fy = b[1]-a[1], b[2]-a[2], d[1]-a[1], d[2]-a[2]
            det = ex*fy-ey*fx
            if abs(det) < 1e-30:
                raise RuntimeError("Degenerate historical gap parent.")
            s = ((x-a[1])*fy-(y-a[2])*fx)/det
            t = (ex*(y-a[2])-ey*(x-a[1]))/det
            shape = ((1-s)*(1-t), s*(1-t), s*t, (1-s)*t)
            velocity = [sum(shape[i]*parent[i][3][d] for i in range(4)) for d in range(3)]
            candidates.append((min(bary), tuple(n[0] for n in parent), velocity))
        if not candidates:
            return None
        # On an internal interface distinct extensions may be discontinuous;
        # avoid arbitrarily selecting a trace when the target is exactly there.
        interior = [c for c in candidates if c[0] > 1e-12]
        if len({c[1] for c in (interior or candidates)}) > 1:
            raise RuntimeError(f"Ambiguous historical gap interface at ({x}, {y}), time {time}.")
        return (interior or candidates)[0][2]


def transfer_background_histories(inner, nodes, time, dt, speed, validity=None, required_steps=(1, 2), gap_history=None, observer=None, *, inner_patches=None, velocities=None, gap_histories=None, displacement_offsets=None):
    """Reuse valid outer histories, otherwise interpolate inside the ALE patch.

    Historical gap snapshots supply the parent Q1 extension when available.
    Values are committed only if every required sample has a numerical source.
    """
    if isinstance(speed, (int, float)):
        speed = (speed, 0.0)
    if any(k not in (1, 2) for k in required_steps):
        raise ValueError("History slots must be 1 or 2.")
    patches = [inner] if inner_patches is None else list(inner_patches)
    speeds = [speed] if velocities is None else list(velocities)
    histories = [gap_history] if gap_histories is None else list(gap_histories)
    if len(patches) != len(speeds) or len(patches) != len(histories):
        raise ValueError("Each patch requires a velocity and gap history.")
    bounds = []
    for patch in patches:
        coordinates = [(node.X, node.Y) for node in patch.Nodes]
        if not coordinates:
            raise ValueError('Historical inner patch has no nodes.')
        xs, ys = zip(*coordinates)
        scale = max(1.0, *(abs(value) for point in coordinates for value in point))
        tolerance = 1e-9 * scale
        bounds.append((min(xs)-tolerance, max(xs)+tolerance,
                       min(ys)-tolerance, max(ys)+tolerance))
    # A locator must not return values from cells deactivated inside a body.
    search_parts = []
    for patch in patches:
        if any(element.IsDefined(KM.ACTIVE) and element.IsNot(KM.ACTIVE) for element in patch.Elements):
            name = 'HistoryTransferActiveElements'
            if patch.HasSubModelPart(name):
                patch.RemoveSubModelPart(name)
            search = patch.CreateSubModelPart(name)
            search.AddElements([element.Id for element in patch.Elements
                                if not element.IsDefined(KM.ACTIVE) or element.Is(KM.ACTIVE)])
            search.AddNodes(sorted({node.Id for element in search.Elements for node in element.GetGeometry()}))
            search_parts.append(search)
        else:
            search_parts.append(patch)
    locators = [KM.BinBasedFastPointLocator2D(patch) for patch in search_parts]
    for locator in locators:
        locator.UpdateSearchDatabase()
    values, missing = {}, []
    counts = {"outer": 0, "inner": 0, "gap": 0, "extrapolation": 0, "missing": 0}
    for node in nodes:
        for k in required_steps:
            if validity is not None and node.Id in validity.slots[k]:
                counts["outer"] += 1
                if observer is not None:
                    observer(node, k, "outer", None, None)
                continue
            matches = []
            for index, (locator, patch_speed, (xmin, xmax, ymin, ymax)) in enumerate(zip(locators, speeds, bounds)):
                offset = (patch_speed[0]*k*dt, patch_speed[1]*k*dt) if displacement_offsets is None else displacement_offsets[index][k]
                x = node.X + offset[0]
                y = node.Y + offset[1]
                # The current mesh bounds enclose the mapped historical query,
                # including rotated patches; the locator performs the exact test.
                if not (xmin <= x <= xmax and ymin <= y <= ymax):
                    continue
                found, shape, element = locator.FindPointOnMesh([x, y, 0.0])
                if found:
                    matches.append((shape, element))
            if len(matches) > 1:
                raise RuntimeError(f"Ambiguous historical inner patches for node {node.Id}, slot {k}.")
            if not matches:
                gap_values = [value for history in histories if history is not None
                              for value in [history.interpolate(k, time-k*dt, node.X, node.Y)]
                              if value is not None]
                if len(gap_values) > 1:
                    raise RuntimeError(f"Ambiguous historical gap patches for node {node.Id}, slot {k}.")
                if gap_values:
                    values[node.Id, k] = gap_values[0]
                    counts["gap"] += 1
                    KM.Logger.PrintWarning("PatchHistoryTransfer",
                        f"Reconstructing historical velocity from the gap for node {node.Id}, "
                        f"slot {k}, time {time-k*dt}, coordinates ({node.X}, {node.Y}). "
                        "Investigate this gap-history sample.")
                    if observer is not None:
                        observer(node, k, "gap", None, None)
                    continue
                missing.append((node.Id, k, time-k*dt, node.X, node.Y))
                counts["missing"] += 1
                if observer is not None:
                    observer(node, k, "missing", None, None)
                continue
            shape, element = matches[0]
            source = list(element.GetGeometry())
            values[node.Id, k] = [sum(shape[i]*source[i].GetSolutionStepValue(KM.VELOCITY, k)[d]
                                       for i in range(4)) for d in range(3)]
            counts["inner"] += 1
            if observer is not None:
                observer(node, k, "inner", element.Id, list(shape))
    if missing:
        raise RuntimeError(f"Missing {counts['missing']} required velocity histories (source counts: {counts}) "
            f"(node, buffer slot, time, x, y): {missing}. Outer history was not valid; "
            "point is outside the historical inner patch and available gap snapshots; "
            "extrapolation is disabled.")
    for node in nodes:
        for k in required_steps:
            if (node.Id, k) in values:
                node.SetSolutionStepValue(KM.VELOCITY, k, values[node.Id, k])
        if 1 in required_steps and (validity is None or node.Id not in validity.slots[1]):
            node.SetSolutionStepValue(KM.VELOCITY, node.GetSolutionStepValue(KM.VELOCITY, 1))
    return counts
