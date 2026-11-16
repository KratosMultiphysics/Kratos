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
