"""Simplified Bishop method for circular slip surfaces on a 2D triangle mesh.

Pure numpy, so every step can be tested without Kratos. The Kratos process that uses it is
bishop_slope_stability_process.py.

Dry slopes only: no pore pressures are considered.

Conventions:
- A circle is (xc, yc, R); the slip surface is its lower branch y_base(x) = yc - sqrt(R^2 - (x - xc)^2).
- The level set of a circle is phi = |x - c| - R (negative inside the circle).
"""

from dataclasses import dataclass

import numpy as np

# Points closer than this (relative to the circle radius) are treated as coincident
_RELATIVE_TOLERANCE = 1.0e-9


@dataclass
class TriangleMesh:
    """Linear triangles used for the geometric operations.

    coordinates:      (N, 2) nodal coordinates
    triangles:        (E, 3) node indices of the corners of each triangle
    edges:            (M, 2) node indices of the unique edges
    edge_triangles:   (M, 2) triangles adjacent to each edge (-1 for the second one of a boundary edge)
    """
    coordinates: np.ndarray
    triangles: np.ndarray
    edges: np.ndarray
    edge_triangles: np.ndarray

    @property
    def is_boundary_edge(self):
        return self.edge_triangles[:, 1] < 0


def build_triangle_mesh(coordinates, triangles):
    coordinates = np.asarray(coordinates, dtype=float)
    triangles = np.asarray(triangles, dtype=int)

    local_edges = np.array([[0, 1], [1, 2], [2, 0]])
    all_edges = np.sort(triangles[:, local_edges].reshape(-1, 2), axis=1)
    owners = np.repeat(np.arange(len(triangles)), 3)
    edges, inverse = np.unique(all_edges, axis=0, return_inverse=True)

    edge_triangles = np.full((len(edges), 2), -1, dtype=int)
    for edge_index, owner in zip(inverse.ravel(), owners):
        slot = 0 if edge_triangles[edge_index, 0] < 0 else 1
        edge_triangles[edge_index, slot] = owner
    return TriangleMesh(coordinates, triangles, edges, edge_triangles)


def circle_level_set(coordinates, center, radius):
    return np.hypot(coordinates[:, 0] - center[0], coordinates[:, 1] - center[1]) - radius


def lower_arc_y(x, center, radius):
    return center[1] - np.sqrt(np.maximum(radius ** 2 - (np.asarray(x) - center[0]) ** 2, 0.0))


def edge_cuts(mesh, phi):
    """Edges whose end nodes have a different sign of phi, and the linearly interpolated cut points."""
    phi1, phi2 = phi[mesh.edges[:, 0]], phi[mesh.edges[:, 1]]
    cut = phi1 * phi2 < 0.0
    t = phi1[cut] / (phi1[cut] - phi2[cut])
    p1, p2 = mesh.coordinates[mesh.edges[cut, 0]], mesh.coordinates[mesh.edges[cut, 1]]
    return np.flatnonzero(cut), p1 + t[:, None] * (p2 - p1)


def exact_edge_circle_intersection(p1, p2, center, radius):
    """Intersection of segment p1-p2 with the circle; the segment is known to cross it once."""
    d = p2 - p1
    f = p1 - np.asarray(center)
    a, b, c = d @ d, 2.0 * (f @ d), f @ f - radius ** 2
    roots = (-b + np.array([-1.0, 1.0]) * np.sqrt(max(b * b - 4.0 * a * c, 0.0))) / (2.0 * a)
    t = roots[np.argmin(np.abs(roots - np.clip(roots, 0.0, 1.0)))]
    return p1 + np.clip(t, 0.0, 1.0) * d


def barycentric_coordinates(mesh, points):
    """(K, E, 3) barycentric coordinates of K points with respect to all E triangles."""
    corners = mesh.coordinates[mesh.triangles]  # (E, 3, 2)
    a, b, c = corners[:, 0], corners[:, 1], corners[:, 2]
    det = (b[:, 0] - a[:, 0]) * (c[:, 1] - a[:, 1]) - (c[:, 0] - a[:, 0]) * (b[:, 1] - a[:, 1])
    px, py = points[:, None, 0], points[:, None, 1]
    l1 = ((px - a[:, 0]) * (c[:, 1] - a[:, 1]) - (c[:, 0] - a[:, 0]) * (py - a[:, 1])) / det
    l2 = ((b[:, 0] - a[:, 0]) * (py - a[:, 1]) - (px - a[:, 0]) * (b[:, 1] - a[:, 1])) / det
    return np.stack([1.0 - l1 - l2, l1, l2], axis=-1)


def locate_points(mesh, points, tolerance=1.0e-9):
    """Index of a triangle containing each point (-1 if outside) and its barycentric coordinates."""
    points = np.atleast_2d(points)
    bary = barycentric_coordinates(mesh, points)
    inside = np.all(bary >= -tolerance, axis=-1)
    found = inside.any(axis=1)
    index = np.where(found, np.argmax(inside, axis=1), -1)
    return index, bary[np.arange(len(points)), np.maximum(index, 0)]


@dataclass
class SlipArc:
    center: np.ndarray
    radius: float
    entry: np.ndarray        # point where the arc leaves the boundary at the smallest x
    exit: np.ndarray         # point where the arc leaves the boundary at the largest x
    cut_edges: np.ndarray    # indices of all edges crossed by the circle
    cut_points: np.ndarray   # (K, 2) linearly interpolated cut points


def find_slip_arc(mesh, center, radius, forbidden_nodes=None):
    """Returns the SlipArc of the circle, or None if the circle is not an admissible slip surface.

    Admissible means: the circle cuts the mesh, no cut triangle has a forbidden node, it crosses
    the boundary exactly twice (so the part inside the domain is a single arc) and both crossings
    are on the lower half of the circle (so the arc is a graph y_base(x) over [x_entry, x_exit]).
    """
    center = np.asarray(center, dtype=float)
    phi = circle_level_set(mesh.coordinates, center, radius)
    # Nodes exactly on the circle are moved inside, so that every crossing is a sign change
    phi = np.where(np.abs(phi) < _RELATIVE_TOLERANCE * radius, -_RELATIVE_TOLERANCE * radius, phi)

    cut_edges, cut_points = edge_cuts(mesh, phi)
    if len(cut_edges) == 0:
        return None

    if forbidden_nodes is not None:
        triangle_phi = phi[mesh.triangles]
        cut_triangles = (triangle_phi.min(axis=1) < 0.0) & (triangle_phi.max(axis=1) > 0.0)
        if np.any(forbidden_nodes[mesh.triangles[cut_triangles]]):
            return None

    boundary_cuts = cut_edges[mesh.is_boundary_edge[cut_edges]]
    if len(boundary_cuts) != 2:
        return None

    crossings = np.array([exact_edge_circle_intersection(mesh.coordinates[mesh.edges[e, 0]],
                                                         mesh.coordinates[mesh.edges[e, 1]], center, radius)
                          for e in boundary_cuts])
    if np.any(crossings[:, 1] >= center[1]):
        return None
    entry, exit_ = crossings[np.argsort(crossings[:, 0])]
    if exit_[0] - entry[0] < _RELATIVE_TOLERANCE * radius:
        return None

    # The lower branch between the crossings must be the part of the circle inside the domain
    x_mid = 0.5 * (entry[0] + exit_[0])
    probe = np.array([[x_mid, lower_arc_y(x_mid, center, radius) + 1.0e-6 * radius]])
    if locate_points(mesh, probe)[0][0] < 0:
        return None

    return SlipArc(center, radius, entry, exit_, cut_edges, cut_points)


def surface_breakpoints(mesh, excluded_nodes=None, angle_tolerance=1.0e-6):
    """x of the boundary nodes (not in excluded_nodes) where the direction of the boundary changes."""
    boundary = mesh.edges[mesh.is_boundary_edge]
    if excluded_nodes is not None:
        boundary = boundary[~np.all(excluded_nodes[boundary], axis=1)]
    neighbours = {}
    for a, b in boundary:
        neighbours.setdefault(a, []).append(b)
        neighbours.setdefault(b, []).append(a)

    breaks = []
    for node, adjacent in neighbours.items():
        if len(adjacent) != 2:
            continue
        p = mesh.coordinates[node]
        d1, d2 = mesh.coordinates[adjacent[0]] - p, mesh.coordinates[adjacent[1]] - p
        sine = abs(d1[0] * d2[1] - d1[1] * d2[0]) / (np.linalg.norm(d1) * np.linalg.norm(d2))
        if sine > angle_tolerance:
            breaks.append(p[0])
    return np.array(sorted(breaks))


def base_breakpoints(mesh, arc, triangle_material):
    """x of the points of the slip arc where the material changes."""
    adjacent = mesh.edge_triangles[arc.cut_edges]
    interior = adjacent[:, 1] >= 0
    different = interior & (triangle_material[adjacent[:, 0]] != triangle_material[np.maximum(adjacent[:, 1], 0)])
    return np.sort(arc.cut_points[different, 0])


def slice_boundaries(x_entry, x_exit, number_of_slices, extra_breakpoints=(), merge_tolerance=1.0e-6):
    """Uniform division of [x_entry, x_exit] plus the extra breakpoints that lie inside it."""
    width = x_exit - x_entry
    points = np.concatenate([np.linspace(x_entry, x_exit, number_of_slices + 1), np.asarray(extra_breakpoints, float)])
    points = np.sort(points[(points >= x_entry) & (points <= x_exit)])
    keep = np.concatenate([[True], np.diff(points) > merge_tolerance * width])
    points = points[keep]
    points[-1] = x_exit
    return points


def column_weights(mesh, triangle_unit_weight, xs, y_base):
    """p(x) = sum over the triangles cut by the vertical x of gamma * (length above y_base).

    xs and y_base are arrays of the same shape; the result has that shape as well.
    """
    xs = np.asarray(xs, dtype=float)
    shape = xs.shape
    xs, y_base = xs.ravel(), np.asarray(y_base, dtype=float).ravel()

    corners = mesh.coordinates[mesh.triangles]  # (E, 3, 2)
    in_range = (corners[:, :, 0].max(axis=1) > xs.min()) & (corners[:, :, 0].min(axis=1) < xs.max())
    corners, gamma = corners[in_range], np.asarray(triangle_unit_weight)[in_range]

    xa, ya = corners[:, :, 0], corners[:, :, 1]
    xb, yb = np.roll(xa, -1, axis=1), np.roll(ya, -1, axis=1)
    x = xs[:, None, None]
    # Half-open test, so a vertical through a vertex counts each crossing once
    crosses = ((xa <= x) & (x < xb)) | ((xb <= x) & (x < xa))
    dx = np.where(xb != xa, xb - xa, 1.0)
    y = ya + (x - xa) * (yb - ya) / dx

    y_low = np.where(crosses, y, np.inf).min(axis=2)
    y_high = np.where(crosses, y, -np.inf).max(axis=2)
    length = np.where(crosses.sum(axis=2) >= 2, y_high - np.maximum(y_low, y_base[:, None]), 0.0)
    return (np.clip(length, 0.0, None) @ gamma).reshape(shape)


@dataclass
class Slices:
    x_left: np.ndarray
    x_right: np.ndarray
    width: np.ndarray        # b
    x_mid: np.ndarray
    y_mid: np.ndarray        # base midpoint
    alpha: np.ndarray        # base inclination, sin(alpha) = (x_mid - xc) / R
    base_length: np.ndarray  # l
    weight: np.ndarray       # W
    cohesion: np.ndarray     # c'
    friction_angle: np.ndarray  # phi' [rad]


def build_slices(center, radius, boundaries, column_weight, base_properties, gauss_points_per_slice=4):
    """column_weight(xs, y_base) -> p; base_properties(x_mid, y_mid) -> (c', phi' [rad])."""
    center = np.asarray(center, dtype=float)
    x_left, x_right = boundaries[:-1], boundaries[1:]
    width = x_right - x_left

    gauss_x, gauss_w = np.polynomial.legendre.leggauss(gauss_points_per_slice)
    xs = 0.5 * (x_left + x_right)[:, None] + 0.5 * width[:, None] * gauss_x[None, :]
    p = column_weight(xs, lower_arc_y(xs, center, radius))
    weight = 0.5 * width * (p @ gauss_w)

    theta = np.arctan2(lower_arc_y(boundaries, center, radius) - center[1], boundaries - center[0])
    base_length = radius * np.abs(np.diff(theta))

    x_mid = 0.5 * (x_left + x_right)
    y_mid = lower_arc_y(x_mid, center, radius)
    alpha = np.arcsin(np.clip((x_mid - center[0]) / radius, -1.0, 1.0))
    cohesion, friction_angle = base_properties(x_mid, y_mid)
    return Slices(x_left, x_right, width, x_mid, y_mid, alpha, base_length, weight,
                  np.asarray(cohesion, float), np.asarray(friction_angle, float))


@dataclass
class BishopResult:
    factor_of_safety: float
    fellenius_factor_of_safety: float
    iterations: int
    m_alpha: np.ndarray
    # +1 if sum(W sin(alpha)) > 0: the mass rotates clockwise about the center, so it slides towards -x
    direction: int


def bishop_factor_of_safety(slices, tolerance=1.0e-4, max_iterations=50, min_m_alpha=0.2):
    """Simplified Bishop, started from Fellenius. Returns None if the circle must be discarded."""
    driving = np.sum(slices.weight * np.sin(slices.alpha))
    if abs(driving) <= 1.0e-12 * max(np.sum(slices.weight), 1.0e-300):
        return None
    # alpha is measured positive in the direction of sliding
    direction = 1 if driving > 0.0 else -1
    alpha = direction * slices.alpha
    driving = abs(driving)

    tan_phi = np.tan(slices.friction_angle)
    w, b, l, c = slices.weight, slices.width, slices.base_length, slices.cohesion

    fellenius = np.sum(c * l + w * np.cos(alpha) * tan_phi) / driving
    numerator_terms = c * b + w * tan_phi

    factor = fellenius if fellenius > 0.0 else 1.0
    for iteration in range(1, max_iterations + 1):
        m_alpha = np.cos(alpha) + np.sin(alpha) * tan_phi / factor
        if np.any(m_alpha <= 0.0):
            return None
        new_factor = np.sum(numerator_terms / m_alpha) / driving
        if abs(new_factor - factor) < tolerance:
            factor = new_factor
            m_alpha = np.cos(alpha) + np.sin(alpha) * tan_phi / factor
            if np.min(m_alpha) < min_m_alpha:
                return None
            return BishopResult(factor, fellenius, iteration, m_alpha, direction)
        factor = new_factor
    return None


def evaluate_circle(mesh, center, radius, number_of_slices, column_weight, base_properties, forbidden_nodes=None,
                    extra_breakpoints=None, gauss_points_per_slice=4, tolerance=1.0e-4, max_iterations=50,
                    min_m_alpha=0.2):
    """Full evaluation of one circle. Returns (BishopResult, Slices, SlipArc) or None.

    extra_breakpoints(arc) -> x of additional slice boundaries for that arc.
    """
    arc = find_slip_arc(mesh, center, radius, forbidden_nodes)
    if arc is None:
        return None
    extra = extra_breakpoints(arc) if extra_breakpoints is not None else ()
    boundaries = slice_boundaries(arc.entry[0], arc.exit[0], number_of_slices, extra)
    slices = build_slices(arc.center, radius, boundaries, column_weight, base_properties, gauss_points_per_slice)
    result = bishop_factor_of_safety(slices, tolerance, max_iterations, min_m_alpha)
    if result is None:
        return None
    return result, slices, arc


def grid_search(evaluate, centers_x, centers_y, radii):
    """evaluate(xc, yc, R) -> factor of safety or None. Returns (min FoS, (xc, yc, R), number of valid circles)."""
    best_factor, best_circle, valid = np.inf, None, 0
    for xc in centers_x:
        for yc in centers_y:
            for radius in radii:
                factor = evaluate(xc, yc, radius)
                if factor is None:
                    continue
                valid += 1
                if factor < best_factor:
                    best_factor, best_circle = factor, (float(xc), float(yc), float(radius))
    return best_factor, best_circle, valid


def refine_with_nelder_mead(evaluate, circle, penalty=1.0e6):
    """Local minimisation of the factor of safety around circle (needs scipy)."""
    from scipy.optimize import minimize

    def objective(v):
        factor = evaluate(v[0], v[1], v[2])
        return penalty if factor is None else factor

    result = minimize(objective, np.asarray(circle, float), method="Nelder-Mead",
                      options={"xatol": 1.0e-3, "fatol": 1.0e-5, "maxiter": 400})
    factor = evaluate(*result.x)
    if factor is None:
        return None
    return factor, tuple(float(v) for v in result.x)
