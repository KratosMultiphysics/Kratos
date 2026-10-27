"""Generates the scripted figures of the documentation pages.

Every figure here is produced from the application's own numpy-only bridge
paths (tessellation, graph extraction, particle graphs, grid sampling, the
benchmark) or from synthetic numpy data, with matplotlib - no torch, no
physicsnemo, no GPU. Hand-authored SVG diagrams live next to the outputs and
are not generated.

    python3 make_docs_figures.py                 # all, into the docs tree
    python3 make_docs_figures.py --list
    python3 make_docs_figures.py --only calibration_views,halo_partition --out /tmp/figs
    python3 make_docs_figures.py --quick         # smaller benchmark run

Each figure is independent: a failure is reported and the others still run.
"""

import argparse
import re
import subprocess
import sys
from pathlib import Path

import numpy

_HERE = Path(__file__).resolve().parent
_APP_DIR = _HERE.parent.parent
_DOCS_DIR = (_APP_DIR.parent.parent / "docs" / "pages" / "Applications"
             / "PhysicsNeMo_Application")

KRATOS_BLUE, PNEMO_GREEN, ACCENT, PINK, GREY = "#1f4e79", "#2e6b2e", "#c27c0e", "#a03060", "#666666"


def _Matplotlib():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as pyplot
    pyplot.rcParams.update({"font.size": 10, "axes.titlesize": 11, "figure.dpi": 130})
    return pyplot


def _Kratos():
    import KratosMultiphysics
    return KratosMultiphysics


# --------------------------------------------------------------------------- helpers

def _TriangleGrid(nx, ny, extent=(1.0, 1.0), jitter=0.0, seed=0):
    """Node coordinates (N, 2) and triangles (T, 3) of a structured triangle mesh."""
    rng = numpy.random.default_rng(seed)
    xs, ys = numpy.linspace(0.0, extent[0], nx + 1), numpy.linspace(0.0, extent[1], ny + 1)
    X, Y = numpy.meshgrid(xs, ys, indexing="ij")
    points = numpy.stack([X.ravel(), Y.ravel()], axis=1)
    if jitter:
        interior = (points[:, 0] > 0) & (points[:, 0] < extent[0]) & (points[:, 1] > 0) & (points[:, 1] < extent[1])
        points[interior] += rng.uniform(-jitter, jitter, size=(interior.sum(), 2))
    def node(i, j): return i * (ny + 1) + j
    triangles = []
    for i in range(nx):
        for j in range(ny):
            a, b, c, d = node(i, j), node(i + 1, j), node(i + 1, j + 1), node(i, j + 1)
            triangles += [(a, b, c), (a, c, d)]
    return points, numpy.array(triangles, dtype=numpy.int64)


def _TriangleModelPart(model, points, triangles, name="Figure", variables=()):
    Kratos = _Kratos()
    model_part = model.CreateModelPart(name)
    model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 2
    for variable in variables:
        model_part.AddNodalSolutionStepVariable(variable)
    model_part.SetBufferSize(1)
    for index, (x, y) in enumerate(points):
        model_part.CreateNewNode(index + 1, float(x), float(y), 0.0)
    properties = model_part.CreateNewProperties(1)
    for index, tri in enumerate(triangles):
        model_part.CreateNewElement("Element2D3N", index + 1, [int(n) + 1 for n in tri], properties)
    return model_part


# --------------------------------------------------------------------------- figures

def make_tessellation_modes(out):
    """Corner tessellation of the four non-simplex geometries by the smallest-id rule."""
    Kratos = _Kratos()
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    from KratosMultiphysics.PhysicsNeMoApplication.bridges.mesh_bridge import tessellation
    pyplot = _Matplotlib()
    GT = Kratos.GeometryData.KratosGeometryType
    hexa = numpy.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]], float)
    prism = numpy.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [0, 1, 1]], float)
    pyramid = numpy.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0.5, 0.5, 1]], float)
    # two hexahedron numberings: the smallest-id rule yields 6 or 5 tetrahedra
    rng = numpy.random.default_rng(3)
    six_ids, five_ids = None, None
    for _ in range(200):
        ids = [int(i) for i in rng.permutation(8) + 1]
        count = len(tessellation.TessellateEntity(GT.Kratos_Hexahedra3D8, ids, {i: hexa[k] for k, i in enumerate(ids)}))
        if count == 6 and six_ids is None: six_ids = ids
        if count == 5 and five_ids is None: five_ids = ids
        if six_ids and five_ids: break
    cases = [("hexahedron, ids give 6 tetrahedra", GT.Kratos_Hexahedra3D8, hexa, six_ids or list(range(1, 9)), [[0, 1, 2, 3], [4, 5, 6, 7], [0, 1, 5, 4], [1, 2, 6, 5], [2, 3, 7, 6], [3, 0, 4, 7]]),
             ("hexahedron, ids give 5 tetrahedra", GT.Kratos_Hexahedra3D8, hexa, five_ids or list(range(1, 9)), [[0, 1, 2, 3], [4, 5, 6, 7], [0, 1, 5, 4], [1, 2, 6, 5], [2, 3, 7, 6], [3, 0, 4, 7]]),
             ("prism, 3 tetrahedra", GT.Kratos_Prism3D6, prism, [4, 1, 6, 2, 5, 3], [[0, 1, 2], [3, 4, 5], [0, 1, 4, 3], [1, 2, 5, 4], [2, 0, 3, 5]]),
             ("pyramid, 2 tetrahedra", GT.Kratos_Pyramid3D5, pyramid, [3, 1, 5, 2, 4], [[0, 1, 2, 3], [0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]])]
    fig = pyplot.figure(figsize=(13, 6.6))
    for column, (title, gtype, corners, ids, faces) in enumerate(cases):
        coords = {node_id: corners[k] for k, node_id in enumerate(ids)}
        simplices = tessellation.TessellateEntity(gtype, ids, coords)
        rows = {node_id: k for k, node_id in enumerate(ids)}
        ax = fig.add_subplot(2, 4, column + 1, projection="3d")
        ax.add_collection3d(Poly3DCollection([corners[f] for f in faces], facecolors="#eaf2fb", edgecolors=KRATOS_BLUE, linewidths=1.2, alpha=0.35))
        for k, node_id in enumerate(ids):
            ax.text(*corners[k], str(node_id), color=KRATOS_BLUE, fontsize=9, fontweight="bold")
        ax.set_title(f"Kratos {title.split(',')[0]}\n(node ids as labelled)", color=KRATOS_BLUE)
        ax2 = fig.add_subplot(2, 4, column + 5, projection="3d")
        cmap = pyplot.get_cmap("Set2")
        for t, simplex in enumerate(simplices):
            pts = numpy.array([corners[rows[n]] for n in simplex])
            centre = pts.mean(axis=0)
            shrunk = centre + 0.82 * (pts - centre)
            tri_faces = [shrunk[[0, 1, 2]], shrunk[[0, 1, 3]], shrunk[[0, 2, 3]], shrunk[[1, 2, 3]]]
            ax2.add_collection3d(Poly3DCollection(tri_faces, facecolors=cmap(t % 8), edgecolors="#333333", linewidths=0.6, alpha=0.85))
        ax2.set_title(f"{len(simplices)} tetrahedra\n(shrunk apart to show the split)", color=PNEMO_GREEN)
        for a in (ax, ax2):
            a.set_xlim(0, 1); a.set_ylim(0, 1); a.set_zlim(0, 1); a.set_axis_off(); a.view_init(elev=22, azim=-58)
    fig.suptitle('Corner tessellation, tessellation_mode = "smallest_id_diagonal": every quadrilateral face is split along the diagonal through its smallest node id, so neighbours agree', fontsize=10)
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    pyplot.close(fig)


def make_higher_order_modes(out):
    """reduce / subdivide / curved on a curved Quadrilateral2D9 and Triangle2D6."""
    Kratos = _Kratos()
    from KratosMultiphysics.PhysicsNeMoApplication.bridges.mesh_bridge import tessellation, curved_tessellation
    pyplot = _Matplotlib()
    GT = Kratos.GeometryData.KratosGeometryType
    quad_local = numpy.array([[-1, -1], [1, -1], [1, 1], [-1, 1], [0, -1], [1, 0], [0, 1], [-1, 0], [0, 0]], float)
    quad_nodes = numpy.array([[0, 0], [1, 0], [1, 1], [0, 1], [0.5, -0.18], [1.18, 0.5], [0.5, 1.18], [-0.18, 0.5], [0.5, 0.5]], float)
    tri_local = numpy.array([[0, 0], [1, 0], [0, 1], [0.5, 0], [0.5, 0.5], [0, 0.5]], float)
    tri_nodes = numpy.array([[0, 0], [1, 0], [0, 1], [0.5, -0.16], [0.62, 0.62], [-0.16, 0.5]], float)
    def mapped(gtype, nodes, local):
        N = curved_tessellation.EvaluateShapeFunctions(gtype, numpy.asarray(local, float))
        N = numpy.asarray(N).reshape(len(local), -1)
        return N @ nodes
    def boundary(gtype, nodes, corners_local):
        segs = []
        for a, b in zip(corners_local, numpy.roll(corners_local, -1, axis=0)):
            t = numpy.linspace(0, 1, 40)[:, None]
            segs.append(mapped(gtype, nodes, a + (b - a) * t))
        return numpy.concatenate(segs)
    rows = [("Quadrilateral2D9", GT.Kratos_Quadrilateral2D9, quad_local, quad_nodes, quad_local[:4]),
            ("Triangle2D6", GT.Kratos_Triangle2D6, tri_local, tri_nodes, tri_local[:3])]
    fig, axes = pyplot.subplots(2, 3, figsize=(12.5, 8))
    for r, (name, gtype, local, nodes, corner_local) in enumerate(rows):
        ids = list(range(1, len(nodes) + 1))
        coords = {i: numpy.array([*nodes[k], 0.0]) for k, i in enumerate(ids)}
        exact = boundary(gtype, nodes, corner_local)
        for c, mode in enumerate(["reduce", "subdivide", "curved"]):
            ax = axes[r, c]
            ax.plot(exact[:, 0], exact[:, 1], color=KRATOS_BLUE, lw=2.0, label="exact geometry")
            if mode in ("reduce", "subdivide"):
                simplices = tessellation.TessellateEntity(gtype, ids, coords, higher_order_mode=mode)
                for simplex in simplices:
                    pts = nodes[[s - 1 for s in simplex] + [simplex[0] - 1]]
                    ax.fill(pts[:, 0], pts[:, 1], facecolor="#eef7ee", edgecolor=PNEMO_GREEN, lw=1.2, alpha=0.9)
                ax.set_title(f'{name}\n"{mode}": {len(simplices)} triangles')
            else:
                level = 2
                n = 2 ** level
                if gtype == GT.Kratos_Quadrilateral2D9:
                    u = numpy.linspace(-1, 1, n + 1)
                    lattice = numpy.array([[a, b] for a in u for b in u])
                    for a in u:
                        ax.plot(*mapped(gtype, nodes, numpy.array([[a, b] for b in numpy.linspace(-1, 1, 40)])).T, color=PNEMO_GREEN, lw=0.9)
                        ax.plot(*mapped(gtype, nodes, numpy.array([[b, a] for b in numpy.linspace(-1, 1, 40)])).T, color=PNEMO_GREEN, lw=0.9)
                    cells = n * n * 2
                else:
                    u = numpy.linspace(0, 1, n + 1)
                    lattice = numpy.array([[a, b] for a in u for b in u if a + b <= 1 + 1e-12])
                    for k in range(n + 1):
                        s = u[k]
                        ax.plot(*mapped(gtype, nodes, numpy.array([[s, b] for b in numpy.linspace(0, 1 - s, 40)])).T, color=PNEMO_GREEN, lw=0.9)
                        ax.plot(*mapped(gtype, nodes, numpy.array([[b, s] for b in numpy.linspace(0, 1 - s, 40)])).T, color=PNEMO_GREEN, lw=0.9)
                        ax.plot(*mapped(gtype, nodes, numpy.array([[b, s - b] for b in numpy.linspace(0, s, 40)])).T, color=PNEMO_GREEN, lw=0.9)
                    cells = n * n
                pts = mapped(gtype, nodes, lattice)
                real = numpy.array([numpy.any(numpy.all(numpy.isclose(lattice, l), axis=1)) for l in local])
                synthetic = numpy.array([not numpy.any(numpy.all(numpy.isclose(local, l), axis=1)) for l in lattice])
                ax.scatter(pts[synthetic, 0], pts[synthetic, 1], s=14, color=ACCENT, zorder=3, label="synthetic points (gather only)")
                ax.set_title(f'{name}\n"curved", levels = {level}: {cells} triangles')
            ax.scatter(nodes[:, 0], nodes[:, 1], s=34, color=KRATOS_BLUE, zorder=4, label="real Kratos nodes")
            ax.set_aspect("equal"); ax.set_xlim(-0.3, 1.3); ax.set_ylim(-0.3, 1.3); ax.set_xticks([]); ax.set_yticks([])
            if r == 0 and c == 2:
                ax.legend(loc="lower right", fontsize=8)
    fig.suptitle('Higher-order geometries: "reduce" drops the mid-side nodes, "subdivide" splits through them (straight edges), "curved" samples the exact geometry on a lattice', fontsize=10)
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    pyplot.close(fig)


def make_mesh_graph(out):
    """The element-edge graph BuildGraph extracts, and the bistride hierarchy on it."""
    Kratos = _Kratos()
    from KratosMultiphysics.PhysicsNeMoApplication.bridges import graph_bridge
    pyplot = _Matplotlib()
    points, triangles = _TriangleGrid(10, 6, extent=(1.6, 1.0), jitter=0.035, seed=1)
    model = Kratos.Model()
    model_part = _TriangleModelPart(model, points, triangles)
    _, edge_index, edge_features, node_ids = graph_bridge.BuildGraph(model_part)
    positions = graph_bridge.NodePositions(model_part)[:, :2]
    fig, axes = pyplot.subplots(1, 3, figsize=(14, 4.2))
    def draw_edges(ax, ei, pos, color, lw):
        a, b = ei
        for i, j in zip(a, b):
            if i < j:
                ax.plot([pos[i, 0], pos[j, 0]], [pos[i, 1], pos[j, 1]], color=color, lw=lw, zorder=1)
    ax = axes[0]
    ax.triplot(points[:, 0], points[:, 1], triangles, color=KRATOS_BLUE, lw=0.8)
    ax.scatter(points[:, 0], points[:, 1], s=12, color=KRATOS_BLUE, zorder=3)
    ax.set_title(f"Kratos mesh: {len(points)} nodes, {len(triangles)} triangles", color=KRATOS_BLUE)
    ax = axes[1]
    draw_edges(ax, edge_index, positions, PNEMO_GREEN, 0.8)
    ax.scatter(positions[:, 0], positions[:, 1], s=12, color=PNEMO_GREEN, zorder=3)
    ax.set_title(f"BuildGraph: edge_index (2, {edge_index.shape[1]}) bidirectional,\nedge_features (E, 4) = relative position + distance", color=PNEMO_GREEN)
    ax = axes[2]
    draw_edges(ax, edge_index, positions, "#bbbbbb", 0.6)
    ax.scatter(positions[:, 0], positions[:, 1], s=10, color="#bbbbbb", zorder=2, label="level 0 (all nodes)")
    try:
        hierarchy = graph_bridge.BuildBistrideHierarchy(edge_index, len(node_ids), graph_bridge.NodePositions(model_part), num_levels=2)
        ms_edges, ms_ids = hierarchy[0], hierarchy[1]
        level1 = numpy.asarray(ms_ids[0])
        ax.scatter(positions[level1, 0], positions[level1, 1], s=26, color=ACCENT, zorder=3, label=f"level 1 ({len(level1)} nodes)")
        if len(ms_ids) > 1:
            level2 = level1[numpy.asarray(ms_ids[1])]
            ax.scatter(positions[level2, 0], positions[level2, 1], s=60, color=PINK, zorder=4, marker="s", label=f"level 2 ({len(level2)} nodes)")
        ax.set_title("BuildBistrideHierarchy: BFS-parity coarsening\n(the multiscale levels of BiStrideMeshGraphNet)")
        ax.legend(loc="upper right", fontsize=8)
    except Exception as error:  # pragma: no cover - documented in the figure itself
        ax.set_title(f"bistride hierarchy unavailable: {type(error).__name__}")
    for ax in axes:
        ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    pyplot.close(fig)


def make_proximity_graph(out):
    """Radius, kNN and periodic radius proximity graphs over a particle cloud."""
    Kratos = _Kratos()
    from KratosMultiphysics.PhysicsNeMoApplication.bridges import particle_bridge
    pyplot = _Matplotlib()
    rng = numpy.random.default_rng(7)
    positions = numpy.zeros((70, 3))
    positions[:, :2] = rng.uniform(0.0, 1.0, size=(70, 2))
    cases = [("radius search, r = 0.16", '{"type": "radius", "radius": 0.16, "max_neighbors": 32, "backend": "numpy"}', False),
             ("k nearest neighbours, k = 4", '{"type": "knn", "radius": 0.16, "max_neighbors": 4, "backend": "numpy"}', False),
             ('periodic radius search, "box_size": [1, 1, 1]', '{"type": "radius", "radius": 0.16, "max_neighbors": 32, "backend": "numpy", "box_size": [1.0, 1.0, 1.0]}', True)]
    fig, axes = pyplot.subplots(1, 3, figsize=(14, 4.6))
    for ax, (title, settings, periodic) in zip(axes, cases):
        edge_index, edge_features = particle_bridge.BuildParticleGraphFromPositions(positions, Kratos.Parameters(settings))
        wrapped = 0
        for i, j in zip(*edge_index):
            if i < j or settings.find('"knn"') >= 0:
                d = positions[j] - positions[i]
                if periodic and numpy.any(numpy.abs(d[:2]) > 0.5):
                    wrapped += 1
                    ax.plot([positions[i, 0], positions[j, 0]], [positions[i, 1], positions[j, 1]], color=ACCENT, lw=0.5, ls=":", zorder=1)
                else:
                    ax.plot([positions[i, 0], positions[j, 0]], [positions[i, 1], positions[j, 1]], color=PNEMO_GREEN, lw=0.6, zorder=1)
        ax.scatter(positions[:, 0], positions[:, 1], s=16, color=KRATOS_BLUE, zorder=3)
        ax.add_patch(pyplot.Rectangle((0, 0), 1, 1, fill=False, ls="--", color="#999999"))
        extra = f", {wrapped} edges wrap through the box (dotted, minimum image)" if periodic else ""
        ax.set_title(f"{title}\n{edge_index.shape[1]} directed edges{extra}", fontsize=9.5)
        ax.set_aspect("equal"); ax.set_xlim(-0.05, 1.05); ax.set_ylim(-0.05, 1.05); ax.set_xticks([]); ax.set_yticks([])
    fig.suptitle("particle_bridge.BuildParticleGraphFromPositions on 70 particles: edge_index (2, E) and edge_features (E, 4) in graph_bridge's convention", fontsize=10)
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    pyplot.close(fig)


def make_grid_sampling(out):
    """A nodal field sampled on a coarse and a fine grid, and the thin-axis squeeze."""
    Kratos = _Kratos()
    from KratosMultiphysics.PhysicsNeMoApplication.bridges import grid_bridge
    pyplot = _Matplotlib()
    points, triangles = _TriangleGrid(14, 9, extent=(1.5, 1.0), jitter=0.02, seed=2)
    model = Kratos.Model()
    model_part = _TriangleModelPart(model, points, triangles, variables=(Kratos.TEMPERATURE,))
    for node in model_part.Nodes:
        node.SetSolutionStepValue(Kratos.TEMPERATURE, numpy.sin(2.2 * node.X) * numpy.cos(3.1 * node.Y) + 0.4 * node.X)
    values = numpy.array([node.GetSolutionStepValue(Kratos.TEMPERATURE) for node in model_part.Nodes])
    fig, axes = pyplot.subplots(1, 3, figsize=(14, 4.2))
    tri = axes[0].tripcolor(points[:, 0], points[:, 1], triangles, values, shading="gouraud", cmap="viridis")
    axes[0].triplot(points[:, 0], points[:, 1], triangles, color="white", lw=0.3, alpha=0.6)
    axes[0].set_title("TEMPERATURE on the Kratos nodes\n(unstructured triangles)", color=KRATOS_BLUE)
    for ax, shape in zip(axes[1:], [(8, 12, 2), (24, 36, 2)]):
        grid, box = grid_bridge.SampleFieldsOnGrid(model_part, [("TEMPERATURE", "node_historical")], shape)
        grid = numpy.asarray(grid)
        planar = grid[0].mean(axis=-1)  # the thin-axis squeeze is a mean over the axis of size 2
        low, high = numpy.asarray(box[0]), numpy.asarray(box[1])
        ax.imshow(planar.T, origin="lower", extent=[low[0], high[0], low[1], high[1]], cmap="viridis", vmin=values.min(), vmax=values.max(), interpolation="nearest")
        ax.set_title(f"SampleFieldsOnGrid, grid_shape {list(shape)}\n(C, D, H, W) = {tuple(grid.shape)}, thin axis squeezed", color=PNEMO_GREEN, fontsize=9.5)
    for ax in axes:
        ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    fig.colorbar(tri, ax=axes, fraction=0.015, pad=0.01)
    fig.savefig(out, bbox_inches="tight")
    pyplot.close(fig)


def make_calibration_views(out):
    """Synthetic ensemble: spread vs error, calibration curve, off-distribution growth."""
    pyplot = _Matplotlib()
    rng = numpy.random.default_rng(11)
    x = numpy.linspace(-2.0, 2.0, 400)
    truth = numpy.sin(2.0 * x) + 0.3 * x
    in_range = numpy.abs(x) <= 1.0
    members = []
    for k in range(6):
        bias = 0.25 * (k - 2.5) * numpy.clip(numpy.abs(x) - 1.0, 0, None) ** 2  # members drift apart off-range
        members.append(truth + 0.05 * rng.standard_normal(x.size) + bias + 0.02 * k * numpy.sin(5 * x))
    members = numpy.array(members)
    mean, std = members.mean(axis=0), members.std(axis=0, ddof=1)
    error = numpy.abs(mean - truth)
    fig, axes = pyplot.subplots(1, 3, figsize=(14, 4.2))
    ax = axes[0]
    ax.fill_between(x, mean - 2 * std, mean + 2 * std, color=PNEMO_GREEN, alpha=0.25, label="ensemble mean +/- 2 std")
    ax.plot(x, truth, color="k", lw=1.2, label="truth")
    ax.plot(x, mean, color=PNEMO_GREEN, lw=1.2, label="ensemble mean")
    ax.axvspan(-1, 1, color=KRATOS_BLUE, alpha=0.06); ax.text(0, ax.get_ylim()[0] + 0.1, "training range", ha="center", color=KRATOS_BLUE, fontsize=9)
    ax.set_title("spread against truth"); ax.legend(fontsize=8, loc="upper left"); ax.set_xlabel("input x")
    ax = axes[1]
    z = (truth - mean) / numpy.maximum(std, 1e-9)
    nominal = numpy.linspace(0.05, 0.99, 40)
    from math import erf, sqrt
    quantile = numpy.array([sqrt(2) * _InverseErf(p) for p in nominal])
    coverage_in = [numpy.mean(numpy.abs(z[in_range]) <= q) for q in quantile]
    overconfident = [numpy.mean(numpy.abs(z[in_range] * 2.0) <= q) for q in quantile]
    ax.plot(nominal, nominal, color="#999999", ls="--", label="perfect calibration")
    ax.plot(nominal, coverage_in, color=PNEMO_GREEN, lw=1.5, label="this ensemble, in range")
    ax.plot(nominal, overconfident, color=PINK, lw=1.5, label="bars half as wide as they should be")
    ax.set_xlabel("nominal coverage"); ax.set_ylabel("observed coverage"); ax.set_title("calibration: is the error bar the right size?"); ax.legend(fontsize=8, loc="upper left")
    ax = axes[2]
    ax.plot(x, error, color=PINK, lw=1.2, label="|error| of the mean")
    ax.plot(x, std, color=PNEMO_GREEN, lw=1.2, label="ensemble std")
    ax.axvspan(-1, 1, color=KRATOS_BLUE, alpha=0.06)
    ax.set_yscale("log"); ax.set_xlabel("input x"); ax.set_title("off-distribution: does the spread grow with the error?"); ax.legend(fontsize=8, loc="upper center")
    fig.suptitle("Three views of an error bar (synthetic six-member ensemble): the spread's direction is right off-range, its size is what calibration metrics check", fontsize=10)
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    pyplot.close(fig)


def _InverseErf(p):
    # Newton on erf for the coverage quantiles; scipy-free on purpose
    from math import erf, exp, pi, sqrt
    y = 0.0
    for _ in range(60):
        y -= (erf(y) - p) / (2.0 / sqrt(pi) * exp(-y * y))
    return y


def make_halo_partition(out):
    """A two-rank partition: owned nodes, ghosts, halo elements, and the truncated neighbourhood."""
    pyplot = _Matplotlib()
    points, triangles = _TriangleGrid(12, 6, extent=(2.0, 1.0), jitter=0.0)
    centroids = points[triangles].mean(axis=1)
    node_rank = (points[:, 0] > 1.0 + 1e-9).astype(int)
    element_rank = (centroids[:, 0] > 1.0).astype(int)
    owned0 = node_rank == 0
    elements0 = element_rank == 0
    ghost0 = numpy.zeros(len(points), bool)
    ghost0[numpy.unique(triangles[elements0])] = True
    ghost0 &= ~owned0
    touches_owned = numpy.any(owned0[triangles], axis=1)
    halo0 = touches_owned & ~elements0
    fig, axes = pyplot.subplots(1, 2, figsize=(14, 4.6))
    for ax, title, show_halo in [(axes[0], "what rank 0 holds: LocalMesh + GhostMesh (no element halo)", False),
                                 (axes[1], "graph_partition_utils.BuildHaloSubgraph: owned nodes + every element touching them", True)]:
        for t, tri in enumerate(triangles):
            poly = points[list(tri) + [tri[0]]]
            if elements0[t]:
                ax.fill(poly[:, 0], poly[:, 1], facecolor="#eaf2fb", edgecolor=KRATOS_BLUE, lw=0.6)
            elif show_halo and halo0[t]:
                ax.fill(poly[:, 0], poly[:, 1], facecolor="#fff4e0", edgecolor=ACCENT, lw=0.8, hatch="///")
            else:
                ax.fill(poly[:, 0], poly[:, 1], facecolor="none", edgecolor="#dddddd", lw=0.5)
        ax.scatter(points[owned0, 0], points[owned0, 1], s=18, color=KRATOS_BLUE, zorder=3, label="owned by rank 0")
        ax.scatter(points[ghost0, 0], points[ghost0, 1], s=22, color=ACCENT, zorder=3, marker="D", label="ghost (owned by rank 1)")
        if show_halo:
            halo_nodes = numpy.zeros(len(points), bool); halo_nodes[numpy.unique(triangles[halo0])] = True
            halo_nodes &= ~owned0 & ~ghost0
            ax.scatter(points[halo_nodes, 0], points[halo_nodes, 1], s=22, color=PINK, zorder=3, marker="^", label="halo node fetched from rank 1")
        ax.axvline(1.0, color="#999999", ls="--", lw=0.8)
        ax.set_title(title, fontsize=9.5); ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([]); ax.legend(fontsize=8, loc="lower left")
    # the interface node whose neighbourhood is truncated without the halo
    interface = numpy.flatnonzero(owned0 & (numpy.abs(points[:, 0] - 1.0) < 1e-9))[3]
    for ax in axes:
        ax.scatter(points[interface, 0], points[interface, 1], s=140, facecolors="none", edgecolors=PINK, lw=1.8, zorder=5)
    axes[0].annotate("an owned interface node:\nits right-hand neighbours live in\nelements rank 0 never sees -\none-hop aggregation is truncated here", xy=points[interface], xytext=(1.25, 0.72), fontsize=8.5, arrowprops=dict(arrowstyle="->", color=PINK), color=PINK)
    axes[1].annotate("with the halo of elements the\nneighbourhood matches the serial run\nat every owned node", xy=points[interface], xytext=(1.25, 0.72), fontsize=8.5, arrowprops=dict(arrowstyle="->", color=PINK), color=PINK)
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    pyplot.close(fig)


def make_benchmark_costs(out, quick=False):
    """Runs benchmark_bridges.py and plots microseconds per entity per path."""
    pyplot = _Matplotlib()
    divisions, grid, repeat = (4, 8, 1) if quick else (12, 16, 2)
    completed = subprocess.run([sys.executable, str(_APP_DIR / "benchmarks" / "benchmark_bridges.py"),
                                "--divisions", str(divisions), "--grid", str(grid), "--repeat", str(repeat)],
                               capture_output=True, text=True, timeout=1800, check=True)
    rows = []
    for line in completed.stdout.splitlines():
        match = re.match(r"^(\S.*?)\s{2,}(\d+)\s+([\d.]+)\s+([\d.]+)\s*$", line)
        if match and not line.startswith("case"):
            rows.append((match.group(1).strip(), int(match.group(2)), float(match.group(4))))
    if not rows:
        raise RuntimeError("no result rows parsed from the benchmark output:\n" + completed.stdout[-1500:])
    rows.sort(key=lambda r: r[2])
    names = [f"{name}  (n = {entities})" for name, entities, _ in rows]
    costs = numpy.array([cost for _, _, cost in rows])
    fig, ax = pyplot.subplots(figsize=(11, 0.42 * len(rows) + 1.6))
    colors = [PNEMO_GREEN if c < 1 else ACCENT if c < 10 else PINK for c in costs]
    ax.barh(names, costs, color=colors)
    for y, c in enumerate(costs):
        ax.text(c * 1.08, y, f"{c:.2f}", va="center", fontsize=8.5)
    ax.set_xscale("log"); ax.set_xlabel("microseconds per entity (log scale)")
    ax.set_title(f"benchmark_bridges.py --divisions {divisions} --grid {grid} --repeat {repeat}  on this machine; green < 1 us, orange < 10 us, pink above", fontsize=9.5)
    ax.grid(axis="x", which="both", alpha=0.3)
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    pyplot.close(fig)


def RasterizeSvgs(root, width=1200):
    """Writes a PNG twin next to every SVG diagram under root (needs cairosvg).

    The pages reference the SVG (crisp at any zoom, renders on the site and on
    GitHub); the PNG is kept alongside for viewers that cannot render SVG.
    Both are committed; regenerate the PNGs here after editing an SVG.
    """
    import cairosvg
    written = []
    for svg in sorted(Path(root).rglob("images/*.svg")):
        png = svg.with_suffix(".png")
        cairosvg.svg2png(url=str(svg), write_to=str(png), output_width=width)
        written.append(png)
        print(f"[ok]   {svg.stem:20s} -> {png}")
    return written


FIGURES = {
    "tessellation_modes": ("Mesh_Bridge/images/tessellation_modes.png", make_tessellation_modes, ("kratos",)),
    "higher_order_modes": ("Mesh_Bridge/images/higher_order_modes.png", make_higher_order_modes, ("kratos",)),
    "mesh_graph": ("Graph_Neural_Networks/images/mesh_graph.png", make_mesh_graph, ("kratos",)),
    "proximity_graph": ("Particle_Methods/images/proximity_graph.png", make_proximity_graph, ("kratos",)),
    "grid_sampling": ("Super_Resolution/images/grid_sampling.png", make_grid_sampling, ("kratos",)),
    "calibration_views": ("Uncertainty/images/calibration_views.png", make_calibration_views, ()),
    "halo_partition": ("Distributed/images/halo_partition.png", make_halo_partition, ()),
    "benchmark_costs": ("General/images/benchmark_costs.png", make_benchmark_costs, ("kratos", "benchmark")),
}


def Generate(names, out_dir, quick=False):
    """Generates the named figures into out_dir. Returns {name: path or exception}."""
    results = {}
    for name in names:
        relative, function, _ = FIGURES[name]
        target = Path(out_dir) / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        try:
            if name == "benchmark_costs":
                function(target, quick=quick)
            else:
                function(target)
            results[name] = target
            print(f"[ok]   {name:20s} -> {target}")
        except Exception as error:  # keep going: figures are independent
            results[name] = error
            print(f"[FAIL] {name:20s} {type(error).__name__}: {error}")
    return results


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out", default=str(_DOCS_DIR), help="root of the documentation pages (default: the docs tree)")
    parser.add_argument("--only", default="", help="comma-separated figure names")
    parser.add_argument("--list", action="store_true", help="list figure names and exit")
    parser.add_argument("--quick", action="store_true", help="smaller, faster benchmark run")
    parser.add_argument("--rasterize", action="store_true", help="also write a PNG twin next to every SVG diagram (needs cairosvg)")
    args = parser.parse_args(argv)
    if args.list:
        for name, (relative, _, needs) in FIGURES.items():
            print(f"{name:20s} {relative:55s} needs: {', '.join(needs) or 'numpy + matplotlib'}")
        return 0
    names = [n.strip() for n in args.only.split(",") if n.strip()] or list(FIGURES)
    unknown = [n for n in names if n not in FIGURES]
    if unknown:
        parser.error(f"unknown figure(s) {unknown}; see --list")
    if args.rasterize:
        RasterizeSvgs(args.out)
        if not args.only:
            names = list(FIGURES)
    results = Generate(names, args.out, quick=args.quick)
    failures = [n for n, r in results.items() if isinstance(r, Exception)]
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
