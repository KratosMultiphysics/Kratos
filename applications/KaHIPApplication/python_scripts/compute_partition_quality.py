#!/usr/bin/env python3
"""
Partition Quality Comparison — KaHIPApplication vs MetisApplication
====================================================================

Reads a Kratos .mdpa mesh file, partitions it with the six KaHIP
preconfigurations (fast / eco / strong and their social-graph variants) and,
when MetisApplication is installed, also with METIS.  Three quality metrics
are measured for each strategy and written to a multi-panel matplotlib figure.

Metrics
-------
edge_cut    lower is better — edges in the nodal-adjacency graph that cross
            partition boundaries (KaHIP's primary optimisation objective).
imbalance   lower is better, 1.0 = perfect balance — max block size divided
            by mean block size.
time_s      lower is better — wall-clock seconds for the Execute() call.

Usage
-----
python compute_partition_quality.py \\
    [--mesh PATH_TO_MDPA]   (default: test_examples/cube.mdpa)
    [--partitions K]        (default: 4)
    [--output FILE.png]     (default: <mesh_dir>/<mesh_stem>_partition_quality.png)
    [--imbalance FLOAT]     (default: 0.03)
    [--seed INT]            (default: 0)
    [--num-trials INT]      (default: 1)
    [--include-social]      (also benchmark social-graph preconfigurations)
    [--dimension INT]       (spatial dimension for METIS; default: 3)
"""

import argparse
import collections
import pathlib
import shutil
import tempfile
import time

import KratosMultiphysics as KM

# Suppress Kratos banner and INFO messages so the comparison table is readable
KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.WARNING)

# ── optional dependencies ─────────────────────────────────────────────────────
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    HAS_PLOT = True
except ImportError:
    HAS_PLOT = False
    print("[WARN] matplotlib/numpy not available — text-only output will be produced.")

try:
    from KratosMultiphysics.KaHIPApplication import KaHIPDivideHeterogeneousInputProcess
    HAS_KAHIP = True
except ImportError:
    HAS_KAHIP = False
    print("[ERROR] KaHIPApplication not available.")

try:
    from KratosMultiphysics.MetisApplication import MetisDivideHeterogeneousInputProcess
    HAS_METIS = True
except ImportError:
    HAS_METIS = False

# ── default mesh ──────────────────────────────────────────────────────────────
_THIS_DIR = pathlib.Path(__file__).parent.resolve()
_DEFAULT_MESH = _THIS_DIR.parent / "test_examples" / "cube.mdpa"

# ── strategy display metadata ─────────────────────────────────────────────────
_PALETTE = {
    "fast":         "#4C72B0",
    "eco":          "#55A868",
    "strong":       "#C44E52",
    "fastsocial":   "#8172B2",
    "ecosocial":    "#CCB974",
    "strongsocial": "#64B5CD",
    "metis":        "#888888",
}

_LABELS = {
    "fast":         "KaHIP Fast",
    "eco":          "KaHIP Eco",
    "strong":       "KaHIP Strong",
    "fastsocial":   "KaHIP FastSocial",
    "ecosocial":    "KaHIP EcoSocial",
    "strongsocial": "KaHIP StrongSocial",
    "metis":        "METIS",
}


# ─────────────────────────────────────────────────────────────────────────────
# Mesh helpers
# ─────────────────────────────────────────────────────────────────────────────

def _read_model_part(mdpa_path: pathlib.Path):
    """Read an mdpa file; return (model, model_part)."""
    stem = str(mdpa_path.with_suffix(""))
    model = KM.Model()
    mp = model.CreateModelPart("mesh")
    mp.AddNodalSolutionStepVariable(KM.PARTITION_INDEX)
    flags = (KM.ModelPartIO.READ
             | KM.ModelPartIO.SKIP_TIMER
             | KM.ModelPartIO.IGNORE_VARIABLES_ERROR)
    KM.ModelPartIO(stem, flags).ReadModelPart(mp)
    return model, mp


def _build_nodal_adjacency(mp: KM.ModelPart) -> dict:
    """
    Element-based nodal adjacency: u and v are neighbours when they share an
    element.  Returns {node_id: frozenset(neighbour_ids)}.
    """
    adj: dict[int, set] = collections.defaultdict(set)
    for elem in mp.Elements:
        ids = [n.Id for n in elem.GetNodes()]
        for i, u in enumerate(ids):
            for v in ids[i + 1:]:
                adj[u].add(v)
                adj[v].add(u)
    return {k: frozenset(v) for k, v in adj.items()}


# ─────────────────────────────────────────────────────────────────────────────
# Quality metrics
# ─────────────────────────────────────────────────────────────────────────────

def _edge_cut(part_map: dict, adj: dict) -> int:
    """Count edges (u < v) where the two endpoints are in different blocks."""
    cut = 0
    for u, neighbours in adj.items():
        p_u = part_map.get(u, -1)
        for v in neighbours:
            if v > u and part_map.get(v, -2) != p_u:
                cut += 1
    return cut


def _imbalance(part_sizes: list) -> float:
    """max(block) / mean(block) among non-empty blocks; 1.0 = perfect."""
    nonempty = [s for s in part_sizes if s > 0]
    if not nonempty:
        return 0.0
    return max(nonempty) / (sum(nonempty) / len(nonempty))


# ─────────────────────────────────────────────────────────────────────────────
# Partitioner runner
# ─────────────────────────────────────────────────────────────────────────────

def _run_and_measure(mdpa_path: pathlib.Path, n_parts: int,
                     make_proc, adj: dict) -> dict:
    """
    Execute a partitioner (callable io → process) in an isolated temp directory,
    then read the resulting per-partition files to compute quality metrics.

    Returns a dict: edge_cut, imbalance, time_s, part_sizes.
    """
    mesh_name = mdpa_path.stem
    flags = (KM.ModelPartIO.READ
             | KM.ModelPartIO.SKIP_TIMER
             | KM.ModelPartIO.IGNORE_VARIABLES_ERROR)

    with tempfile.TemporaryDirectory(prefix="kratos_part_") as tmpdir:
        tmp = pathlib.Path(tmpdir)

        # The partitioner writes output relative to the IO stem, so copy the
        # mesh into the temp directory to keep all output there.
        shutil.copy(mdpa_path, tmp / mdpa_path.name)

        io = KM.ModelPartIO(str(tmp / mesh_name), flags)

        t0 = time.perf_counter()
        proc = make_proc(io)
        proc.Execute()
        elapsed = time.perf_counter() - t0

        # Collect node → partition mapping from the per-partition output files.
        # Every partition file marks local (owner) nodes with
        # PARTITION_INDEX == file_index; ghost nodes also appear but with their
        # true owner index.  We take the first occurrence for each node ID —
        # all files agree on the value.
        part_map: dict[int, int] = {}
        part_sizes = [0] * n_parts

        part_dir = tmp / f"{mesh_name}_partitioned"
        for i in range(n_parts):
            part_stem = str(part_dir / f"{mesh_name}_{i}")
            part_model = KM.Model()
            part_mp = part_model.CreateModelPart(f"p{i}")
            part_mp.AddNodalSolutionStepVariable(KM.PARTITION_INDEX)
            try:
                KM.ModelPartIO(part_stem, flags).ReadModelPart(part_mp)
            except Exception:
                continue

            for node in part_mp.Nodes:
                nid = node.Id
                p = int(node.GetSolutionStepValue(KM.PARTITION_INDEX))
                if nid not in part_map:
                    part_map[nid] = p
                if p == i:          # count only locally-owned nodes
                    part_sizes[i] += 1

    return {
        "edge_cut":   _edge_cut(part_map, adj),
        "imbalance":  _imbalance(part_sizes),
        "time_s":     elapsed,
        "part_sizes": part_sizes,
    }


# ─────────────────────────────────────────────────────────────────────────────
# Factory functions
# ─────────────────────────────────────────────────────────────────────────────

def _kahip_factory(mode: str, n_parts: int, imbalance: float,
                   seed: int, num_trials: int):
    """Return a callable (io) → KaHIPDivideHeterogeneousInputProcess."""
    settings = KM.Parameters(f"""{{
        "preconfiguration":      "{mode}",
        "imbalance":             {imbalance},
        "seed":                  {seed},
        "suppress_output":       true,
        "num_trials":            {num_trials},
        "verbosity":             0,
        "synchronize_conditions": false
    }}""")
    return lambda io: KaHIPDivideHeterogeneousInputProcess(io, n_parts, settings)


def _metis_factory(n_parts: int, dimension: int):
    """Return a callable (io) → MetisDivideHeterogeneousInputProcess."""
    return lambda io: MetisDivideHeterogeneousInputProcess(
        io, n_parts, dimension, 0, False)


# ─────────────────────────────────────────────────────────────────────────────
# Plotting
# ─────────────────────────────────────────────────────────────────────────────

_METRICS = [
    ("edge_cut",   "Edge cut  (lower = better)",              "Cut edges"),
    ("imbalance",  "Load imbalance  (lower = better, 1 = perfect)",
                                                              "max / mean block size"),
    ("time_s",     "Partition time  (lower = better)",        "Seconds"),
]


def _plot(results: dict, n_parts: int, mesh_name: str,
          output_path: pathlib.Path) -> None:
    strategies = list(results.keys())
    x = np.arange(len(strategies))
    colors = [_PALETTE.get(s, "#AAAAAA") for s in strategies]
    tick_labels = [_LABELS.get(s, s) for s in strategies]

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    fig.suptitle(
        f"Mesh partition quality — {mesh_name}  ({n_parts} partitions)",
        fontsize=13, fontweight="bold", y=1.03)

    for ax, (key, title, ylabel) in zip(axes, _METRICS):
        values = [results[s][key] for s in strategies]
        max_val = max(values) if values else 1.0
        bars = ax.bar(x, values, width=0.55, color=colors,
                      edgecolor="white", linewidth=0.8)

        # Annotate each bar with its numeric value
        for bar, val in zip(bars, values):
            label = (f"{val:,}" if key == "edge_cut"
                     else f"{val:.4f}" if key == "imbalance"
                     else f"{val:.2f}s")
            ax.text(bar.get_x() + bar.get_width() / 2,
                    bar.get_height() + max_val * 0.015,
                    label, ha="center", va="bottom", fontsize=7.5)

        ax.set_title(title, fontsize=9.5, pad=8)
        ax.set_ylabel(ylabel, fontsize=9)
        ax.set_xticks(x)
        ax.set_xticklabels(tick_labels, rotation=38, ha="right", fontsize=8)
        ax.yaxis.grid(True, linestyle="--", alpha=0.45)
        ax.set_axisbelow(True)
        # Give a little headroom above the tallest bar for the annotation
        ax.set_ylim(0, max_val * 1.18)

    # Legend
    handles = [
        plt.Rectangle((0, 0), 1, 1, color=_PALETTE.get(s, "#AAAAAA"),
                       label=_LABELS.get(s, s))
        for s in strategies
    ]
    fig.legend(handles=handles, loc="lower center",
               ncol=len(strategies), bbox_to_anchor=(0.5, -0.08),
               fontsize=8, frameon=False)

    plt.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    print(f"[OK] Plot saved → {output_path}")


# ─────────────────────────────────────────────────────────────────────────────
# Entry point
# ─────────────────────────────────────────────────────────────────────────────

def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Compare mesh partition quality: KaHIP strategies vs METIS.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--mesh", type=pathlib.Path, default=_DEFAULT_MESH,
                   metavar="PATH",
                   help="Path to the .mdpa mesh file")
    p.add_argument("--partitions", type=int, default=4, metavar="K",
                   help="Number of partitions")
    p.add_argument("--output", type=pathlib.Path,
                   default=None,
                   metavar="FILE",
                   help="Output PNG filename (default: <mesh_dir>/<mesh_stem>_partition_quality.png)")
    p.add_argument("--imbalance", type=float, default=0.03,
                   help="Allowed imbalance fraction passed to KaHIP (e.g. 0.03 = 3%%)")
    p.add_argument("--seed", type=int, default=0,
                   help="Random seed")
    p.add_argument("--num-trials", type=int, default=1, dest="num_trials",
                   help="KaHIP num_trials: run this many times, keep best result")
    p.add_argument("--include-social", action="store_true",
                   help="Also benchmark social-graph preconfigurations "
                        "(fastsocial / ecosocial / strongsocial)")
    p.add_argument("--dimension", type=int, default=3, choices=[2, 3],
                   help="Spatial dimension forwarded to METIS")
    return p.parse_args()


def main() -> None:
    args = _parse_args()

    mdpa_path = args.mesh.resolve()
    if not mdpa_path.exists():
        raise FileNotFoundError(f"Mesh file not found: {mdpa_path}")

    # Default output: alongside the mesh file
    output_path = (args.output.resolve() if args.output is not None
                   else mdpa_path.parent / f"{mdpa_path.stem}_partition_quality.png")

    n_parts = args.partitions
    if n_parts < 2:
        raise ValueError("--partitions must be >= 2")

    print("=" * 60)
    print("  KaHIP Partition Quality Comparison")
    print("=" * 60)
    print(f"  Mesh       : {mdpa_path}")
    print(f"  Partitions : {n_parts}")
    print(f"  Imbalance  : {args.imbalance}")
    print(f"  Seed       : {args.seed}")
    print(f"  Trials     : {args.num_trials}")
    print("=" * 60)
    print()

    # ── read mesh once and build nodal adjacency ──────────────────────────────
    print("Reading mesh …")
    _, mp = _read_model_part(mdpa_path)
    adj = _build_nodal_adjacency(mp)
    n_edges = sum(len(v) for v in adj.values()) // 2
    print(f"  {mp.NumberOfNodes():,} nodes  |  "
          f"{mp.NumberOfElements():,} elements  |  "
          f"{n_edges:,} nodal edges")
    print()

    # ── build strategy list ───────────────────────────────────────────────────
    strategies: dict[str, callable] = {}

    if HAS_KAHIP:
        modes = ["fast", "eco", "strong"]
        if args.include_social:
            modes += ["fastsocial", "ecosocial", "strongsocial"]
        for mode in modes:
            strategies[mode] = _kahip_factory(
                mode, n_parts, args.imbalance, args.seed, args.num_trials)
    else:
        print("[ERROR] KaHIPApplication not found — no KaHIP strategies available.\n")

    if HAS_METIS:
        strategies["metis"] = _metis_factory(n_parts, args.dimension)
        print("[INFO] MetisApplication found — METIS will be included in comparison.")
    else:
        print("[INFO] MetisApplication not found — METIS comparison skipped.")
    print()

    if not strategies:
        raise RuntimeError("No partitioners available. "
                           "Ensure KaHIPApplication (and optionally MetisApplication) "
                           "is built and importable.")

    # ── run each strategy ─────────────────────────────────────────────────────
    col = max(len(s) for s in strategies) + 2
    header = (f"{'Strategy':<{col}}  {'Edge cut':>12}  "
              f"{'Imbalance':>12}  {'Time (s)':>10}")
    print(header)
    print("-" * len(header))

    results: dict[str, dict] = {}

    for name, factory in strategies.items():
        label = _LABELS.get(name, name)
        print(f"  Running {label} …", end="", flush=True)
        try:
            m = _run_and_measure(mdpa_path, n_parts, factory, adj)
        except Exception as exc:
            print(f"\r[FAIL] {label}: {exc}")
            continue

        results[name] = m
        print(f"\r{name:<{col}}  {m['edge_cut']:>12,}  "
              f"{m['imbalance']:>12.4f}  {m['time_s']:>10.3f}")

    print()

    if not results:
        print("[ERROR] All strategies failed — nothing to plot.")
        return

    # ── print summary ─────────────────────────────────────────────────────────
    best_cut = min(results, key=lambda s: results[s]["edge_cut"])
    best_bal = min(results, key=lambda s: results[s]["imbalance"])
    best_spd = min(results, key=lambda s: results[s]["time_s"])
    print(f"  Best edge cut  : {_LABELS.get(best_cut, best_cut)}"
          f" ({results[best_cut]['edge_cut']:,} cut edges)")
    print(f"  Best balance   : {_LABELS.get(best_bal, best_bal)}"
          f" (imbalance = {results[best_bal]['imbalance']:.4f})")
    print(f"  Fastest        : {_LABELS.get(best_spd, best_spd)}"
          f" ({results[best_spd]['time_s']:.3f}s)")
    print()

    # ── plot ──────────────────────────────────────────────────────────────────
    if HAS_PLOT:
        _plot(results, n_parts, mdpa_path.name, output_path)
    else:
        print("[INFO] Install matplotlib + numpy to generate the comparison plot.")
        print("       pip install matplotlib numpy")


if __name__ == "__main__":
    main()
