"""Benchmarks the application's bridge hot paths on large meshes.

The roadmap kept a "C++ acceleration of the gather/scatter paths if
profiling on large meshes ever demands it" bullet; this script provides the
profiling. It times, on structured meshes of a configurable size, every
per-entity path a training/inference step exercises:

- tessellation + provenance construction (tet and hex meshes, both modes),
- nodal gather/scatter round trip through the provenance map,
- element-field scatter-back (ScatterFieldBack, mean reduction),
- grid sampling + grid scatter (SampleFieldsOnGrid / ScatterGridToNodes),
  with --sampling both comparing the Kratos locator against physicsnemo's BVH,
- ROM unknowns gather/scatter (rom_bridge, VariableUtils-backed),
- the serial path of distributed_utils.GatherFieldToRank0.

Run (after the usual PYTHONPATH/LD_LIBRARY_PATH exports):

    python3 benchmarks/benchmark_bridges.py --divisions 32 --grid 64

Plain wall-clock timing (median of --repeat runs), printed as a table with
per-entity costs. Pure Kratos + numpy - no torch/physicsnemo needed.
"""

import argparse
import statistics
import time

import numpy

import KratosMultiphysics as Kratos
from KratosMultiphysics.PhysicsNeMoApplication.bridges import graph_bridge
from KratosMultiphysics.PhysicsNeMoApplication.bridges import grid_bridge
from KratosMultiphysics.PhysicsNeMoApplication.bridges import particle_bridge
from KratosMultiphysics.PhysicsNeMoApplication.bridges import rom_bridge
from KratosMultiphysics.PhysicsNeMoApplication.distributed import distributed_utils
from KratosMultiphysics.PhysicsNeMoApplication.bridges.mesh_bridge import domain_mesh_builder
from KratosMultiphysics.PhysicsNeMoApplication.utilities import array_backend_utils

_HEX_CORNERS = ((0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
                (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1))
_HEX_TO_TETS = ((0, 1, 2, 6), (0, 2, 3, 6), (0, 3, 7, 6),
                (0, 7, 4, 6), (0, 4, 5, 6), (0, 5, 1, 6))


def _CreateStructuredModelPart(model, name, divisions, element_type):
    """Structured unit-cube mesh of tets ("tet") or hexes ("hex")."""
    model_part = model.CreateModelPart(name)
    model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
    model_part.AddNodalSolutionStepVariable(Kratos.TEMPERATURE)
    properties = model_part.CreateNewProperties(1)

    n = divisions + 1
    for i in range(n):
        for j in range(n):
            for k in range(n):
                model_part.CreateNewNode(
                    i * n * n + j * n + k + 1, i / divisions, j / divisions, k / divisions)

    def node_id(i, j, k):
        return i * n * n + j * n + k + 1

    element_id = 0
    for i in range(divisions):
        for j in range(divisions):
            for k in range(divisions):
                corners = [node_id(i + di, j + dj, k + dk) for di, dj, dk in _HEX_CORNERS]
                if element_type == "hex":
                    element_id += 1
                    model_part.CreateNewElement("Element3D8N", element_id, corners, properties)
                else:
                    for tet in _HEX_TO_TETS:
                        element_id += 1
                        model_part.CreateNewElement(
                            "Element3D4N", element_id, [corners[c] for c in tet], properties)

    for node in model_part.Nodes:
        node.SetSolutionStepValue(Kratos.PRESSURE, 1.0 + node.X + 2.0 * node.Y - node.Z)
    return model_part


def _Time(function, repeat):
    timings = []
    for _ in range(repeat):
        start = time.perf_counter()
        result = function()
        timings.append(time.perf_counter() - start)
    return statistics.median(timings), result


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--divisions", type=int, default=32,
                        help="cells per axis of the unit-cube meshes (default 32)")
    parser.add_argument("--grid", type=int, default=64,
                        help="lattice points per axis for the grid-sampling case (default 64)")
    parser.add_argument("--repeat", type=int, default=3,
                        help="repetitions per case; the median is reported (default 3)")
    parser.add_argument("--sampling", choices=("kratos", "physicsnemo", "both"),
                        default="kratos",
                        help="which grid-sampling backend to time: the Kratos point "
                             "locator, physicsnemo's mesh BVH, or \"both\" with the "
                             "speedup (default kratos)")
    parser.add_argument("--backend", choices=("numpy", "cupy", "both"), default="numpy",
                        help="array backend for the paths that accept one; \"both\" times "
                             "each backend and reports the speedup (default numpy)")
    arguments = parser.parse_args()

    backends = ("numpy", "cupy") if arguments.backend == "both" else (arguments.backend,)
    if "cupy" in backends:
        if not array_backend_utils.IsCuPyAvailable():
            print("cupy was requested but no CUDA device answered; "
                  "the cupy column will repeat the numpy path.\n", flush=True)
        else:
            # Create the CUDA context and prime the allocator up front. The
            # first device touch in a process costs ~100 ms, and letting it
            # land inside a timed case understates that case by ~2x.
            import cupy
            cupy.asarray(numpy.zeros(1024)).sum()
            cupy.cuda.Stream.null.synchronize()

    model = Kratos.Model()
    print(f"building structured meshes (divisions = {arguments.divisions}) ...", flush=True)
    tet_part = _CreateStructuredModelPart(model, "TetPart", arguments.divisions, "tet")
    hex_part = _CreateStructuredModelPart(model, "HexPart", arguments.divisions, "hex")
    n_nodes = tet_part.NumberOfNodes()
    n_tets = tet_part.NumberOfElements()
    n_hexes = hex_part.NumberOfElements()
    print(f"  {n_nodes} nodes, {n_tets} tetrahedra / {n_hexes} hexahedra\n", flush=True)

    rows = []

    def record(name, entity_count, seconds):
        rows.append((name, entity_count, seconds, None))

    def record_per_backend(name, entity_count, call):
        """Times a backend-accepting path once per requested backend.

        The thresholds that keep small problems on the CPU are dropped for
        the duration, so a "cupy" column is really the device path rather
        than numpy wearing its name.
        """
        timings = {}
        saved = array_backend_utils.DEFAULT_SIZE_THRESHOLD
        saved_rom = rom_bridge._GPU_BASIS_THRESHOLD
        array_backend_utils.DEFAULT_SIZE_THRESHOLD = 0
        rom_bridge._GPU_BASIS_THRESHOLD = 0
        try:
            for backend in backends:
                # CuPy kernels are asynchronous: without the synchronize the
                # timer would measure the launch, not the work.
                if backend == "cupy" and array_backend_utils.IsCuPyAvailable():
                    import cupy

                    def timed(backend=backend):
                        result = call(backend)
                        cupy.cuda.Stream.null.synchronize()
                        return result
                else:
                    def timed(backend=backend):
                        return call(backend)
                timed()  # warm up: kernel load, basis upload
                timings[backend], _ = _Time(timed, arguments.repeat)
        finally:
            array_backend_utils.DEFAULT_SIZE_THRESHOLD = saved
            rom_bridge._GPU_BASIS_THRESHOLD = saved_rom
        rows.append((name, entity_count, timings[backends[0]],
                     timings.get("cupy") if len(backends) > 1 else None))

    # --- tessellation + provenance ------------------------------------------
    seconds, tet_provenance = _Time(lambda: domain_mesh_builder.BuildProvenance(tet_part),
                                    arguments.repeat)
    record("BuildProvenance tet (smallest_id)", n_tets, seconds)
    seconds, _ = _Time(lambda: domain_mesh_builder.BuildProvenance(hex_part), arguments.repeat)
    record("BuildProvenance hex (smallest_id)", n_hexes, seconds)
    seconds, _ = _Time(lambda: domain_mesh_builder.BuildProvenance(hex_part, tessellation_mode="fan"),
                       arguments.repeat)
    record("BuildProvenance hex (fan)", n_hexes, seconds)

    # --- nodal gather/scatter through the provenance map --------------------
    node_ids = [node.Id for node in tet_part.Nodes]
    nodal_field = numpy.random.default_rng(0).standard_normal(n_nodes)
    seconds, gathered = _Time(
        lambda: tet_provenance.GatherNodalField(node_ids, nodal_field), arguments.repeat)
    record("GatherNodalField", n_nodes, seconds)
    seconds, _ = _Time(
        lambda: tet_provenance.ScatterNodalField(node_ids, gathered), arguments.repeat)
    record("ScatterNodalField", n_nodes, seconds)

    # --- full nodal scatter-back (tensor adaptors + provenance) -------------
    seconds, _ = _Time(
        lambda: domain_mesh_builder.ScatterFieldBack(
            tet_provenance, gathered, tet_part, Kratos.TEMPERATURE, "node_historical"),
        arguments.repeat)
    record("ScatterFieldBack nodal", n_nodes, seconds)

    # --- element-field scatter-back (per-entity reduction path) -------------
    cell_field = numpy.random.default_rng(1).standard_normal(tet_provenance.number_of_cells)
    seconds, _ = _Time(
        lambda: domain_mesh_builder.ScatterFieldBack(
            tet_provenance, cell_field, tet_part, Kratos.TEMPERATURE, "element"),
        arguments.repeat)
    record("ScatterFieldBack element", n_tets, seconds)

    # --- graph edge features (recomputed every step on a deforming mesh) ----
    _, edge_index, _, _ = graph_bridge.BuildGraph(tet_part)
    record_per_backend(
        "ComputeEdgeFeatures", edge_index.shape[1],
        lambda backend: graph_bridge.ComputeEdgeFeatures(tet_part, edge_index, backend=backend))

    # --- grid sampling + scatter --------------------------------------------
    grid_shape = (arguments.grid,) * 3
    grid = bounding_box = None
    try:
        seconds, (grid, bounding_box) = _Time(
            lambda: grid_bridge.SampleFieldsOnGrid(
                tet_part, [("PRESSURE", "node_historical")], grid_shape),
            arguments.repeat)
        record(f"SampleFieldsOnGrid {arguments.grid}^3", int(numpy.prod(grid_shape)), seconds)

        if arguments.sampling in ("physicsnemo", "both"):
            # The question item 20 asks: is the mesh's own BVH faster than
            # the point locator? Measure it, do not assume it - the default
            # backend stays "kratos" until a number says otherwise.
            try:
                bvh_seconds, _ = _Time(
                    lambda: grid_bridge.SampleFieldsOnGrid(
                        tet_part, [("PRESSURE", "node_historical")], grid_shape,
                        bounding_box=bounding_box, backend="physicsnemo"),
                    arguments.repeat)
                record(f"SampleFieldsOnGrid {arguments.grid}^3 (BVH)",
                       int(numpy.prod(grid_shape)), bvh_seconds)
                ratio = seconds / bvh_seconds if bvh_seconds > 0 else float("inf")
                verdict = "the BVH wins" if ratio > 1.0 else "the locator wins"
                print(f"  grid sampling: locator {seconds * 1e3:.1f} ms vs BVH "
                      f"{bvh_seconds * 1e3:.1f} ms -> {verdict} "
                      f"({ratio:.2f}x)\n", flush=True)
            except ImportError as error:
                print(f"  skipping the BVH sampling comparison: {error}\n", flush=True)
    except TypeError as error:
        # The vectorized locator is a compiled-core entry point, and its
        # signature has drifted across core builds. Skipping the case keeps
        # every other measurement in this run rather than losing the lot.
        print(f"  skipping SampleFieldsOnGrid: {error.__class__.__name__} from "
              f"BinBasedFastPointLocator3D.VectorizedFind (core build mismatch)\n", flush=True)

    if grid is not None:
        record_per_backend(
            "ScatterGridToNodes", n_nodes,
            lambda backend: grid_bridge.ScatterGridToNodes(
                grid, bounding_box, tet_part, [("TEMPERATURE", "node_historical")],
                backend=backend))

    # --- ROM gather/scatter (VariableUtils-backed) --------------------------
    basis = rom_bridge.RomBasis(
        phi=numpy.zeros((n_nodes, 1)),
        node_ids=numpy.array(node_ids, dtype=numpy.int64),
        nodal_unknowns=("PRESSURE",),
        singular_values=None)
    seconds, unknowns = _Time(
        lambda: rom_bridge.GatherUnknownsVector(tet_part, basis), arguments.repeat)
    record("ROM GatherUnknownsVector", n_nodes, seconds)
    seconds, _ = _Time(
        lambda: rom_bridge.ScatterUnknownsVector(tet_part, basis, unknowns), arguments.repeat)
    record("ROM ScatterUnknownsVector", n_nodes, seconds)

    # --- ROM projection (the dense basis GEMV, where the GPU pays) ----------
    dense_basis = rom_bridge.RomBasis(
        phi=numpy.random.default_rng(2).standard_normal((n_nodes, 64)),
        node_ids=numpy.array(node_ids, dtype=numpy.int64),
        nodal_unknowns=("PRESSURE",),
        singular_values=None)
    snapshot = numpy.random.default_rng(3).standard_normal(n_nodes)
    record_per_backend(
        f"ROM ProjectToReducedSpace ({n_nodes}x64)", n_nodes,
        lambda backend: rom_bridge.ProjectToReducedSpace(dense_basis, snapshot, backend=backend))

    # --- particle proximity graph (the quadratic path) ----------------------
    particle_positions = numpy.random.default_rng(4).random((2000, 3))
    record_per_backend(
        "BuildParticleGraph radius (N=2000)", 2000,
        lambda backend: particle_bridge.BuildParticleGraphFromPositions(
            particle_positions,
            Kratos.Parameters('{"type" : "radius", "radius" : 0.05, "backend" : "%s"}' % backend)))

    # --- distributed gather, serial path ------------------------------------
    seconds, _ = _Time(
        lambda: distributed_utils.GatherFieldToRank0(tet_part, "PRESSURE", "node_historical"),
        arguments.repeat)
    record("GatherFieldToRank0 (serial)", n_nodes, seconds)

    comparing = len(backends) > 1
    header = f"{'case':40s} {'entities':>10s} {'seconds':>10s} {'us/entity':>10s}"
    if comparing:
        header += f" {'cupy s':>10s} {'speedup':>9s}"
    print(header)
    print("-" * (len(header) + 1))
    for name, entity_count, seconds, cupy_seconds in rows:
        line = (f"{name:40s} {entity_count:10d} {seconds:10.3f} "
                f"{1e6 * seconds / entity_count:10.2f}")
        if comparing:
            line += (f" {cupy_seconds:10.3f} {seconds / cupy_seconds:8.2f}x"
                     if cupy_seconds else f" {'-':>10s} {'numpy only':>9s}")
        print(line)


if __name__ == "__main__":
    main()
