"""Delegation to physicsnemo-curator: Kratos as an ETL source.

Bridges Kratos model parts into ``physicsnemo_curator`` (the
``physicsnemo-curator`` package), whose pipelines are built as
``Source -> Filter -> Sink``:

- **Source**: ``CreateKratosMeshSource`` yields one
  ``physicsnemo.mesh.Mesh`` per item, tessellated from a Kratos model part
  through the mesh bridge, so any Kratos solve becomes a curator data
  source. Only the source side is ours - the sinks ship upstream.
- **Sinks**: ``CreateZarrSink`` / ``CreateVtuSink`` wrap curator's
  ``MeshZarrSink`` / ``MeshVTUSink``, which write AI-ready Zarr stores and
  VTK unstructured grids. Note they are typed on ``physicsnemo.mesh.Mesh``,
  NOT ``DomainMesh`` (only upstream's ``MeshSink`` accepts one), which is
  why the source yields ``BuildMesh(...)``; and the VTK sink writes .vtu,
  not .vtp.
- **Pipeline**: ``RunCuratorPipeline`` composes source/filters/sink and
  drives every item, returning the pipeline summary.

physicsnemo-curator is an optional runtime dependency, imported lazily with
actionable error messages. It is published on no package index, so it must
be installed from source, and its maturin build downloads a large Rust
toolchain - see the Mesh Bridge documentation page.
"""

import KratosMultiphysics as Kratos
from KratosMultiphysics.PhysicsNeMoApplication.bridges.mesh_bridge import domain_mesh_builder

_source_class = None

_INSTALL_HINT = (
    "Install it from source with e.g. 'pip install "
    "git+https://github.com/NVIDIA/physicsnemo-curator' (it is not published to PyPI "
    "under any name, and its maturin build downloads a large Rust toolchain).")


def _TryImportCurator():
    try:
        from physicsnemo_curator.core.base import Param, Source
        return Source, Param
    except ImportError as e:
        raise ImportError(
            "PhysicsNeMoApplication.curator_bridge requires physicsnemo-curator, which could "
            "not be imported. " + _INSTALL_HINT) from e


def _TryImportCuratorMeshSinks():
    """Curator's mesh sinks, which need the package's "mesh" extra.

    Kept separate from _TryImportCurator because the two fail
    independently: the source side needs only physicsnemo_curator.core.base,
    while importing anything under domains.mesh pulls the filter chain and
    therefore pyarrow, which the core install does not bring in.
    """
    try:
        from physicsnemo_curator.domains.mesh.sinks.mesh_vtu import MeshVTUSink
        from physicsnemo_curator.domains.mesh.sinks.mesh_zarr import MeshZarrSink
        return MeshZarrSink, MeshVTUSink
    except ImportError as e:
        raise ImportError(
            "PhysicsNeMoApplication.curator_bridge requires physicsnemo-curator's mesh sinks, "
            "which could not be imported (they additionally need the \"mesh\" extra - notably "
            "'pip install pyarrow'). " + _INSTALL_HINT) from e


def _NormalizeItems(model_parts):
    """Accepts a sequence of model parts or a (callable, count) pair.

    The callable form lets a trajectory be produced lazily - the item is
    built only when the pipeline reaches it, so a long series never has to
    exist in memory at once.
    """
    if isinstance(model_parts, tuple) and len(model_parts) == 2 and callable(model_parts[0]):
        factory, count = model_parts
        if count < 0:
            raise ValueError(f"Item count must be >= 0 [ count = {count} ].")
        return factory, int(count)
    items = list(model_parts)
    return (lambda index: items[index]), len(items)


def _GetKratosMeshSourceClass():
    global _source_class
    if _source_class is not None:
        return _source_class

    Source, Param = _TryImportCurator()

    class KratosMeshSource(Source["Mesh"]):
        """Curator source yielding tessellated Kratos model parts.

        Each index yields exactly one physicsnemo.mesh.Mesh built by the
        mesh bridge, with the requested Kratos fields attached as point/cell
        data. All tessellation knobs of BuildMesh apply, so higher-order and
        curved meshes reach curator sinks with the same fidelity they reach
        physicsnemo training.
        """

        name = "Kratos Mesh"
        description = (
            "Tessellate Kratos model parts into physicsnemo Mesh objects via the mesh bridge")

        @classmethod
        def params(cls):
            return [
                Param(name="source_container",
                      description="Kratos container to tessellate",
                      type=str, default="Elements", choices=["Elements", "Conditions"]),
                Param(name="tessellation_mode",
                      description="Quadrilateral/hexahedron splitting rule",
                      type=str, default="smallest_id_diagonal",
                      choices=["smallest_id_diagonal", "fan"]),
                Param(name="higher_order_mode",
                      description="Treatment of higher-order geometries",
                      type=str, default="reduce",
                      choices=["reduce", "subdivide", "curved"]),
                Param(name="curved_refinement_levels",
                      description="Refinement levels used by the curved mode",
                      type=int, default=2),
            ]

        def __init__(self, model_parts, field_specs=(), source_container="Elements",
                     tessellation_mode="smallest_id_diagonal", higher_order_mode="reduce",
                     curved_refinement_levels=2):
            self._item, self._count = _NormalizeItems(model_parts)
            self._field_specs = list(field_specs)
            self._source_container = source_container
            self._tessellation_mode = tessellation_mode
            self._higher_order_mode = higher_order_mode
            self._curved_refinement_levels = curved_refinement_levels
            # One cache per source: the tessellation is reused while the
            # mesh is unchanged, the field data is always collected fresh.
            self._provenance_cache = domain_mesh_builder.ProvenanceCache(
                source_container, tessellation_mode, higher_order_mode,
                curved_refinement_levels)

        def __len__(self) -> int:
            return self._count

        def __getitem__(self, index):
            if not 0 <= index < self._count:
                raise IndexError(f"Source index {index} out of range [ len = {self._count} ].")
            model_part = self._item(index)
            mesh, _ = domain_mesh_builder.BuildMesh(
                model_part, self._field_specs, self._source_container,
                self._tessellation_mode, self._higher_order_mode,
                self._curved_refinement_levels,
                provenance=self._provenance_cache.Get(model_part))
            yield mesh

    _source_class = KratosMeshSource
    return _source_class


def CreateKratosMeshSource(model_parts, field_specs=(), source_container="Elements",
                           tessellation_mode="smallest_id_diagonal",
                           higher_order_mode="reduce", curved_refinement_levels=2):
    """Creates a curator Source over Kratos model parts (lazily importing curator).

    Args:
        model_parts: A sequence of model parts, or a (callable, count) pair
            where callable(index) returns the model part for that index
            (lazy - nothing is built until the pipeline asks for it).
        field_specs: Iterable of (variable, data_location) pairs, as
            BuildMesh takes them.
        source_container / tessellation_mode / higher_order_mode /
            curved_refinement_levels: see BuildProvenance.

    Returns:
        A Source instance ready for .filter(...)/.write(...) or
        RunCuratorPipeline.
    """
    return _GetKratosMeshSourceClass()(
        model_parts, field_specs, source_container, tessellation_mode,
        higher_order_mode, curved_refinement_levels)


def CreateZarrSink(output_dir, compression_level: int = 3, chunk_size_mb: float = 1.0,
                   naming_template=None):
    """Curator's MeshZarrSink writing one .zarr store per item.

    The store carries the mesh points as "mesh_pos" with a leading time
    axis plus every attached field. Upstream also computes edge
    connectivity from global_data["edges"]; a bridge-built mesh carries
    none, so the sink logs that it found no edges and writes without them -
    run curator's EdgeComputeFilter first if the consumer needs them.
    """
    MeshZarrSink, _ = _TryImportCuratorMeshSinks()
    return MeshZarrSink(output_dir=str(output_dir), compression_level=compression_level,
                        chunk_size_mb=chunk_size_mb, naming_template=naming_template)


def CreateVtuSink(output_dir, naming_template=None):
    """Curator's MeshVTUSink writing one .vtu file per item (VTU, not VTP)."""
    _, MeshVTUSink = _TryImportCuratorMeshSinks()
    return MeshVTUSink(output_dir=str(output_dir), naming_template=naming_template)


def WriteMeshToCuratorSink(sink, mesh, index: int) -> list:
    """Writes one mesh straight to a curator sink, outside a pipeline.

    The sinks name their output from the index, so an in-loop exporter
    passes the Kratos step and gets one store/file per exported step.

    Returns:
        The list of paths the sink wrote.
    """
    return sink(iter([mesh]), index)


def RunCuratorPipeline(source, sink, filters=()) -> dict:
    """Builds and drives a curator pipeline over every source item.

    Args:
        source: A Source (e.g. from CreateKratosMeshSource).
        sink: A Sink (e.g. from CreateZarrSink / CreateVtuSink).
        filters: Optional iterable of curator Filters applied in order.

    Returns:
        The pipeline's summary dict.
    """
    pipeline = source
    for item in filters:
        pipeline = pipeline.filter(item)
    pipeline = pipeline.write(sink)
    # Pipeline is indexable but not iterable: processing is per-item and lazy.
    for index in range(len(pipeline)):
        pipeline[index]
    return pipeline.summary()
