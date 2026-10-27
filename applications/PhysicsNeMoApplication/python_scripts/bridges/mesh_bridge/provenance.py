"""Field-provenance map between a Kratos model part and its tessellation.

In the "reduce" and "subdivide" higher-order modes the tessellation only
uses existing Kratos nodes (see tessellation.py), so every simplex point
corresponds to exactly one node and nodal gather/scatter is an exact
bijection. The "curved" mode (curved_tessellation.py) additionally samples
the isoparametric geometry at SYNTHETIC points: those carry
``point_provenance == -1`` (always a trailing block after the real, id-
sorted points) plus parent/local-coordinate/shape-weight arrays -
GatherNodalField interpolates their values through the parent's shape
functions, while ScatterNodalField writes real nodes only (synthetic rows
are dropped; the round trip stays exact on real nodes). Cell fields map
many simplex sub-cells back to one source entity and need a reduction.

Pure Python + numpy: this module never imports torch or physicsnemo.
"""

from dataclasses import dataclass, field

import numpy


def _RowsOf(reference_ids, query_ids) -> numpy.ndarray:
    """query_ids -> row indices into reference_ids, vectorized.

    Stands in for the ``{id: row}`` dict plus ``numpy.fromiter`` generator
    this module used to build on every gather: both are interpreter-level
    loops over numpy scalars, and both gather paths run every step. Raises
    KeyError on an id absent from reference_ids, exactly as the dict lookup
    it replaces did.
    """
    reference = numpy.asarray(reference_ids, dtype=numpy.int64).ravel()
    query = numpy.asarray(query_ids, dtype=numpy.int64)
    if reference.size == 0:
        if query.size:
            raise KeyError(int(query.ravel()[0]))
        return numpy.empty(query.shape, dtype=numpy.int64)

    # Not assuming sortedness: Kratos containers are id-sorted in practice,
    # but the ids handed in here are an arbitrary caller-supplied order.
    order = numpy.argsort(reference, kind="stable")
    sorted_reference = reference[order]
    position = numpy.searchsorted(sorted_reference, query)
    clipped = numpy.minimum(position, reference.size - 1)
    missing = sorted_reference[clipped] != query
    if missing.any():
        raise KeyError(int(query[missing].ravel()[0]))
    return order[clipped]


@dataclass
class MeshProvenanceMap:
    """Flat-array provenance of a tessellated model part.

    Attributes:
        simplex_points: (P, 3) float64 point coordinates.
        simplex_cells: (C, k+1) int64 connectivity referencing simplex_points
            rows (k = 2 for triangles, 3 for tetrahedra).
        point_provenance: (P,) int64, the Kratos node id of each simplex
            point, or -1 for synthetic (curved-mode) points, which always
            form a trailing block.
        cell_provenance: (C, 2) int64, (source_entity_id, sub_cell_index) of
            each simplex cell.
        source_container: "Elements" or "Conditions".
        synthetic_parent_ids: (S,) int64 source entity id per synthetic
            point (None when the tessellation has no synthetic points).
        synthetic_local_coordinates: (S, 3) float64 local coordinates of
            each synthetic point in its parent geometry (padded to 3).
        synthetic_node_ids: (S, n_max) int64 parent node ids per synthetic
            point; rows are padded by repeating the row's first id.
        synthetic_weights: (S, n_max) float64 shape-function weights aligned
            with synthetic_node_ids; padding entries carry weight 0.
    """
    simplex_points: numpy.ndarray
    simplex_cells: numpy.ndarray
    point_provenance: numpy.ndarray
    cell_provenance: numpy.ndarray
    source_container: str
    synthetic_parent_ids: numpy.ndarray = None
    synthetic_local_coordinates: numpy.ndarray = None
    synthetic_node_ids: numpy.ndarray = None
    synthetic_weights: numpy.ndarray = None
    _node_id_to_point: dict = field(default_factory=dict, repr=False)
    _sorted_real_ids: numpy.ndarray = field(default=None, repr=False)
    _sorted_real_rows: numpy.ndarray = field(default=None, repr=False)

    def __post_init__(self):
        if not self._node_id_to_point:
            self._node_id_to_point = {
                int(node_id): index for index, node_id in enumerate(self.point_provenance)
                if node_id >= 0
            }
        if self._sorted_real_ids is None:
            # The id -> point-row map as sorted arrays, so the per-step
            # lookups below are a searchsorted rather than a Python loop.
            # Synthetic points (id -1) are excluded: they are not addressable
            # by node id at all.
            provenance = numpy.asarray(self.point_provenance, dtype=numpy.int64)
            real_rows = numpy.flatnonzero(provenance >= 0)
            real_ids = provenance[real_rows]
            order = numpy.argsort(real_ids, kind="stable")
            self._sorted_real_ids = real_ids[order]
            self._sorted_real_rows = real_rows[order]

    def _LookupPointRows(self, node_ids):
        """node ids -> (point rows, found mask), vectorized against the cache."""
        query = numpy.asarray(node_ids, dtype=numpy.int64).ravel()
        if self._sorted_real_ids.size == 0:
            return (numpy.zeros(query.shape, dtype=numpy.int64),
                    numpy.zeros(query.shape, dtype=bool))
        position = numpy.searchsorted(self._sorted_real_ids, query)
        clipped = numpy.minimum(position, self._sorted_real_ids.size - 1)
        found = self._sorted_real_ids[clipped] == query
        return self._sorted_real_rows[clipped], found

    @property
    def number_of_synthetic_points(self) -> int:
        return 0 if self.synthetic_parent_ids is None else len(self.synthetic_parent_ids)

    @property
    def number_of_points(self) -> int:
        return self.simplex_points.shape[0]

    @property
    def number_of_cells(self) -> int:
        return self.simplex_cells.shape[0]

    def GetPointIndices(self, node_ids) -> numpy.ndarray:
        """Maps Kratos node ids to simplex point row indices."""
        rows, found = self._LookupPointRows(node_ids)
        if not found.all():
            raise KeyError(int(numpy.asarray(node_ids, dtype=numpy.int64).ravel()[
                numpy.flatnonzero(~found)[0]]))
        return rows

    def GatherNodalField(self, node_ids, field_array) -> numpy.ndarray:
        """Reorders a per-node field (aligned with node_ids) onto simplex points.

        Real points take their node's value exactly; synthetic (curved-mode)
        points are interpolated through their parent's shape functions -
        exact for fields the parent geometry interpolates exactly. The
        result dtype is promoted to float when synthetic points exist.

        Args:
            node_ids: The node ids the rows of field_array correspond to (e.g.
                the ids of the model part's Nodes container, in order).
            field_array: (N, ...) array of per-node values.

        Returns:
            (P, ...) array aligned with simplex_points.
        """
        field_array = numpy.asarray(field_array)

        synthetic_count = self.number_of_synthetic_points
        if synthetic_count == 0:
            return field_array[_RowsOf(node_ids, self.point_provenance)]

        real_count = self.number_of_points - synthetic_count
        gather = _RowsOf(node_ids, self.point_provenance[:real_count])

        # synthetic block: weights @ parent-node values, fully vectorized
        weight_rows = _RowsOf(node_ids, self.synthetic_node_ids.ravel()).reshape(
            self.synthetic_node_ids.shape)                            # (S, n_max)
        contributions = field_array[weight_rows]                      # (S, n_max, ...)
        weights = self.synthetic_weights.reshape(
            self.synthetic_weights.shape + (1,) * (field_array.ndim - 1))
        synthetic_values = (contributions * weights).sum(axis=1)      # (S, ...)

        result = numpy.empty(
            (self.number_of_points,) + field_array.shape[1:],
            dtype=numpy.result_type(field_array.dtype, numpy.float64))
        result[:real_count] = field_array[gather]
        result[real_count:] = synthetic_values
        return result

    def ScatterNodalField(self, node_ids, point_field) -> numpy.ndarray:
        """Inverse of GatherNodalField: simplex-point field -> per-node field.

        Exact (no averaging): every corner node maps to exactly one simplex
        point. Nodes of the model part that are not part of any tessellated
        entity keep a zero value.

        Args:
            node_ids: The node ids the output rows should correspond to.
            point_field: (P, ...) array aligned with simplex_points.

        Returns:
            (N, ...) array aligned with node_ids.
        """
        point_field = numpy.asarray(point_field)
        result = numpy.zeros((len(node_ids),) + point_field.shape[1:], dtype=point_field.dtype)
        rows, found = self._LookupPointRows(node_ids)
        # Nodes belonging to no tessellated entity are simply not written,
        # so they keep the zero the result was allocated with.
        result[found] = point_field[rows[found]]
        return result

    def AggregateCellField(self, cell_field, reduction="mean", cell_weights=None):
        """Reduces a per-simplex-cell field to one value per source entity.

        Args:
            cell_field: (C, ...) array aligned with simplex_cells.
            reduction: "mean", "weighted_mean" (needs cell_weights, e.g.
                sub-cell volumes) or "first" (categorical labels).
            cell_weights: (C,) weights for "weighted_mean".

        Returns:
            (entity_ids, values): (E,) int64 source entity ids (sorted) and
            (E, ...) reduced values aligned with them.
        """
        cell_field = numpy.asarray(cell_field)
        source_ids = self.cell_provenance[:, 0]
        entity_ids, inverse = numpy.unique(source_ids, return_inverse=True)

        if reduction == "first":
            first_rows = numpy.zeros(len(entity_ids), dtype=numpy.int64)
            # reversed so the first occurrence wins
            first_rows[inverse[::-1]] = numpy.arange(len(source_ids) - 1, -1, -1)
            return entity_ids, cell_field[first_rows]

        extra_axes = (1,) * (cell_field.ndim - 1)
        if reduction == "mean":
            weights = numpy.ones(len(source_ids), dtype=numpy.float64)
        elif reduction == "weighted_mean":
            if cell_weights is None:
                raise ValueError("\"weighted_mean\" reduction requires cell_weights.")
            weights = numpy.asarray(cell_weights, dtype=numpy.float64)
        else:
            raise ValueError(f"Unsupported reduction \"{reduction}\". Use \"mean\", \"weighted_mean\" or \"first\".")

        # bincount per channel rather than numpy.add.at over the whole field:
        # add.at is numpy's slowest reduction path (an unbuffered ufunc loop)
        # and this runs per step. Channel counts are small (1, 3, 9), so the
        # loop is a handful of C calls - measured 1.3x on scalars and 5.4x on
        # vectors against add.at, where folding the channel into the bin index
        # instead loses on scalar fields.
        inverse = inverse.ravel()
        weight_sums = numpy.bincount(inverse, weights=weights, minlength=len(entity_ids))
        weighted = (cell_field * weights.reshape((-1,) + extra_axes)).reshape(
            len(source_ids), -1)
        if weighted.shape[1]:
            sums = numpy.stack(
                [numpy.bincount(inverse, weights=weighted[:, channel],
                                minlength=len(entity_ids))
                 for channel in range(weighted.shape[1])], axis=1)
        else:
            sums = numpy.zeros((len(entity_ids), 0), dtype=numpy.float64)
        sums = sums.reshape((len(entity_ids),) + cell_field.shape[1:])
        return entity_ids, sums / weight_sums.reshape((-1,) + extra_axes)

    def ComputeSimplexMeasures(self) -> numpy.ndarray:
        """Areas (triangles) or volumes (tetrahedra) of all simplex cells."""
        pts = self.simplex_points[self.simplex_cells]
        if self.simplex_cells.shape[1] == 3:
            cross = numpy.cross(pts[:, 1] - pts[:, 0], pts[:, 2] - pts[:, 0])
            return 0.5 * numpy.linalg.norm(cross, axis=-1)
        if self.simplex_cells.shape[1] == 4:
            v = pts[:, 1:] - pts[:, :1]
            return numpy.abs(numpy.linalg.det(v)) / 6.0
        raise RuntimeError(f"Unsupported simplex size: {self.simplex_cells.shape[1]}.")
