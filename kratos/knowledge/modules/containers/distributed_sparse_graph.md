---
type: C++ Module
title: distributed_sparse_graph
resource: containers/distributed_sparse_graph.h#L1-L473
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/iostream
    resolved_by: tree-sitter
    confidence: exact
  - target: external/mutex
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-ublas-interface-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-serializer-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-parallel-environment-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/containers-distributed-numbering-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-communication-coloring-utilities-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/containers-sparse-graph-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/containers-sparse-contiguous-row-graph-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-parallel-utilities-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/unordered-map
    resolved_by: tree-sitter
    confidence: exact
  - target: external/unordered-set
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-define-h
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [DistributedSparseGraph](../../classes/containers/distributed_sparse_graph/DistributedSparseGraph.md)
- [DistributedSparseGraph](../../functions/containers/distributed_sparse_graph/DistributedSparseGraph/DistributedSparseGraph.md)
- [pGetComm](../../functions/containers/distributed_sparse_graph/DistributedSparseGraph/pGetComm.md)
- [Size](../../functions/containers/distributed_sparse_graph/DistributedSparseGraph/Size.md)
- [LocalSize](../../functions/containers/distributed_sparse_graph/DistributedSparseGraph/LocalSize.md)
- [Has](../../functions/containers/distributed_sparse_graph/DistributedSparseGraph/Has.md)
- [ComputeLocalMinMaxColumnIndex](../../functions/containers/distributed_sparse_graph/DistributedSparseGraph/ComputeLocalMinMaxColumnIndex.md)
- [ComputeMaxGlobalColumnIndex](../../functions/containers/distributed_sparse_graph/DistributedSparseGraph/ComputeMaxGlobalColumnIndex.md)
- [Clear](../../functions/containers/distributed_sparse_graph/DistributedSparseGraph/Clear.md)
- [AddEntry](../../functions/containers/distributed_sparse_graph/DistributedSparseGraph/AddEntry.md)
- [AddEntries](../../functions/containers/distributed_sparse_graph/AddEntries.md)
- [AddEntries](../../functions/containers/distributed_sparse_graph/AddEntries-2.md)
- [AddEntries](../../functions/containers/distributed_sparse_graph/AddEntries-3.md)
- [AddEntries](../../functions/containers/distributed_sparse_graph/AddEntries-4.md)
- [Finalize](../../functions/containers/distributed_sparse_graph/DistributedSparseGraph/Finalize.md)
- [Info](../../functions/containers/distributed_sparse_graph/DistributedSparseGraph/Info.md)
- [PrintInfo](../../functions/containers/distributed_sparse_graph/DistributedSparseGraph/PrintInfo.md)
- [PrintData](../../functions/containers/distributed_sparse_graph/DistributedSparseGraph/PrintData.md)
- [save](../../functions/containers/distributed_sparse_graph/DistributedSparseGraph/save.md)
- [load](../../functions/containers/distributed_sparse_graph/DistributedSparseGraph/load.md)
- [DistributedSparseGraph](../../functions/containers/distributed_sparse_graph/DistributedSparseGraph/DistributedSparseGraph-2.md)

# Imports

- `<iostream>`
- `<mutex>`
- `"includes/ublas_interface.h"`
- `"includes/serializer.h"`
- `"includes/parallel_environment.h"`
- `"containers/distributed_numbering.h"`
- `"utilities/communication_coloring_utilities.h"`
- `"containers/sparse_graph.h"`
- `"containers/sparse_contiguous_row_graph.h"`
- `"utilities/parallel_utilities.h"`
- `<unordered_map>`
- `<unordered_set>`
- `"includes/define.h"`

# Member of

- [KratosMultiphysics](../../packages/KratosMultiphysics.md)