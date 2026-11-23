---
type: C++ Module
title: distributed_csr_matrix
resource: containers/distributed_csr_matrix.h#L1-L1147
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/iostream
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-define-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/containers-csr-matrix-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/containers-distributed-sparse-graph-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/containers-distributed-vector-importer-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/containers-distributed-vector-exporter-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-key-hash-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-atomic-utilities-h
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [DistributedCsrMatrix](../../classes/containers/distributed_csr_matrix/DistributedCsrMatrix.md)
- [DistributedCsrMatrix](../../functions/containers/distributed_csr_matrix/DistributedCsrMatrix/DistributedCsrMatrix.md)
- [DistributedCsrMatrix](../../functions/containers/distributed_csr_matrix/DistributedCsrMatrix/DistributedCsrMatrix-2.md)
- [DistributedCsrMatrix](../../functions/containers/distributed_csr_matrix/DistributedCsrMatrix/DistributedCsrMatrix-3.md)
- [DistributedCsrMatrix](../../functions/containers/distributed_csr_matrix/DistributedCsrMatrix/DistributedCsrMatrix-4.md)
- [DistributedCsrMatrix](../../functions/containers/distributed_csr_matrix/DistributedCsrMatrix/DistributedCsrMatrix-5.md)
- [Clear](../../functions/containers/distributed_csr_matrix/DistributedCsrMatrix/Clear.md)
- [SetValue](../../functions/containers/distributed_csr_matrix/DistributedCsrMatrix/SetValue.md)
- [local_size1](../../functions/containers/distributed_csr_matrix/DistributedCsrMatrix/local_size1.md)
- [size2](../../functions/containers/distributed_csr_matrix/DistributedCsrMatrix/size2.md)
- [local_nnz](../../functions/containers/distributed_csr_matrix/DistributedCsrMatrix/local_nnz.md)
- [pGetComm](../../functions/containers/distributed_csr_matrix/DistributedCsrMatrix/pGetComm.md)
- [GetOffDiagonalBlockLocalId](../../functions/containers/distributed_csr_matrix/DistributedCsrMatrix/GetOffDiagonalBlockLocalId.md)
- [GetOffDiaGlobalId](../../functions/containers/distributed_csr_matrix/DistributedCsrMatrix/GetOffDiaGlobalId.md)
- [GetDiagonalIndex2DataInGlobalNumbering](../../functions/containers/distributed_csr_matrix/DistributedCsrMatrix/GetDiagonalIndex2DataInGlobalNumbering.md)
- [GetOffDiagonalIndex2DataInGlobalNumbering](../../functions/containers/distributed_csr_matrix/DistributedCsrMatrix/GetOffDiagonalIndex2DataInGlobalNumbering.md)
- [ApplyHomogeneousDirichlet](../../functions/containers/distributed_csr_matrix/DistributedCsrMatrix/ApplyHomogeneousDirichlet.md)
- [SpMV](../../functions/containers/distributed_csr_matrix/DistributedCsrMatrix/SpMV.md)
- [SpMV](../../functions/containers/distributed_csr_matrix/DistributedCsrMatrix/SpMV-2.md)
- [TransposeSpMV](../../functions/containers/distributed_csr_matrix/DistributedCsrMatrix/TransposeSpMV.md)
- [TransposeSpMV](../../functions/containers/distributed_csr_matrix/DistributedCsrMatrix/TransposeSpMV-2.md)
- [NormFrobenius](../../functions/containers/distributed_csr_matrix/DistributedCsrMatrix/NormFrobenius.md)
- [BeginAssemble](../../functions/containers/distributed_csr_matrix/DistributedCsrMatrix/BeginAssemble.md)
- [FinalizeAssemble](../../functions/containers/distributed_csr_matrix/DistributedCsrMatrix/FinalizeAssemble.md)
- [Assemble](../../functions/containers/distributed_csr_matrix/Assemble.md)
- [AssembleEntry](../../functions/containers/distributed_csr_matrix/DistributedCsrMatrix/AssembleEntry.md)
- [Assemble](../../functions/containers/distributed_csr_matrix/Assemble-2.md)
- [ToMap](../../functions/containers/distributed_csr_matrix/DistributedCsrMatrix/ToMap.md)
- [ToSerialCSR](../../functions/containers/distributed_csr_matrix/DistributedCsrMatrix/ToSerialCSR.md)
- [NormDiagonal](../../functions/containers/distributed_csr_matrix/DistributedCsrMatrix/NormDiagonal.md)
- [MaxDiagonal](../../functions/containers/distributed_csr_matrix/DistributedCsrMatrix/MaxDiagonal.md)
- [MinDiagonal](../../functions/containers/distributed_csr_matrix/DistributedCsrMatrix/MinDiagonal.md)
- [Info](../../functions/containers/distributed_csr_matrix/DistributedCsrMatrix/Info.md)
- [PrintInfo](../../functions/containers/distributed_csr_matrix/DistributedCsrMatrix/PrintInfo.md)
- [PrintData](../../functions/containers/distributed_csr_matrix/DistributedCsrMatrix/PrintData.md)
- [BinarySearch](../../functions/containers/distributed_csr_matrix/BinarySearch.md)
- [PrepareNonLocalCommunications](../../functions/containers/distributed_csr_matrix/DistributedCsrMatrix/PrepareNonLocalCommunications.md)
- [save](../../functions/containers/distributed_csr_matrix/DistributedCsrMatrix/save.md)
- [load](../../functions/containers/distributed_csr_matrix/DistributedCsrMatrix/load.md)
- [ReconstructDirectAccessVectors](../../functions/containers/distributed_csr_matrix/DistributedCsrMatrix/ReconstructDirectAccessVectors.md)

# Imports

- `<iostream>`
- `"includes/define.h"`
- `"containers/csr_matrix.h"`
- `"containers/distributed_sparse_graph.h"`
- `"containers/distributed_vector_importer.h"`
- `"containers/distributed_vector_exporter.h"`
- `"includes/key_hash.h"`
- `"utilities/atomic_utilities.h"`

# Member of

- [KratosMultiphysics](../../packages/KratosMultiphysics.md)