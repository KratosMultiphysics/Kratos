---
type: C++ Module
title: csr_matrix
resource: containers/csr_matrix.h#L1-L1045
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/iostream
    resolved_by: tree-sitter
    confidence: exact
  - target: external/limits
    resolved_by: tree-sitter
    confidence: exact
  - target: external/span-span-hpp
    resolved_by: tree-sitter
    confidence: exact
  - target: external/containers-sparse-contiguous-row-graph-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/containers-system-vector-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-parallel-utilities-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-atomic-utilities-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-reduction-utilities-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-key-hash-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-parallel-environment-h
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

- [CsrMatrix](../../classes/containers/csr_matrix/CsrMatrix.md)
- [CsrMatrix](../../functions/containers/csr_matrix/CsrMatrix/CsrMatrix.md)
- [CsrMatrix](../../functions/containers/csr_matrix/CsrMatrix/CsrMatrix-2.md)
- [CsrMatrix](../../functions/containers/csr_matrix/CsrMatrix.md)
- [CsrMatrix](../../functions/containers/csr_matrix/CsrMatrix/CsrMatrix-3.md)
- [CsrMatrix](../../functions/containers/csr_matrix/CsrMatrix/CsrMatrix-4.md)
- [CsrMatrix](../../functions/containers/csr_matrix/CsrMatrix/CsrMatrix-5.md)
- [Clear](../../functions/containers/csr_matrix/CsrMatrix/Clear.md)
- [pGetComm](../../functions/containers/csr_matrix/CsrMatrix/pGetComm.md)
- [SetValue](../../functions/containers/csr_matrix/CsrMatrix/SetValue.md)
- [size1](../../functions/containers/csr_matrix/CsrMatrix/size1.md)
- [size2](../../functions/containers/csr_matrix/CsrMatrix/size2.md)
- [nnz](../../functions/containers/csr_matrix/CsrMatrix/nnz.md)
- [IsOwnerOfData](../../functions/containers/csr_matrix/CsrMatrix/IsOwnerOfData.md)
- [SetIsOwnerOfData](../../functions/containers/csr_matrix/CsrMatrix/SetIsOwnerOfData.md)
- [SetColSize](../../functions/containers/csr_matrix/CsrMatrix/SetColSize.md)
- [SetRowSize](../../functions/containers/csr_matrix/CsrMatrix/SetRowSize.md)
- [ComputeColSize](../../functions/containers/csr_matrix/CsrMatrix/ComputeColSize.md)
- [AssignIndex1Data](../../functions/containers/csr_matrix/CsrMatrix/AssignIndex1Data.md)
- [AssignIndex2Data](../../functions/containers/csr_matrix/CsrMatrix/AssignIndex2Data.md)
- [AssignValueData](../../functions/containers/csr_matrix/CsrMatrix/AssignValueData.md)
- [ResizeIndex1Data](../../functions/containers/csr_matrix/CsrMatrix/ResizeIndex1Data.md)
- [ResizeIndex2Data](../../functions/containers/csr_matrix/CsrMatrix/ResizeIndex2Data.md)
- [ResizeValueData](../../functions/containers/csr_matrix/CsrMatrix/ResizeValueData.md)
- [CheckColSize](../../functions/containers/csr_matrix/CsrMatrix/CheckColSize.md)
- [FindValueIndex](../../functions/containers/csr_matrix/CsrMatrix/FindValueIndex.md)
- [Has](../../functions/containers/csr_matrix/CsrMatrix/Has.md)
- [SpMV](../../functions/containers/csr_matrix/SpMV.md)
- [SpMV](../../functions/containers/csr_matrix/SpMV-2.md)
- [TransposeSpMV](../../functions/containers/csr_matrix/TransposeSpMV.md)
- [TransposeSpMV](../../functions/containers/csr_matrix/TransposeSpMV-2.md)
- [NormFrobenius](../../functions/containers/csr_matrix/CsrMatrix/NormFrobenius.md)
- [reserve](../../functions/containers/csr_matrix/CsrMatrix/reserve.md)
- [ToMap](../../functions/containers/csr_matrix/CsrMatrix/ToMap.md)
- [BeginAssemble](../../functions/containers/csr_matrix/CsrMatrix/BeginAssemble.md)
- [FinalizeAssemble](../../functions/containers/csr_matrix/CsrMatrix/FinalizeAssemble.md)
- [Assemble](../../functions/containers/csr_matrix/Assemble.md)
- [SafeAssemble](../../functions/containers/csr_matrix/SafeAssemble.md)
- [Assemble](../../functions/containers/csr_matrix/Assemble-2.md)
- [AssembleEntry](../../functions/containers/csr_matrix/CsrMatrix/AssembleEntry.md)
- [ApplyHomogeneousDirichlet](../../functions/containers/csr_matrix/ApplyHomogeneousDirichlet.md)
- [NormDiagonal](../../functions/containers/csr_matrix/CsrMatrix/NormDiagonal.md)
- [MaxDiagonal](../../functions/containers/csr_matrix/CsrMatrix/MaxDiagonal.md)
- [MinDiagonal](../../functions/containers/csr_matrix/CsrMatrix/MinDiagonal.md)
- [Info](../../functions/containers/csr_matrix/CsrMatrix/Info.md)
- [PrintInfo](../../functions/containers/csr_matrix/CsrMatrix/PrintInfo.md)
- [PrintData](../../functions/containers/csr_matrix/CsrMatrix/PrintData.md)
- [BinarySearch](../../functions/containers/csr_matrix/BinarySearch.md)
- [save](../../functions/containers/csr_matrix/CsrMatrix/save.md)
- [load](../../functions/containers/csr_matrix/CsrMatrix/load.md)

# Imports

- `<iostream>`
- `<limits>`
- `<span/span.hpp>`
- `"containers/sparse_contiguous_row_graph.h"`
- `"containers/system_vector.h"`
- `"utilities/parallel_utilities.h"`
- `"utilities/atomic_utilities.h"`
- `"utilities/reduction_utilities.h"`
- `"includes/key_hash.h"`
- `"includes/parallel_environment.h"`
- `"includes/define.h"`

# Member of

- [KratosMultiphysics](../../packages/KratosMultiphysics.md)