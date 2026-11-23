---
type: C++ Module
title: kratos_space
resource: spaces/kratos_space.h#L1-L645
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/string
    resolved_by: tree-sitter
    confidence: exact
  - target: external/iostream
    resolved_by: tree-sitter
    confidence: exact
  - target: external/cstddef
    resolved_by: tree-sitter
    confidence: exact
  - target: external/numeric
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-define-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-ublas-interface-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-matrix-market-interface-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-dof-updater-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/containers-csr-matrix-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/containers-system-vector-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/containers-distributed-csr-matrix-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-parallel-utilities-h
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

- [KratosSpace](../../classes/spaces/kratos_space/KratosSpace.md)
- [KratosSpace](../../functions/spaces/kratos_space/KratosSpace/KratosSpace.md)
- [CreateEmptyMatrixPointer](../../functions/spaces/kratos_space/KratosSpace/CreateEmptyMatrixPointer.md)
- [CreateEmptyVectorPointer](../../functions/spaces/kratos_space/KratosSpace/CreateEmptyVectorPointer.md)
- [Size](../../functions/spaces/kratos_space/KratosSpace/Size.md)
- [Size1](../../functions/spaces/kratos_space/KratosSpace/Size1.md)
- [Size2](../../functions/spaces/kratos_space/KratosSpace/Size2.md)
- [GetColumn](../../functions/spaces/kratos_space/GetColumn.md)
- [SetColumn](../../functions/spaces/kratos_space/SetColumn.md)
- [Copy](../../functions/spaces/kratos_space/KratosSpace/Copy.md)
- [Copy](../../functions/spaces/kratos_space/KratosSpace/Copy-2.md)
- [Dot](../../functions/spaces/kratos_space/KratosSpace/Dot.md)
- [TwoNorm](../../functions/spaces/kratos_space/KratosSpace/TwoNorm.md)
- [TwoNorm](../../functions/spaces/kratos_space/KratosSpace/TwoNorm-2.md)
- [TwoNorm](../../functions/spaces/kratos_space/KratosSpace/TwoNorm-3.md)
- [TwoNorm](../../functions/spaces/kratos_space/KratosSpace/TwoNorm-4.md)
- [JacobiNorm](../../functions/spaces/kratos_space/KratosSpace/JacobiNorm.md)
- [JacobiNorm](../../functions/spaces/kratos_space/KratosSpace/JacobiNorm-2.md)
- [Mult](../../functions/spaces/kratos_space/KratosSpace/Mult.md)
- [Mult](../../functions/spaces/kratos_space/KratosSpace/Mult-2.md)
- [TransposeMult](../../functions/spaces/kratos_space/TransposeMult.md)
- [GraphDegree](../../functions/spaces/kratos_space/KratosSpace/GraphDegree.md)
- [GraphNeighbors](../../functions/spaces/kratos_space/KratosSpace/GraphNeighbors.md)
- [InplaceMult](../../functions/spaces/kratos_space/KratosSpace/InplaceMult.md)
- [Assign](../../functions/spaces/kratos_space/KratosSpace/Assign.md)
- [UnaliasedAdd](../../functions/spaces/kratos_space/KratosSpace/UnaliasedAdd.md)
- [ScaleAndAdd](../../functions/spaces/kratos_space/KratosSpace/ScaleAndAdd.md)
- [ScaleAndAdd](../../functions/spaces/kratos_space/KratosSpace/ScaleAndAdd-2.md)
- [RowDot](../../functions/spaces/kratos_space/KratosSpace/RowDot.md)
- [SetValue](../../functions/spaces/kratos_space/KratosSpace/SetValue.md)
- [Set](../../functions/spaces/kratos_space/KratosSpace/Set.md)
- [Resize](../../functions/spaces/kratos_space/KratosSpace/Resize.md)
- [Resize](../../functions/spaces/kratos_space/KratosSpace/Resize-2.md)
- [Resize](../../functions/spaces/kratos_space/KratosSpace/Resize-3.md)
- [Resize](../../functions/spaces/kratos_space/KratosSpace/Resize-4.md)
- [Clear](../../functions/spaces/kratos_space/KratosSpace/Clear.md)
- [Clear](../../functions/spaces/kratos_space/KratosSpace/Clear-2.md)
- [ResizeData](../../functions/spaces/kratos_space/ResizeData.md)
- [ResizeData](../../functions/spaces/kratos_space/KratosSpace/ResizeData.md)
- [ResizeData](../../functions/spaces/kratos_space/KratosSpace/ResizeData-2.md)
- [SetToZero](../../functions/spaces/kratos_space/SetToZero.md)
- [SetToZero](../../functions/spaces/kratos_space/KratosSpace/SetToZero.md)
- [SetToZero](../../functions/spaces/kratos_space/KratosSpace/SetToZero-2.md)
- [AssembleLHS](../../functions/spaces/kratos_space/AssembleLHS.md)
- [Info](../../functions/spaces/kratos_space/KratosSpace/Info.md)
- [PrintInfo](../../functions/spaces/kratos_space/KratosSpace/PrintInfo.md)
- [PrintData](../../functions/spaces/kratos_space/KratosSpace/PrintData.md)
- [IsDistributed](../../functions/spaces/kratos_space/KratosSpace/IsDistributed.md)
- [GetValue](../../functions/spaces/kratos_space/KratosSpace/GetValue.md)
- [GatherValues](../../functions/spaces/kratos_space/KratosSpace/GatherValues.md)
- [WriteMatrixMarketMatrix](../../functions/spaces/kratos_space/WriteMatrixMarketMatrix.md)
- [WriteMatrixMarketVector](../../functions/spaces/kratos_space/WriteMatrixMarketVector.md)
- [CreateDofUpdater](../../functions/spaces/kratos_space/KratosSpace/CreateDofUpdater.md)
- [KratosSpace](../../functions/spaces/kratos_space/KratosSpace/KratosSpace-2.md)

# Imports

- `<string>`
- `<iostream>`
- `<cstddef>`
- `<numeric>`
- `"includes/define.h"`
- `"includes/ublas_interface.h"`
- `"includes/matrix_market_interface.h"`
- `"utilities/dof_updater.h"`
- `"containers/csr_matrix.h"`
- `"containers/system_vector.h"`
- `"containers/distributed_csr_matrix.h"`
- `"utilities/parallel_utilities.h"`
- `"utilities/atomic_utilities.h"`

# Member of

- [KratosMultiphysics](../../packages/KratosMultiphysics.md)