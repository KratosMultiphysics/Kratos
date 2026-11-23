---
type: C++ Module
title: ublas_space
resource: spaces/ublas_space.h#L1-L1148
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/numeric
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-define-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-process-info-h
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
  - target: external/utilities-parallel-utilities-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-reduction-utilities-h
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [MultValueNoAdd](../../classes/spaces/ublas_space/MultValueNoAdd.md)
- [MultValueNoAdd](../../functions/spaces/ublas_space/MultValueNoAdd/MultValueNoAdd.md)
- [MultAndAddValue](../../classes/spaces/ublas_space/MultAndAddValue.md)
- [MultAndAddValue](../../functions/spaces/ublas_space/MultAndAddValue/MultAndAddValue.md)
- [UblasSpace](../../functions/spaces/ublas_space/UblasSpace.md)
- [CreateEmptyMatrixPointer](../../functions/spaces/ublas_space/CreateEmptyMatrixPointer.md)
- [CreateEmptyVectorPointer](../../functions/spaces/ublas_space/CreateEmptyVectorPointer.md)
- [IsNull](../../functions/spaces/ublas_space/IsNull.md)
- [Size](../../functions/spaces/ublas_space/Size.md)
- [Size1](../../functions/spaces/ublas_space/Size1.md)
- [Size2](../../functions/spaces/ublas_space/Size2.md)
- [GetColumn](../../functions/spaces/ublas_space/GetColumn.md)
- [SetColumn](../../functions/spaces/ublas_space/SetColumn.md)
- [Copy](../../functions/spaces/ublas_space/Copy.md)
- [Copy](../../functions/spaces/ublas_space/Copy-2.md)
- [Dot](../../functions/spaces/ublas_space/Dot.md)
- [TwoNorm](../../functions/spaces/ublas_space/TwoNorm.md)
- [TwoNorm](../../functions/spaces/ublas_space/TwoNorm-2.md)
- [TwoNorm](../../functions/spaces/ublas_space/TwoNorm-3.md)
- [JacobiNorm](../../functions/spaces/ublas_space/JacobiNorm.md)
- [JacobiNorm](../../functions/spaces/ublas_space/JacobiNorm-2.md)
- [Mult](../../functions/spaces/ublas_space/Mult.md)
- [Mult](../../functions/spaces/ublas_space/Mult-2.md)
- [TransposeMult](../../functions/spaces/ublas_space/TransposeMult.md)
- [GraphDegree](../../functions/spaces/ublas_space/GraphDegree.md)
- [InplaceMult](../../functions/spaces/ublas_space/InplaceMult.md)
- [Assign](../../functions/spaces/ublas_space/Assign.md)
- [UnaliasedAdd](../../functions/spaces/ublas_space/UnaliasedAdd.md)
- [ScaleAndAdd](../../functions/spaces/ublas_space/ScaleAndAdd.md)
- [ScaleAndAdd](../../functions/spaces/ublas_space/ScaleAndAdd-2.md)
- [RowDot](../../functions/spaces/ublas_space/RowDot.md)
- [SetValue](../../functions/spaces/ublas_space/SetValue.md)
- [Set](../../functions/spaces/ublas_space/Set.md)
- [Resize](../../functions/spaces/ublas_space/Resize.md)
- [Resize](../../functions/spaces/ublas_space/Resize-2.md)
- [Resize](../../functions/spaces/ublas_space/Resize-3.md)
- [Resize](../../functions/spaces/ublas_space/Resize-4.md)
- [Clear](../../functions/spaces/ublas_space/Clear.md)
- [Clear](../../functions/spaces/ublas_space/Clear-2.md)
- [ResizeData](../../functions/spaces/ublas_space/ResizeData.md)
- [ResizeData](../../functions/spaces/ublas_space/ResizeData-2.md)
- [ResizeData](../../functions/spaces/ublas_space/ResizeData-3.md)
- [SetToZero](../../functions/spaces/ublas_space/SetToZero.md)
- [SetToZero](../../functions/spaces/ublas_space/SetToZero-2.md)
- [SetToZero](../../functions/spaces/ublas_space/SetToZero-3.md)
- [AssembleLHS](../../functions/spaces/ublas_space/AssembleLHS.md)
- [CheckAndCorrectZeroDiagonalValues](../../functions/spaces/ublas_space/CheckAndCorrectZeroDiagonalValues.md)
- [GetScaleNorm](../../functions/spaces/ublas_space/GetScaleNorm.md)
- [GetDiagonalNorm](../../functions/spaces/ublas_space/GetDiagonalNorm.md)
- [GetAveragevalueDiagonal](../../functions/spaces/ublas_space/GetAveragevalueDiagonal.md)
- [GetMaxDiagonal](../../functions/spaces/ublas_space/GetMaxDiagonal.md)
- [GetMinDiagonal](../../functions/spaces/ublas_space/GetMinDiagonal.md)
- [Info](../../functions/spaces/ublas_space/Info.md)
- [PrintInfo](../../functions/spaces/ublas_space/PrintInfo.md)
- [PrintData](../../functions/spaces/ublas_space/PrintData.md)
- [IsDistributed](../../functions/spaces/ublas_space/IsDistributed.md)
- [FastestDirectSolverList](../../functions/spaces/ublas_space/FastestDirectSolverList.md)
- [GetValue](../../functions/spaces/ublas_space/GetValue.md)
- [GatherValues](../../functions/spaces/ublas_space/GatherValues.md)
- [WriteMatrixMarketMatrix](../../functions/spaces/ublas_space/WriteMatrixMarketMatrix.md)
- [WriteMatrixMarketVector](../../functions/spaces/ublas_space/WriteMatrixMarketVector.md)
- [CreateDofUpdater](../../functions/spaces/ublas_space/CreateDofUpdater.md)
- [ParallelProductNoAdd](../../functions/spaces/ublas_space/ParallelProductNoAdd.md)
- [CreatePartition](../../functions/spaces/ublas_space/CreatePartition.md)
- [partial_product_no_add](../../functions/spaces/ublas_space/partial_product_no_add.md)

# Imports

- `<numeric>`
- `"includes/define.h"`
- `"includes/process_info.h"`
- `"includes/ublas_interface.h"`
- `"includes/matrix_market_interface.h"`
- `"utilities/dof_updater.h"`
- `"utilities/parallel_utilities.h"`
- `"utilities/reduction_utilities.h"`

# Member of

- [KratosMultiphysics](../../packages/KratosMultiphysics.md)