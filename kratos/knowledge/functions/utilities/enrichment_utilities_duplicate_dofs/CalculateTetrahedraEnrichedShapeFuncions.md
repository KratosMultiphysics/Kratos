---
type: C++ Function
title: CalculateTetrahedraEnrichedShapeFuncions
resource: utilities/enrichment_utilities_duplicate_dofs.h#L89-L488
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/amatrix_interface/prod
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/amatrix_interface/norm_2
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/split_tetrahedra/TetrahedraSplit/TetrahedraSplitMode
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/split_tetrahedra/TetrahedraSplit/Split_Tetrahedra
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/split_tetrahedra/TetrahedraSplit/TetrahedraGetNewConnectivityGID
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static int CalculateTetrahedraEnrichedShapeFuncions(TMatrixType const& rPoints, TGradientType const& DN_DX, TVectorType& rDistances, TVectorType& rVolumes, TMatrixType& rShapeFunctionValues, TVectorType& rPartitionsSign, std::vector<TMatrixType>& rGradientsValue, TMatrixType& NEnriched)`

# Calls

- [prod](../../../functions/includes/amatrix_interface/prod.md)
- [norm_2](../../../functions/includes/amatrix_interface/norm_2.md)
- [TetrahedraSplitMode](../../../functions/utilities/split_tetrahedra/TetrahedraSplit/TetrahedraSplitMode.md)
- [Split_Tetrahedra](../../../functions/utilities/split_tetrahedra/TetrahedraSplit/Split_Tetrahedra.md)
- [TetrahedraGetNewConnectivityGID](../../../functions/utilities/split_tetrahedra/TetrahedraSplit/TetrahedraGetNewConnectivityGID.md)