---
type: C++ Method
title: CalculateDiscontinuousShapeFunctions
resource: utilities/discont_utils.h#L73-L430
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

`static int CalculateDiscontinuousShapeFunctions( BoundedMatrix<double,(3+1), 3 >& rPoints, BoundedMatrix<double, (3+1), 3 >& DN_DX, array_1d<double,(3+1)>& rDistances, array_1d<double,(3*(3-1))>& rVolumes, BoundedMatrix<double, 3*(3-1), (3+1) >& rShapeFunctionValues, array_1d<double,(3*(3-1))>& rPartitionsSign, std::vector<Matrix>& rGradientsValue, BoundedMatrix<double,3*(3-1), (3+1)>& Nenriched, array_1d<double,6>& edge_areas )`

# Calls

- [prod](../../../../functions/includes/amatrix_interface/prod.md)
- [norm_2](../../../../functions/includes/amatrix_interface/norm_2.md)
- [TetrahedraSplitMode](../../../../functions/utilities/split_tetrahedra/TetrahedraSplit/TetrahedraSplitMode.md)
- [Split_Tetrahedra](../../../../functions/utilities/split_tetrahedra/TetrahedraSplit/Split_Tetrahedra.md)
- [TetrahedraGetNewConnectivityGID](../../../../functions/utilities/split_tetrahedra/TetrahedraSplit/TetrahedraGetNewConnectivityGID.md)