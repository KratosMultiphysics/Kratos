---
type: C++ Function
title: KRATOS_TEST_CASE_IN_SUITE
resource: tests/cpp_tests/containers/test_edge_based_data_structure.cpp#L43-L113
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/containers/edge_based_data_structure/EdgeBasedDataStructure/CalculateEdgeDataStructure
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/edge_based_data_structure/EdgeData/IsBoundary
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/edge_based_data_structure/EdgeData/GetOffDiagonalConsistentMass
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/edge_based_data_structure/EdgeData/GetOffDiagonalLaplacian
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/edge_based_data_structure/EdgeBasedDataStructure/GetMassMatrixDiagonal
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/edge_based_data_structure/EdgeBasedDataStructure/GetLumpedMassMatrixDiagonal
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/containers/test_edge_based_data_structure/CheckIndicesVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/edge_based_data_structure/EdgeBasedDataStructure/NumberOfEdges
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`KRATOS_TEST_CASE_IN_SUITE(EdgeBasedDataStructure2D, KratosCoreFastSuite)`

# Calls

- [CalculateEdgeDataStructure](../../../../../functions/containers/edge_based_data_structure/EdgeBasedDataStructure/CalculateEdgeDataStructure.md)
- [IsBoundary](../../../../../functions/containers/edge_based_data_structure/EdgeData/IsBoundary.md)
- [GetOffDiagonalConsistentMass](../../../../../functions/containers/edge_based_data_structure/EdgeData/GetOffDiagonalConsistentMass.md)
- [GetOffDiagonalLaplacian](../../../../../functions/containers/edge_based_data_structure/EdgeData/GetOffDiagonalLaplacian.md)
- [GetMassMatrixDiagonal](../../../../../functions/containers/edge_based_data_structure/EdgeBasedDataStructure/GetMassMatrixDiagonal.md)
- [GetLumpedMassMatrixDiagonal](../../../../../functions/containers/edge_based_data_structure/EdgeBasedDataStructure/GetLumpedMassMatrixDiagonal.md)
- [CheckIndicesVector](../../../../../functions/tests/cpp_tests/containers/test_edge_based_data_structure/CheckIndicesVector.md)
- [NumberOfEdges](../../../../../functions/containers/edge_based_data_structure/EdgeBasedDataStructure/NumberOfEdges.md)