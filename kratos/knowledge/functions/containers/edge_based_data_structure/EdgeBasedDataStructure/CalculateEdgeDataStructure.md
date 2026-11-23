---
type: C++ Method
title: CalculateEdgeDataStructure
resource: containers/edge_based_data_structure.h#L219-L232
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/containers/edge_based_data_structure/EdgeBasedDataStructure/CheckGeometryTypes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/edge_based_data_structure/EdgeBasedDataStructure/FillEdgesSparseGraph
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/edge_based_data_structure/EdgeBasedDataStructure/CalculateEdgeDataValues
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/flux_corrected_transport_convection_process/FluxCorrectedTransportConvectionProcess/ExecuteInitialize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/containers/test_edge_based_data_structure/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CalculateEdgeDataStructure(const ModelPart& rModelPart)`

# Calls

- [CheckGeometryTypes](../../../../functions/containers/edge_based_data_structure/EdgeBasedDataStructure/CheckGeometryTypes.md)
- [FillEdgesSparseGraph](../../../../functions/containers/edge_based_data_structure/EdgeBasedDataStructure/FillEdgesSparseGraph.md)
- [CalculateEdgeDataValues](../../../../functions/containers/edge_based_data_structure/EdgeBasedDataStructure/CalculateEdgeDataValues.md)

# Called by

- [ExecuteInitialize](../../../../functions/processes/flux_corrected_transport_convection_process/FluxCorrectedTransportConvectionProcess/ExecuteInitialize.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/containers/test_edge_based_data_structure/KRATOS_TEST_CASE_IN_SUITE.md)