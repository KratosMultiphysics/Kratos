---
type: C++ Method
title: GetLumpedMassMatrixDiagonal
resource: containers/edge_based_data_structure.h#L288-L291
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/processes/flux_corrected_transport_convection_process/ExplicitSolve
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/flux_corrected_transport_convection_process/FluxCorrectedTransportConvectionProcess/EvaluateLimiter
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/containers/test_edge_based_data_structure/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`double GetLumpedMassMatrixDiagonal(IndexType I) const`

# Called by

- [ExplicitSolve](../../../../functions/processes/flux_corrected_transport_convection_process/ExplicitSolve.md)
- [EvaluateLimiter](../../../../functions/processes/flux_corrected_transport_convection_process/FluxCorrectedTransportConvectionProcess/EvaluateLimiter.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/containers/test_edge_based_data_structure/KRATOS_TEST_CASE_IN_SUITE.md)