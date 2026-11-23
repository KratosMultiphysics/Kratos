---
type: C++ Method
title: GetOffDiagonalConsistentMass
resource: containers/edge_based_data_structure.h#L88-L91
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/processes/flux_corrected_transport_convection_process/ExplicitSolve
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/flux_corrected_transport_convection_process/FluxCorrectedTransportConvectionProcess/CalculateResidual
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/flux_corrected_transport_convection_process/FluxCorrectedTransportConvectionProcess/CalculateAntidiffusiveEdgeContributions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/containers/test_edge_based_data_structure/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`double GetOffDiagonalConsistentMass() const`

# Called by

- [ExplicitSolve](../../../../functions/processes/flux_corrected_transport_convection_process/ExplicitSolve.md)
- [CalculateResidual](../../../../functions/processes/flux_corrected_transport_convection_process/FluxCorrectedTransportConvectionProcess/CalculateResidual.md)
- [CalculateAntidiffusiveEdgeContributions](../../../../functions/processes/flux_corrected_transport_convection_process/FluxCorrectedTransportConvectionProcess/CalculateAntidiffusiveEdgeContributions.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/containers/test_edge_based_data_structure/KRATOS_TEST_CASE_IN_SUITE.md)