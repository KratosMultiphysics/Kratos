---
type: C++ Method
title: GetMatrixRow
resource: solving_strategies/strategies/butcher_tableau.h#L117-L120
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/solving_strategies/strategies/butcher_tableau/ButcherTableau/GetMatrixRowBegin
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/strategies/butcher_tableau/ButcherTableau/GetMatrixRowEnd
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/flux_corrected_transport_convection_process/FluxCorrectedTransportConvectionProcess/SolveSubStep
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::tuple<RowType::const_iterator, RowType::const_iterator> GetMatrixRow(const unsigned int SubStepIndex) const`

# Calls

- [GetMatrixRowBegin](../../../../../functions/solving_strategies/strategies/butcher_tableau/ButcherTableau/GetMatrixRowBegin.md)
- [GetMatrixRowEnd](../../../../../functions/solving_strategies/strategies/butcher_tableau/ButcherTableau/GetMatrixRowEnd.md)

# Called by

- [SolveSubStep](../../../../../functions/processes/flux_corrected_transport_convection_process/FluxCorrectedTransportConvectionProcess/SolveSubStep.md)