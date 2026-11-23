---
type: C++ Method
title: CalculateSystemContributions
resource: solving_strategies/schemes/residualbased_incrementalupdate_static_scheme_slip.h#L189-L203
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CalculateSystemContributions(Element& rCurrentElement, LocalSystemMatrixType& LHS_Contribution, LocalSystemVectorType& RHS_Contribution, Element::EquationIdVectorType& EquationId, const ProcessInfo& CurrentProcessInfo) override`

# Calls

- [GetGeometry](../../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)