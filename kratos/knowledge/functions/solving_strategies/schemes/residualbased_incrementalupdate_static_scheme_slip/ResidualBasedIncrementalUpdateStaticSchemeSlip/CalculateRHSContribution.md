---
type: C++ Method
title: CalculateRHSContribution
resource: solving_strategies/schemes/residualbased_incrementalupdate_static_scheme_slip.h#L206-L219
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CalculateRHSContribution(Element& rCurrentElement, LocalSystemVectorType& RHS_Contribution, Element::EquationIdVectorType& EquationId, const ProcessInfo& CurrentProcessInfo) override`

# Calls

- [GetGeometry](../../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)