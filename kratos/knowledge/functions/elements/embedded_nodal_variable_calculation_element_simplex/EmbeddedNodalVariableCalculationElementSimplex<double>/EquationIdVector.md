---
type: C++ Method
title: EquationIdVector
resource: elements/embedded_nodal_variable_calculation_element_simplex.cpp#L186-L198
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/node/GetDofPosition
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void EmbeddedNodalVariableCalculationElementSimplex<double>::EquationIdVector( EquationIdVectorType &rResult, const ProcessInfo &rCurrentProcessInfo) const`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [GetDofPosition](../../../../functions/includes/node/GetDofPosition.md)