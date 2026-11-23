---
type: C++ Method
title: GetDofList
resource: elements/embedded_nodal_variable_calculation_element_simplex.cpp#L218-L229
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void EmbeddedNodalVariableCalculationElementSimplex<double>::GetDofList( DofsVectorType &rElementalDofList, const ProcessInfo &rCurrentProcessInfo) const`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)