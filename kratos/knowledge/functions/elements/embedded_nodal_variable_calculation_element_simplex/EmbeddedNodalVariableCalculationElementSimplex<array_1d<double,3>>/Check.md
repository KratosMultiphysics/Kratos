---
type: C++ Method
title: Check
resource: elements/embedded_nodal_variable_calculation_element_simplex.cpp#L59-L77
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/node/SolutionStepsDataHas
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`int EmbeddedNodalVariableCalculationElementSimplex<array_1d<double,3>>::Check(const ProcessInfo &rCurrentProcessInfo) const`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [SolutionStepsDataHas](../../../../functions/includes/node/SolutionStepsDataHas.md)