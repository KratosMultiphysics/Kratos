---
type: C++ Method
title: CalculateRightHandSide
resource: elements/embedded_nodal_variable_calculation_element_simplex.cpp#L131-L154
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/elements/embedded_nodal_variable_calculation_element_simplex/EmbeddedNodalVariableCalculationElementSimplex<TVarType>/GetDistanceBasedShapeFunctionValues
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void EmbeddedNodalVariableCalculationElementSimplex<double>::CalculateRightHandSide( VectorType &rRigthHandSideVector, const ProcessInfo &rCurrentProcessInfo)`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [GetDistanceBasedShapeFunctionValues](../../../../functions/elements/embedded_nodal_variable_calculation_element_simplex/EmbeddedNodalVariableCalculationElementSimplex<TVarType>/GetDistanceBasedShapeFunctionValues.md)