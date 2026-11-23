---
type: C++ Method
title: GetDistanceBasedShapeFunctionValues
resource: elements/embedded_nodal_variable_calculation_element_simplex.cpp#L248-L255
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/elements/embedded_nodal_variable_calculation_element_simplex/EmbeddedNodalVariableCalculationElementSimplex<double>/CalculateLeftHandSide
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/elements/embedded_nodal_variable_calculation_element_simplex/EmbeddedNodalVariableCalculationElementSimplex<array_1d<double,3>>/CalculateLeftHandSide
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/elements/embedded_nodal_variable_calculation_element_simplex/EmbeddedNodalVariableCalculationElementSimplex<double>/CalculateRightHandSide
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/elements/embedded_nodal_variable_calculation_element_simplex/EmbeddedNodalVariableCalculationElementSimplex<array_1d<double, 3>>/CalculateRightHandSide
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`const array_1d<double, 2> EmbeddedNodalVariableCalculationElementSimplex<TVarType>::GetDistanceBasedShapeFunctionValues()`

# Called by

- [CalculateLeftHandSide](../../../../functions/elements/embedded_nodal_variable_calculation_element_simplex/EmbeddedNodalVariableCalculationElementSimplex<double>/CalculateLeftHandSide.md)
- [CalculateLeftHandSide](../../../../functions/elements/embedded_nodal_variable_calculation_element_simplex/EmbeddedNodalVariableCalculationElementSimplex<array_1d<double,3>>/CalculateLeftHandSide.md)
- [CalculateRightHandSide](../../../../functions/elements/embedded_nodal_variable_calculation_element_simplex/EmbeddedNodalVariableCalculationElementSimplex<double>/CalculateRightHandSide.md)
- [CalculateRightHandSide](../../../../functions/elements/embedded_nodal_variable_calculation_element_simplex/EmbeddedNodalVariableCalculationElementSimplex<array_1d<double, 3>>/CalculateRightHandSide.md)