---
type: C++ Method
title: GetValueFromTable
resource: sources/table_accessor.cpp#L44-L78
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/node/SolutionStepsDataHas
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/sources/table_accessor/TableAccessor/GetValue
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`double TableAccessor::GetValueFromTable( const Variable<double> &rIndependentVariable, const Variable<double> &rDependentVariable, const Properties& rProperties, const GeometryType& rGeometry, const Vector& rShapeFunctionVector, const ProcessInfo& rProcessInfo) const`

# Calls

- [SolutionStepsDataHas](../../../../functions/includes/node/SolutionStepsDataHas.md)

# Called by

- [GetValue](../../../../functions/sources/table_accessor/TableAccessor/GetValue.md)