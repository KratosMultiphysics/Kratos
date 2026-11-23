---
type: C++ Method
title: InternalAssignValueScalar
resource: processes/assign_scalar_field_to_entities_process.cpp#L546-L595
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/function_parser_utility/DependsOnSpace
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/assign_scalar_field_to_entities_process/AssignScalarFieldToEntitiesProcess<TEntity, THistorical>/Execute
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void AssignScalarFieldToEntitiesProcess<TEntity, THistorical>::InternalAssignValueScalar( const Variable<double>& rVar, const double Time )`

# Calls

- [DependsOnSpace](../../../../functions/utilities/function_parser_utility/DependsOnSpace.md)

# Called by

- [Execute](../../../../functions/processes/assign_scalar_field_to_entities_process/AssignScalarFieldToEntitiesProcess<TEntity, THistorical>/Execute.md)