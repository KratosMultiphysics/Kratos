---
type: C++ Method
title: GetCoordinatesEntity
resource: processes/assign_scalar_input_to_entities_process.cpp#L178-L181
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

`array_1d<double, 3> AssignScalarInputToEntitiesProcess<Condition, AssignScalarInputToEntitiesProcessSettings::SaveAsNonHistoricalVariable>:: GetCoordinatesEntity(const IndexType Id)`

# Calls

- [GetGeometry](../../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)