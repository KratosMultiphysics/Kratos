---
type: C++ Method
title: CallSpecializedConvertDistributedValuesToPoint
resource: utilities/variable_redistribution_utility.cpp#L235-L280
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/geometry/PointsNumber
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void VariableRedistributionUtility::CallSpecializedConvertDistributedValuesToPoint( ModelPart& rModelPart, TContainerType& rEntitiesContainer, const Variable<TValueType>& rDistributedVariable, const Variable<TValueType>& rPointVariable)`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [PointsNumber](../../../../functions/geometries/geometry/PointsNumber.md)