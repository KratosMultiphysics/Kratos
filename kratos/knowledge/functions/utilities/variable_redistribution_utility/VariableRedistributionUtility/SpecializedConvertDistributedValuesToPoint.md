---
type: C++ Method
title: SpecializedConvertDistributedValuesToPoint
resource: utilities/variable_redistribution_utility.cpp#L345-L387
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

`void VariableRedistributionUtility::SpecializedConvertDistributedValuesToPoint( ModelPart& rModelPart, TContainerType& rEntitiesContainer, const Variable< TValueType >& rDistributedVariable, const Variable< TValueType >& rPointVariable)`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)