---
type: C++ Method
title: ComputeNodalSizes
resource: utilities/variable_redistribution_utility.cpp#L492-L514
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
  called_by:
  - target: functions/utilities/variable_redistribution_utility/VariableRedistributionUtility/SpecializedDistributePointValues
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void VariableRedistributionUtility::ComputeNodalSizes( ModelPart& rModelPart, TContainerType& rEntitiesContainer)`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [PointsNumber](../../../../functions/geometries/geometry/PointsNumber.md)

# Called by

- [SpecializedDistributePointValues](../../../../functions/utilities/variable_redistribution_utility/VariableRedistributionUtility/SpecializedDistributePointValues.md)