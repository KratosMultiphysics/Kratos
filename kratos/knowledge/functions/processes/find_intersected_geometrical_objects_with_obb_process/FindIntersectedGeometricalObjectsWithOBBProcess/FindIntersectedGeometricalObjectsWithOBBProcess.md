---
type: C++ Method
title: FindIntersectedGeometricalObjectsWithOBBProcess
resource: processes/find_intersected_geometrical_objects_with_obb_process.cpp#L30-L53
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/model_part/ModelPart/CreateNewProperties
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/RecursivelyValidateAndAssignDefaults
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetBool
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetDouble
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/find_intersected_geometrical_objects_with_obb_process/ConvertIntersection
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`FindIntersectedGeometricalObjectsWithOBBProcess::FindIntersectedGeometricalObjectsWithOBBProcess( ModelPart& rModelPartIntersected, ModelPart& rModelPartIntersecting, const double BoundingBoxFactor, const Flags Options ) : BaseType(rModelPartIntersected, rModelPartIntersecting, Options), mBoundingBoxFactor(BoundingBoxFactor)`

# Calls

- [CreateNewProperties](../../../../functions/sources/model_part/ModelPart/CreateNewProperties.md)
- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [RecursivelyValidateAndAssignDefaults](../../../../functions/sources/kratos_parameters/Parameters/RecursivelyValidateAndAssignDefaults.md)
- [GetBool](../../../../functions/sources/kratos_parameters/Parameters/GetBool.md)
- [GetDouble](../../../../functions/sources/kratos_parameters/Parameters/GetDouble.md)
- [ConvertIntersection](../../../../functions/processes/find_intersected_geometrical_objects_with_obb_process/ConvertIntersection.md)