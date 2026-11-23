---
type: C++ Method
title: FindIntersectedGeometricalObjectsProcess
resource: processes/find_intersected_geometrical_objects_process.cpp#L36-L45
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/RecursivelyValidateAndAssignDefaults
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetBool
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`FindIntersectedGeometricalObjectsProcess::FindIntersectedGeometricalObjectsProcess( ModelPart& rModelPartIntersected, ModelPart& rModelPartIntersecting, const Flags Options ) : mrModelPartIntersected(rModelPartIntersected), mrModelPartIntersecting(rModelPartIntersecting), mOptions(Options), mpOctree(new OctreeType())`

# Calls

- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [RecursivelyValidateAndAssignDefaults](../../../../functions/sources/kratos_parameters/Parameters/RecursivelyValidateAndAssignDefaults.md)
- [GetBool](../../../../functions/sources/kratos_parameters/Parameters/GetBool.md)