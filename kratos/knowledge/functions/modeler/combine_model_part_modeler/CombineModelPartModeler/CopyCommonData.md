---
type: C++ Method
title: CopyCommonData
resource: modeler/combine_model_part_modeler.cpp#L111-L126
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/model_part/IsSubModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/model_part/pGetNodalSolutionStepVariablesList
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CombineModelPartModeler::CopyCommonData( ModelPart& rCombinedModelPart) const`

# Calls

- [IsSubModelPart](../../../../functions/includes/model_part/IsSubModelPart.md)
- [pGetNodalSolutionStepVariablesList](../../../../functions/includes/model_part/pGetNodalSolutionStepVariablesList.md)