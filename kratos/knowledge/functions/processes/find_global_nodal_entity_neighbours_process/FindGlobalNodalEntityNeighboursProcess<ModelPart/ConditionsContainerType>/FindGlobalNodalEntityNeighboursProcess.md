---
type: C++ Method
title: FindGlobalNodalEntityNeighboursProcess
resource: processes/find_global_nodal_entity_neighbours_process.cpp#L83-L89
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/project/Project/GetModel
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/model_part/FullName
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`FindGlobalNodalEntityNeighboursProcess<ModelPart::ConditionsContainerType>::FindGlobalNodalEntityNeighboursProcess( ModelPart& rModelPart) : mrModel(rModelPart.GetModel()), mrOutputVariable(NEIGHBOUR_CONDITIONS)`

# Calls

- [GetModel](../../../../../functions/python_scripts/project/Project/GetModel.md)
- [FullName](../../../../../functions/includes/model_part/FullName.md)