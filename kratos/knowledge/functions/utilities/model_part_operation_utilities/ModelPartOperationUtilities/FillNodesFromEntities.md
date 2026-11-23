---
type: C++ Method
title: FillNodesFromEntities
resource: utilities/model_part_operation_utilities.cpp#L85-L96
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

`void ModelPartOperationUtilities::FillNodesFromEntities( std::vector<ModelPart::NodeType*>& rOutput, typename std::vector<TEntityType*>::iterator pEntityBegin, typename std::vector<TEntityType*>::iterator pEntityEnd)`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)