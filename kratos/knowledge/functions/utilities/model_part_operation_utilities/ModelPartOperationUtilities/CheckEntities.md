---
type: C++ Method
title: CheckEntities
resource: utilities/model_part_operation_utilities.cpp#L260-L291
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

`void ModelPartOperationUtilities::CheckEntities( std::vector<IndexType>& rEntityIdsWithIssues, const TContainerType& rCheckEntities, const std::set<CNodePointersType>& rMainEntities)`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)