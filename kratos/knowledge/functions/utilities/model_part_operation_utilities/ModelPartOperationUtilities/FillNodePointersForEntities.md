---
type: C++ Method
title: FillNodePointersForEntities
resource: utilities/model_part_operation_utilities.cpp#L41-L60
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/utilities/model_part_operation_utilities/ModelPartOperationUtilities/CheckValidityOfModelPartsForOperations
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/model_part_operation_utilities/ModelPartOperation
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ModelPartOperationUtilities::FillNodePointersForEntities( std::set<CNodePointersType>& rOutput, const TContainerType& rEntityContainer)`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)

# Called by

- [CheckValidityOfModelPartsForOperations](../../../../functions/utilities/model_part_operation_utilities/ModelPartOperationUtilities/CheckValidityOfModelPartsForOperations.md)
- [ModelPartOperation](../../../../functions/utilities/model_part_operation_utilities/ModelPartOperation.md)